use clap::Parser;
use polars::prelude::*;
use rayon::prelude::*;
use std::fs;
use std::path::Path;
use tf_binding_rs::occupancy::{read_pwm_to_ewm, total_landscape};
use tf_binding_rs::types::EWMCollection;

#[derive(thiserror::Error, Debug)]
pub enum ScannerError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Polars error: {0}")]
    Polars(#[from] PolarsError),

    #[error("Invalid motif name format")]
    InvalidMotifFormat,

    #[error("Missing sequence column in input file")]
    MissingSequenceColumn,

    #[error("PWM processing error: {0}")]
    PwmError(String),
}

#[derive(Parser)]
#[command(
    name = "motif-scanner",
    about = "Scans DNA sequences for transcription factor binding motifs and calculates predicted occupancy scores",
    long_about = "A tool for analyzing DNA sequences to identify potential transcription factor binding sites. \
                  It processes sequence data, calculates motif positions and their predicted occupancy scores, \
                  and outputs a table of binding sites with their positions, orientations, and occupancy metrics.",
    author = "Jiayu Huang | WUSTL Cohen Lab",
    version,
    after_help = "Example usage:\n    \
                  motif-scanner data.csv motifs.meme results.parquet --cutoff 0.3 --mu 12\n    \
                  motif-scanner sequences.csv pwm.meme output.csv",
    color = clap::ColorChoice::Always
)]
#[derive(Debug)]
struct Args {
    /// Path to input data file (CSV format)
    /// Must contain a 'sequence' column with DNA sequences
    #[arg(value_name = "DATA_FILE")]
    data_file: String,

    /// Path to .meme format file containing Position Weight Matrices (PWMs)
    /// for the motifs to be scanned
    #[arg(value_name = "PWM_FILE")]
    pwm_file: String,

    /// Path for output file (supports .csv or .parquet format)
    /// Will create output directory if it doesn't exist
    #[arg(value_name = "OUTPUT_FILE")]
    output_file: String,

    /// Minimum predicted occupancy threshold
    /// Only motifs with occupancy scores above this value will be included
    /// in the output
    #[arg(long, default_value = "0.2")]
    cutoff: f64,

    /// Predicted affinity parameter (mu) of transcription factor to motif
    /// Higher values indicate stronger binding affinity
    #[arg(long, default_value = "9")]
    mu: i32,
}

trait UnzipN<A, B, C, D, E, F> {
    fn unzip_n_vec(self) -> (Vec<A>, Vec<B>, Vec<C>, Vec<D>, Vec<E>, Vec<F>);
}

impl<I, A, B, C, D, E, F> UnzipN<A, B, C, D, E, F> for I
where
    I: Iterator<Item = (A, B, C, D, E, F)>,
{
    fn unzip_n_vec(self) -> (Vec<A>, Vec<B>, Vec<C>, Vec<D>, Vec<E>, Vec<F>) {
        let mut a_vec = Vec::new();
        let mut b_vec = Vec::new();
        let mut c_vec = Vec::new();
        let mut d_vec = Vec::new();
        let mut e_vec = Vec::new();
        let mut f_vec = Vec::new();

        for (a, b, c, d, e, f) in self {
            a_vec.push(a);
            b_vec.push(b);
            c_vec.push(c);
            d_vec.push(d);
            e_vec.push(e);
            f_vec.push(f);
        }

        (a_vec, b_vec, c_vec, d_vec, e_vec, f_vec)
    }
}

fn process_sequences(
    df: &DataFrame,
    ewm: &EWMCollection,
    mu: f64,
    cutoff: f64,
) -> Result<DataFrame, ScannerError> {
    let sequences = df
        .column("sequence")
        .map_err(|_| ScannerError::MissingSequenceColumn)?;

    let total_seqs = sequences.len();
    println!("{} sequences to scan", total_seqs);

    // convert ChunkedArray<String> to Vec<String> for parallel processing
    let sequences_vec: Vec<_> = sequences.str()?.into_iter().collect();

    // Parallel processing of sequences
    let results: Vec<_> = sequences_vec
        .into_par_iter()
        .enumerate()
        .filter_map(|(idx, seq)| {
            seq.map(|sequence| {
                let landscape = match total_landscape(sequence, ewm, mu) {
                    Ok(l) => l,
                    Err(_) => return Vec::new(),
                };

                let n_positions = landscape.height();
                let mut local_results = Vec::new();

                // Iterate through each position in the landscape
                for pos in 0..n_positions {
                    // Iterate through each motif in the EWM collection
                    for (motif_id, motif_df) in ewm.iter() {
                        // Check both forward and reverse strands
                        for strand in ["F", "R"] {
                            let col_name = format!("{}_{}", motif_id, strand);

                            // Get the column for this motif+strand from the landscape
                            if let Ok(motif_col) = landscape.column(&col_name) {
                                if let Ok(occ) = motif_col.get(pos).unwrap().try_extract::<f64>() {
                                    if occ > cutoff {
                                        local_results.push((
                                            idx as i32,
                                            pos as i32,
                                            motif_id.split('_').next().unwrap().to_string(),
                                            strand.to_string(),
                                            motif_df.height() as i32,
                                            occ,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
                local_results
            })
        })
        .flatten()
        .collect();

    // Unzip results into separate vectors
    let (labels, positions, motifs, strands, lengths, occupancies): (
        Vec<i32>,
        Vec<i32>,
        Vec<String>,
        Vec<String>,
        Vec<i32>,
        Vec<f64>,
    ) = results.into_iter().unzip_n_vec();

    let df = DataFrame::new(vec![
        Column::new("label".into(), labels),
        Column::new("position".into(), positions),
        Column::new("motif".into(), motifs),
        Column::new("strand".into(), strands),
        Column::new("length".into(), lengths),
        Column::new("occupancy".into(), occupancies),
    ])?;

    Ok(df)
}

fn save_results(df: &mut DataFrame, output_file: &str) -> Result<(), ScannerError> {
    match Path::new(output_file)
        .extension()
        .and_then(|ext| ext.to_str())
    {
        Some("parquet") => {
            let mut file = std::fs::File::create(output_file)?;
            ParquetWriter::new(&mut file)
                .with_compression(ParquetCompression::Snappy)
                .finish(df)?;
        }
        _ => {
            let mut file = std::fs::File::create(output_file)?;
            CsvWriter::new(&mut file).include_header(true).finish(df)?;
        }
    }

    Ok(())
}

fn main() -> Result<(), ScannerError> {
    let start_time = std::time::Instant::now();

    let args = Args::parse();

    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(&args.output_file).parent() {
        fs::create_dir_all(parent)?;
    }

    let df = LazyCsvReader::new(&args.data_file)
        .with_has_header(true)
        .finish()?
        .filter(col("sequence").str().contains(lit("N"), false).not())
        .filter(col("sequence").str().contains(lit("Y"), false).not())
        .collect()?;

    // read pwm file and convert to ewm
    let ewm = read_pwm_to_ewm(&args.pwm_file).map_err(|e| ScannerError::PwmError(e.to_string()))?;

    let mut results_df = process_sequences(&df, &ewm, args.mu as f64, args.cutoff)?;

    let elapsed = start_time.elapsed();
    println!(
        "Total execution time: {:.4} minutes",
        elapsed.as_secs_f64() / 60.0
    );

    // save results
    save_results(&mut results_df, &args.output_file)?;

    Ok(())
}
