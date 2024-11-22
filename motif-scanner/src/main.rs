use clap::Parser;
use polars::prelude::*;
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

fn process_sequences(
    df: &DataFrame,
    ewm: &EWMCollection,
    mu: f64,
    cutoff: f64,
) -> Result<DataFrame, ScannerError> {
    let sequences = df
        .column("sequence")
        .map_err(|_| ScannerError::MissingSequenceColumn)?;

    let mut labels: Vec<i32> = Vec::new();
    let mut positions: Vec<i32> = Vec::new();
    let mut motifs: Vec<String> = Vec::new();
    let mut strands: Vec<String> = Vec::new();
    let mut lengths: Vec<i32> = Vec::new();
    let mut occupancies: Vec<f64> = Vec::new();

    let total_seqs = sequences.len();
    println!("{} sequences to scan", total_seqs);

    for (idx, seq) in sequences.str()?.into_iter().enumerate() {
        if let Some(sequence) = seq {
            let landscape = total_landscape(sequence, ewm, mu).map_err(|_| {
                ScannerError::PwmError(format!(
                    "Error calculating occupancy landscape for sequence {}",
                    idx
                ))
            })?;

            // Get the height (number of positions) of the landscape DataFrame
            let n_positions = landscape.height();

            // Iterate through each motif in the EWM collection
            for (motif_id, motif_df) in ewm.iter() {
                // Check both forward and reverse strands
                for strand in ["F", "R"] {
                    let col_name = format!("{}_{}", motif_id, strand);

                    // Get the column for this motif+strand from the landscape
                    if let Ok(motif_col) = landscape.column(&col_name) {
                        // Iterate through positions
                        for pos in 0..n_positions {
                            if let Ok(occ) = motif_col.get(pos).unwrap().try_extract::<f64>() {
                                if occ > cutoff {
                                    labels.push(idx as i32);
                                    positions.push(pos as i32);
                                    motifs.push(motif_id.split('_').next().unwrap().to_string());
                                    strands.push(strand.to_string());
                                    lengths.push(motif_df.height() as i32);
                                    occupancies.push(occ);
                                }
                            }
                        }
                    }
                }
            }
        }

        if (idx + 1) % 5000 == 0 {
            println!("\t{} / {} sequences scanned", idx + 1, total_seqs);
        }
    }

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

    let results_df = process_sequences(&df, &ewm, args.mu as f64, args.cutoff)?;

    let elapsed = start_time.elapsed();
    println!(
        "Total execution time: {:.4} minutes",
        elapsed.as_secs_f64() / 60.0
    );

    Ok(())
}
