use clap::Parser;

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

fn main() {
    let args = Args::parse();
    println!("{:?}", args);
}
