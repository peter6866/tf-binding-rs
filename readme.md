# TF Binding Analysis Tools

This workspace contains tools for analyzing transcription factor (TF) binding sites in DNA sequences:

- **[tf-binding-rs](tf-binding-rs/)**: A Rust library for TF binding site prediction and sequence analysis
- **[motif-scanner](motif-scanner/)**: A command-line tool for scanning DNA sequences for TF binding sites

## 🧬 tf-binding-rs

[<img alt="github" src="https://img.shields.io/badge/github-peter6866/tf--binding--rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/peter6866/tf-binding-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/tf-binding-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/tf-binding-rs)

A Rust library providing efficient implementations for:

- FASTA file manipulation and sequence processing
- Position Weight Matrix (PWM) handling
- Energy Weight Matrix (EWM) conversion
- TF binding site occupancy prediction
- Multi-TF occupancy analysis

```rust
use tf_binding_rs::occupancy;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ewms = occupancy::read_pwm_to_ewm("motifs.meme")?;
    let sequence = "ATCGATCGTAGCTACGT";
    let landscape = occupancy::total_landscape(sequence, &ewms, -3.0)?;
    println!("Binding landscape:\n{}", landscape);
    Ok(())
}
```

## 🔍 motif-scanner

A command-line tool for scanning DNA sequences and predicting TF binding sites. Features:

- Batch processing of sequence files
- Occupancy score calculation
- Multiple output formats (CSV, Parquet)
- Filtering by occupancy threshold

```bash
motif-scanner input.csv motifs.meme output.csv --cutoff 0.2 --mu 9.0
```

## Installation

### From Source

```bash
# Clone the repository
git clone https://github.com/peter6866/tf-binding-rs
cd tf-binding-rs

# Build both the library and the scanner
cargo build --release --workspace

# Install the motif-scanner binary
cargo install --path motif-scanner
```

### From crates.io

```bash
# Install just the motif-scanner tool
cargo install motif-scanner

# For library usage, add to your Cargo.toml:
[dependencies]
tf-binding-rs = "0.1.4"
```

## Documentation

- [tf-binding-rs API Documentation](https://docs.rs/tf-binding-rs)
- [motif-scanner Usage Guide](motif-scanner/README.md)
