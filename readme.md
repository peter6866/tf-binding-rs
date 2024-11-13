# tf-binding-rs (In Development)

[<img alt="github" src="https://img.shields.io/badge/github-peter6866/tf--binding--rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/peter6866/tf-binding-rs)
[<img alt="crates.io" src="https://img.shields.io/crates/v/tf-binding-rs.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/tf-binding-rs)

A Rust library for predicting transcription factor (TF) binding site occupancy in DNA sequences. This toolkit provides efficient implementations for:

- FASTA file manipulation and sequence processing
- Position Weight Matrix (PWM) handling and Energy Weight Matrix (EWM) conversion
- TF binding site occupancy prediction using biophysical models
- Information content and binding site diversity calculations

Built with performance in mind, this library offers a fast and memory-efficient alternative to traditional Python implementations for genomic analysis.

## Features

- 🧬 Fast FASTA file reading and writing
- 📊 PWM/EWM-based binding site scoring
- 🔍 Efficient sequence scanning for binding sites
- 📈 Occupancy landscape calculation
- 🧮 Statistical and thermodynamic calculations

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
tf-binding-rs = "0.1.1"
```

Or install using cargo:

```bash
cargo add tf-binding-rs
```

## Examples

### Reading FASTA Files

```rust
use tf_binding_rs::fasta;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read sequences from a FASTA file
    let sequences = fasta::read_fasta("path/to/sequences.fasta")?;

    // Print sequence information
    println!("Number of sequences: {}", sequences.height());

    // Calculate GC content
    let gc_stats = fasta::gc_content(&sequences)?;
    println!("GC content analysis: {:?}", gc_stats);

    Ok(())
}
```

### Working with PWM Files

```rust
use tf_binding_rs::occupancy;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read PWM motifs from MEME format file
    let pwm_collection = occupancy::read_pwm_files("path/to/motifs.meme")?;

    // Process each motif
    for (motif_id, pwm) in pwm_collection {
        println!("Processing motif: {}", motif_id);
        println!("Matrix dimensions: {:?}", pwm.shape());
    }

    Ok(())
}
```

## Use Cases

- Genomic sequence analysis
- TF binding site prediction
- Regulatory sequence characterization
- High-throughput DNA sequence processing

## Documentation

For detailed API documentation, visit [docs.rs/tf-binding-rs](https://docs.rs/tf-binding-rs)
