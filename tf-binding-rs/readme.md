# tf-binding-rs (In Development)

[![github](https://img.shields.io/badge/github-peter6866/tf--binding--rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github)](https://github.com/peter6866/tf-binding-rs)
[![crates.io](https://img.shields.io/crates/v/tf-binding-rs.svg?style=for-the-badge&color=fc8d62&logo=rust)](https://crates.io/crates/tf-binding-rs)

A Rust library for predicting transcription factor (TF) binding site occupancy in DNA sequences. This toolkit provides efficient implementations for:

- FASTA file manipulation and sequence processing
- Position Weight Matrix (PWM) handling and Energy Weight Matrix (EWM) conversion
- TF binding site occupancy prediction using statistical thermodynamics
- Binding energy landscape and occupancy probability calculations
- Multi-TF occupancy analysis

## Features

- 🧬 Fast FASTA file reading and writing
- 📊 PWM/EWM-based binding site analysis
- 🔍 Efficient sequence scanning with energy matrices
- 📈 Occupancy landscape calculation for multiple TFs
- 🧮 Statistical thermodynamics-based predictions

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

### Working with PWMs and Energy Matrices

```rust
use tf_binding_rs::occupancy;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read PWMs and convert to Energy Weight Matrices
    let ewm_collection = occupancy::read_pwm_to_ewm("path/to/motifs.meme")?;

    // Calculate binding landscape for a sequence
    let sequence = "ATCGATCGTAGCTACGT";
    let mu = -3.0; // chemical potential

    // Get occupancy predictions for all TFs
    let occupancy_landscape = occupancy::total_landscape(
        &sequence,
        &ewm_collection,
        mu
    )?;

    println!("Occupancy predictions:\n{}", occupancy_landscape);
    Ok(())
}
```

## Use Cases

- Genomic sequence analysis
- TF binding site prediction and quantification
- Multi-factor binding landscape analysis
- Regulatory sequence characterization
- Statistical thermodynamics of protein-DNA interactions

## Documentation

For detailed API documentation, visit [docs.rs/tf-binding-rs](https://docs.rs/tf-binding-rs)
