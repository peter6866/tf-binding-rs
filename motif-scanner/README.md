# motif-scanner

[![github](https://img.shields.io/badge/github-peter6866/tf--binding--rs-8da0cb?style=for-the-badge&labelColor=555555&logo=github)](https://github.com/peter6866/tf-binding-rs)
[![crates.io](https://img.shields.io/crates/v/motif-scanner.svg?style=for-the-badge&color=fc8d62&logo=rust)](https://crates.io/crates/motif-scanner)

A command-line tool for scanning DNA sequences and predicting transcription factor binding sites.

## Features

- 🧬 Batch processing of sequence files
- 📊 PWM/EWM-based binding site analysis
- 🔍 Configurable occupancy threshold filtering
- 📈 Multiple output formats (CSV, Parquet)
- ⚡ Parallel processing for large datasets

## Installation

### From crates.io

```bash
cargo install motif-scanner
```

### From Source

```bash
git clone https://github.com/peter6866/tf-binding-rs
cd tf-binding-rs
cargo install --path motif-scanner
```

## Usage

Basic usage:

```bash
motif-scanner input.csv motifs.meme output.csv
```

With options:

```bash
motif-scanner input.csv motifs.meme output.parquet --cutoff 0.3 --mu 12
```

### Arguments

- `DATA_FILE`: Input CSV file containing sequences (must have a 'sequence' column)
- `PWM_FILE`: MEME format file containing Position Weight Matrices
- `OUTPUT_FILE`: Path for output file (.csv or .parquet format)
- `--cutoff`: Minimum occupancy threshold (default: 0.2)
- `--mu`: Chemical potential parameter (default: 9)

### Input Format

The input CSV file must contain a column named 'sequence' with DNA sequences:

```csv
id,sequence
seq1,ATCGATCGTGCTAGCTA
seq2,GCTAGCTAGCTAGCTAG
```

### Output Format

The tool generates a table with the following columns:

- `label`: Sequence index from input file
- `position`: Position of the binding site
- `motif`: Name of the transcription factor
- `strand`: Binding strand (F/R)
- `length`: Length of the motif
- `occupancy`: Predicted occupancy score

## Example

```bash
# Scan sequences with default parameters
motif-scanner sequences.csv pwm.meme results.csv

# Use stricter threshold and higher chemical potential
motif-scanner sequences.csv pwm.meme results.parquet --cutoff 0.4 --mu 15

# Process and save as Parquet format
motif-scanner data.csv motifs.meme output.parquet
```

## Performance

The tool uses parallel processing for efficient scanning of large sequence datasets. Memory usage scales with the number of input sequences and motifs being scanned.
