use crate::error::{MotifError, Result};
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Reads sequences from a FASTA format file and converts them into a Polars DataFrame.
///
/// # Arguments
/// * `filename` - Path to the FASTA file to read
///
/// # Returns
/// * `Result<DataFrame>` - A DataFrame with two columns:
///   - "label": The sequence identifiers (without '>' prefix)
///   - "sequence": The corresponding DNA/RNA sequences in uppercase
///
/// # Errors
/// * Returns `MotifError::InvalidFileFormat` if no sequences are found
/// * Returns `MotifError::DataError` if DataFrame creation fails
/// * Returns `std::io::Error` for file reading issues
pub fn read_fasta(filename: &str) -> Result<DataFrame> {
    let mut sequences: Vec<(String, String)> = Vec::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut current_header = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.starts_with('>') {
            if !current_header.is_empty() {
                sequences.push((current_header, current_sequence.to_uppercase()));
                current_sequence.clear();
            }
            current_header = line[1..].to_string();
        } else if !line.is_empty() {
            current_sequence.push_str(line);
        }
    }

    if !current_header.is_empty() {
        sequences.push((current_header, current_sequence.to_uppercase()));
    }

    if sequences.is_empty() {
        return Err(MotifError::InvalidFileFormat("No sequences found".into()));
    }

    let (labels, sequences): (Vec<String>, Vec<String>) = sequences.into_iter().unzip();
    let df = DataFrame::new(vec![
        Column::new("label".into(), labels),
        Column::new("sequence".into(), sequences),
    ])
    .map_err(|_| MotifError::DataError("Failed to create DataFrame".into()))?;

    Ok(df)
}

/// Writes sequences from a Polars DataFrame to a FASTA format file.
///
/// # Arguments
/// * `df` - DataFrame containing sequences with "label" and "sequence" columns
/// * `filename` - Path where the FASTA file should be written
///
/// # Returns
/// * `Result<()>` - Unit type if successful
///
/// # Errors
/// * Returns `MotifError::DataError` if required columns are missing
/// * Returns `MotifError::Io` for file writing issues
pub fn write_fasta(df: &DataFrame, filename: &str) -> Result<()> {
    let labels = df
        .column("label")
        .map_err(|e| MotifError::DataError(e.to_string()))?
        .str()
        .unwrap();
    let sequences = df
        .column("sequence")
        .map_err(|e| MotifError::DataError(e.to_string()))?
        .str()
        .unwrap();

    let mut file = File::create(filename).map_err(MotifError::Io)?;

    for idx in 0..df.height() {
        let label = labels.get(idx).unwrap();
        let sequence = sequences.get(idx).unwrap();

        writeln!(file, ">{}", label).map_err(MotifError::Io)?;
        writeln!(file, "{}", sequence).map_err(MotifError::Io)?;
    }

    Ok(())
}

/// Generates the reverse complement of a DNA sequence.
///
/// # Arguments
/// * `sequence` - Input DNA sequence string
///
/// # Returns
/// * `String` - The reverse complement sequence where:
///   - A ↔ T
///   - C ↔ G
///
/// # Panics
/// * Panics if the input sequence contains characters other than A, T, C, or G
pub fn rev_comp(sequence: &str) -> String {
    let compliment = HashMap::from([('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]);
    sequence.chars().rev().map(|c| compliment[&c]).collect()
}

/// Calculates the GC content for each sequence in the input DataFrame.
///
/// # Arguments
/// * `df` - DataFrame containing sequences with "label" and "sequence" columns
///
/// # Returns
/// * `Result<DataFrame>` - A DataFrame with:
///   - Original labels
///   - "gc_content": Fraction of G and C bases in each sequence
///
/// # Errors
/// * Returns `MotifError::DataError` if required columns are missing or DataFrame creation fails
pub fn gc_content(df: &DataFrame) -> Result<DataFrame> {
    let sequences = df
        .column("sequence")
        .map_err(|e| MotifError::DataError(e.to_string()))?
        .str()
        .unwrap();

    let gc_content: Vec<f64> = sequences
        .into_iter()
        .map(|seq| {
            let seq = seq.unwrap();
            let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
            gc_count / seq.len() as f64
        })
        .collect();

    let labels = df
        .column("label")
        .map_err(|e| MotifError::DataError(e.to_string()))?;

    let new_df = DataFrame::new(vec![
        labels.clone(),
        Column::new("gc_content".into(), gc_content),
    ])
    .map_err(|e| MotifError::DataError(e.to_string()))?;

    Ok(new_df)
}

/// Identifies sequences containing specified restriction sites.
///
/// # Arguments
/// * `df` - DataFrame containing sequences with "label" and "sequence" columns
/// * `restrictions` - Slice of restriction site patterns to search for
///
/// # Returns
/// * `Result<DataFrame>` - A DataFrame with:
///   - Original labels
///   - "has_restriction_sites": Boolean indicating if any restriction site was found
///
/// # Errors
/// * Returns `MotifError::DataError` if required columns are missing or DataFrame creation fails
pub fn has_restriction_sites(df: &DataFrame, restrictions: &[&str]) -> Result<DataFrame> {
    let restrictions_set: HashSet<String> = restrictions.iter().map(|r| r.to_string()).collect();

    let sequences = df
        .column("sequence")
        .map_err(|e| MotifError::DataError(e.to_string()))?
        .str()
        .unwrap();

    let mask: Vec<bool> = sequences
        .into_iter()
        .map(|seq| {
            let seq = seq.unwrap();
            restrictions_set.iter().any(|r| seq.contains(r))
        })
        .collect();

    let labels = df
        .column("label")
        .map_err(|e| MotifError::DataError(e.to_string()))?;

    let new_df = DataFrame::new(vec![
        labels.clone(),
        Column::new("has_restriction_sites".into(), mask),
    ])
    .map_err(|e| MotifError::DataError(e.to_string()))?;

    Ok(new_df)
}
