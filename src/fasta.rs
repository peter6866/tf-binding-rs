use crate::error::{MotifError, Result};
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Read a FASTA file and return sequences as a Polars Dataframe
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

/// Write sequences from a Polars DataFrame to a FASTA file
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

/// Take the reverse compliment of a sequence
pub fn rev_comp(sequence: &str) -> String {
    let compliment = HashMap::from([('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]);
    sequence.chars().rev().map(|c| compliment[&c]).collect()
}

/// Calculate the GC content of every sequence in a Polars DataFrame
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

/// Generate a boolean mask indicating which sequences contain restriction sites
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
