use crate::error::MotifError;
use crate::types::*;
use polars::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::iter::Peekable;

/// Advances the iterator until a MOTIF line is found
fn skip_until_motif<I>(lines: &mut Peekable<I>)
where
    I: Iterator<Item = Result<String, std::io::Error>>,
{
    while let Some(Ok(line)) = lines.peek() {
        if line.starts_with("MOTIF") {
            break;
        }
        lines.next();
    }
}

/// Parses a single PWM from the iterator
fn parse_pwm<I>(lines: &mut I) -> Result<Option<(String, PWM)>, MotifError>
where
    I: Iterator<Item = Result<String, std::io::Error>>,
{
    // Get motif ID from MOTIF line
    let motif_line = match lines.next() {
        Some(Ok(line)) if line.starts_with("MOTIF") => line,
        _ => return Ok(None),
    };

    let motif_id = motif_line
        .split_whitespace()
        .nth(1)
        .ok_or_else(|| MotifError::InvalidFileFormat("Missing motif ID".into()))?
        .to_string();

    // Skip header lines
    for _ in 0..2 {
        lines.next();
    }

    // Read PWM rows until we hit a non-PWM line
    let pwm_rows: Vec<Vec<f64>> = lines
        .take_while(|line| {
            line.as_ref()
                .map(|l| l.starts_with(|c: char| c.is_whitespace() || c == '0' || c == '1'))
                .unwrap_or(false)
        })
        .map(|line| {
            let line = line.map_err(|e| MotifError::Io(e))?;
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| MotifError::InvalidFileFormat(format!("Invalid PWM value: {}", e)))?;

            Ok(values)
        })
        .collect::<Result<Vec<_>, MotifError>>()?;

    if pwm_rows.is_empty() {
        return Err(MotifError::InvalidFileFormat("Empty PWM".into()));
    }

    // Create PWM DataFrame
    let pwm = DataFrame::new(vec![
        Column::new(
            "A".into(),
            pwm_rows.iter().map(|row| row[0]).collect::<Vec<f64>>(),
        ),
        Column::new(
            "C".into(),
            pwm_rows.iter().map(|row| row[1]).collect::<Vec<f64>>(),
        ),
        Column::new(
            "G".into(),
            pwm_rows.iter().map(|row| row[2]).collect::<Vec<f64>>(),
        ),
        Column::new(
            "T".into(),
            pwm_rows.iter().map(|row| row[3]).collect::<Vec<f64>>(),
        ),
    ])
    .map_err(|e| MotifError::DataError(e.to_string()))?;

    Ok(Some((motif_id, pwm)))
}

/// Reads Position Weight Matrices (PWMs) from a MEME format file
///
/// This function parses a MEME-formatted file containing one or more Position Weight Matrices,
/// each identified by a unique motif ID. The PWMs represent DNA binding motifs where each position
/// contains probabilities for the four nucleotides (A, C, G, T).
///
/// # Arguments
/// * `filename` - Path to the MEME format file to read
///
/// # Returns
/// * `Result<PWMCollection, MotifError>` - A HashMap where keys are motif IDs and values are their corresponding PWMs
///
/// # Errors
/// * `MotifError::Io` - If the file cannot be opened or read
/// * `MotifError::InvalidFileFormat` - If the file format is invalid or no PWMs are found
/// * `MotifError::DataError` - If there are issues creating the PWM DataFrame
///
/// # Example
/// ```ignore
/// use tf_binding_rs::occupancy::read_pwm_files;
///
/// let pwms = read_pwm_files("path/to/motifs.meme").unwrap();
/// for (motif_id, pwm) in pwms {
///     println!("Found motif: {}", motif_id);
/// }
/// ```
///
/// # Format
/// The input file should be in MEME format, where each PWM is preceded by a "MOTIF" line
/// containing the motif ID, followed by the matrix values.
pub fn read_pwm_files(filename: &str) -> Result<PWMCollection, MotifError> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines().peekable();
    let mut pwms = HashMap::new();

    // Skip header until first MOTIF
    skip_until_motif(&mut lines);

    // Parse all PWMs
    while let Some((id, pwm)) = parse_pwm(&mut lines)? {
        pwms.insert(id, pwm);
        skip_until_motif(&mut lines);
    }

    if pwms.is_empty() {
        return Err(MotifError::InvalidFileFormat("No PWMs found".into()));
    }

    Ok(pwms)
}
