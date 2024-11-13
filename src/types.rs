use polars::prelude::*;
use std::collections::HashMap;

/// Represents a Position Weight Matrix (PWM)
/// Stored as a DataFrame with columns A, C, G, T
pub type PWM = DataFrame;

/// Collection of PWMs indexed by motif ID
pub type PWMCollection = HashMap<String, PWM>;
