use std::io;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MotifError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    #[error("Invalid sequence at position {position}: {message}")]
    InvalidSequence { position: usize, message: String },

    #[error("Invalid PWM format: {0}")]
    InvalidPwm(String),

    #[error("Invalid file format: {0}")]
    InvalidFileFormat(String),

    #[error("Data error: {0}")]
    DataError(String),

    #[error("Invalid parameter: {name} = {value}, {message}")]
    InvalidParameter {
        name: String,
        value: String,
        message: String,
    },

    #[error("Invalid input: {0}")]
    InvalidInput(String),
}

// // Type alias for Result with MotifError
// pub type Result<T> = std::result::Result<T, MotifError>;

impl MotifError {
    /// Create a new InvalidSequence error
    pub fn invalid_sequence(position: usize, message: impl Into<String>) -> Self {
        MotifError::InvalidSequence {
            position,
            message: message.into(),
        }
    }

    /// Create a new InvalidPwm error
    pub fn invalid_pwm(message: impl Into<String>) -> Self {
        MotifError::InvalidPwm(message.into())
    }

    /// Create a new InvalidParameter error
    pub fn invalid_parameter(
        name: impl Into<String>,
        value: impl ToString,
        message: impl Into<String>,
    ) -> Self {
        MotifError::InvalidParameter {
            name: name.into(),
            value: value.to_string(),
            message: message.into(),
        }
    }
}
