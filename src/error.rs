//! Error types for PDBRust library

use std::fmt;
use std::io;
use std::num::{ParseFloatError, ParseIntError};
use thiserror::Error;

/// Errors that can occur when working with PDB files.
#[derive(Debug)]
pub enum PdbError {
    /// An IO error occurred while reading or writing a file.
    IoError(io::Error),
    /// A record in the PDB file was invalid or malformed.
    InvalidRecord(String),
    /// An error occurred while parsing a number.
    ParseError(String),
}

impl fmt::Display for PdbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PdbError::IoError(e) => write!(f, "IO error: {}", e),
            PdbError::InvalidRecord(e) => write!(f, "Invalid record: {}", e),
            PdbError::ParseError(e) => write!(f, "Parse error: {}", e),
        }
    }
}

impl std::error::Error for PdbError {}

impl From<io::Error> for PdbError {
    fn from(err: io::Error) -> Self {
        PdbError::IoError(err)
    }
}

impl From<ParseIntError> for PdbError {
    fn from(err: ParseIntError) -> Self {
        PdbError::ParseError(format!("Failed to parse integer: {}", err))
    }
}

impl From<ParseFloatError> for PdbError {
    fn from(err: ParseFloatError) -> Self {
        PdbError::ParseError(format!("Failed to parse float: {}", err))
    }
}