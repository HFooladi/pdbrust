//! Error types for PDBRust library

use std::fmt;
use std::io;
use std::num::{ParseFloatError, ParseIntError};

/// Errors that can occur when working with PDB files.
#[derive(Debug)]
pub enum PdbError {
    /// An IO error occurred while reading or writing a file.
    IoError(io::Error),
    /// A record in the PDB file was invalid or malformed.
    InvalidRecord(String),
    /// An error occurred while parsing a number.
    ParseError(String),
    /// Structures have mismatched atom counts for alignment/RMSD.
    AtomCountMismatch {
        /// Expected number of atoms (from target structure)
        expected: usize,
        /// Found number of atoms (from mobile structure)
        found: usize,
    },
    /// No atoms found for the specified selection criteria.
    NoAtomsSelected(String),
    /// Insufficient atoms for the operation (need at least 3 for alignment).
    InsufficientAtoms(String),
}

impl fmt::Display for PdbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PdbError::IoError(e) => write!(f, "IO error: {}", e),
            PdbError::InvalidRecord(e) => write!(f, "Invalid record: {}", e),
            PdbError::ParseError(e) => write!(f, "Parse error: {}", e),
            PdbError::AtomCountMismatch { expected, found } => {
                write!(
                    f,
                    "Atom count mismatch: expected {} atoms, found {}",
                    expected, found
                )
            }
            PdbError::NoAtomsSelected(msg) => write!(f, "No atoms selected: {}", msg),
            PdbError::InsufficientAtoms(msg) => write!(f, "Insufficient atoms: {}", msg),
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
