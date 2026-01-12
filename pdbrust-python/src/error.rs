//! Error handling and Python exception mapping

use pdbrust::PdbError;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::PyErr;

/// Convert a PdbError to a Python exception
pub fn convert_error(err: PdbError) -> PyErr {
    match err {
        PdbError::IoError(e) => PyIOError::new_err(e.to_string()),
        PdbError::InvalidRecord(msg) => PyValueError::new_err(msg),
        PdbError::ParseError(msg) => PyValueError::new_err(format!("Parse error: {}", msg)),
        PdbError::AtomCountMismatch { expected, found } => PyValueError::new_err(format!(
            "Atom count mismatch: expected {} atoms, found {}",
            expected, found
        )),
        PdbError::NoAtomsSelected(msg) => {
            PyValueError::new_err(format!("No atoms selected: {}", msg))
        }
        PdbError::InsufficientAtoms(msg) => {
            PyValueError::new_err(format!("Insufficient atoms: {}", msg))
        }
    }
}

/// Convert RCSB download errors to Python exceptions
#[cfg(feature = "rcsb")]
pub fn convert_download_error(err: pdbrust::rcsb::DownloadError) -> PyErr {
    use pdbrust::rcsb::DownloadError;
    match err {
        DownloadError::NotFound(id) => PyValueError::new_err(format!("PDB ID not found: {}", id)),
        DownloadError::RequestFailed(msg) => PyRuntimeError::new_err(msg),
        DownloadError::ParseError(e) => convert_error(e),
        DownloadError::IoError(e) => PyIOError::new_err(e.to_string()),
    }
}

/// Convert RCSB search errors to Python exceptions
#[cfg(feature = "rcsb")]
pub fn convert_search_error(err: pdbrust::rcsb::SearchError) -> PyErr {
    use pdbrust::rcsb::SearchError;
    match err {
        SearchError::RequestFailed(msg) => PyRuntimeError::new_err(msg),
        SearchError::ParseError(msg) => PyValueError::new_err(format!("Parse error: {}", msg)),
        SearchError::ApiError(msg) => PyRuntimeError::new_err(format!("API error: {}", msg)),
        SearchError::InvalidQuery(msg) => PyValueError::new_err(format!("Invalid query: {}", msg)),
    }
}
