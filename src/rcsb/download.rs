//! RCSB PDB file download utilities.
//!
//! This module provides functionality for downloading structure files
//! from the RCSB Protein Data Bank.

use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;

use crate::core::PdbStructure;
use crate::parser::{parse_mmcif_string, parse_pdb_string};
use crate::PdbError;

use super::RCSB_DOWNLOAD_URL;

/// Supported file formats for download.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    /// PDB format (.pdb)
    Pdb,
    /// mmCIF format (.cif)
    Cif,
}

impl FileFormat {
    /// Get the file extension for this format.
    pub fn extension(&self) -> &'static str {
        match self {
            FileFormat::Pdb => "pdb",
            FileFormat::Cif => "cif",
        }
    }

    /// Get the compressed file extension for this format.
    pub fn compressed_extension(&self) -> &'static str {
        match self {
            FileFormat::Pdb => "pdb.gz",
            FileFormat::Cif => "cif.gz",
        }
    }
}

impl fmt::Display for FileFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FileFormat::Pdb => write!(f, "PDB"),
            FileFormat::Cif => write!(f, "mmCIF"),
        }
    }
}

/// Errors that can occur during download operations.
#[derive(Debug)]
pub enum DownloadError {
    /// HTTP request failed
    RequestFailed(String),
    /// PDB ID not found
    NotFound(String),
    /// Failed to parse the downloaded file
    ParseError(PdbError),
    /// I/O error
    IoError(std::io::Error),
}

impl fmt::Display for DownloadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DownloadError::RequestFailed(msg) => write!(f, "Download failed: {}", msg),
            DownloadError::NotFound(pdb_id) => write!(f, "PDB ID not found: {}", pdb_id),
            DownloadError::ParseError(err) => write!(f, "Parse error: {}", err),
            DownloadError::IoError(err) => write!(f, "I/O error: {}", err),
        }
    }
}

impl std::error::Error for DownloadError {}

impl From<reqwest::Error> for DownloadError {
    fn from(err: reqwest::Error) -> Self {
        DownloadError::RequestFailed(err.to_string())
    }
}

impl From<PdbError> for DownloadError {
    fn from(err: PdbError) -> Self {
        DownloadError::ParseError(err)
    }
}

impl From<std::io::Error> for DownloadError {
    fn from(err: std::io::Error) -> Self {
        DownloadError::IoError(err)
    }
}

/// Build the download URL for a PDB structure.
fn build_download_url(pdb_id: &str, format: FileFormat) -> String {
    let pdb_id_upper = pdb_id.to_uppercase();
    match format {
        FileFormat::Pdb => format!("{}/{}.pdb", RCSB_DOWNLOAD_URL, pdb_id_upper),
        FileFormat::Cif => format!("{}/{}.cif", RCSB_DOWNLOAD_URL, pdb_id_upper),
    }
}

/// Download a PDB file as a string.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `format` - The desired file format
///
/// # Returns
///
/// The file contents as a string.
///
/// # Errors
///
/// Returns a `DownloadError` if the download fails.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_pdb_string, FileFormat};
///
/// let content = download_pdb_string("1UBQ", FileFormat::Pdb)?;
/// println!("Downloaded {} bytes", content.len());
/// ```
pub fn download_pdb_string(pdb_id: &str, format: FileFormat) -> Result<String, DownloadError> {
    let url = build_download_url(pdb_id, format);
    let client = reqwest::blocking::Client::new();

    let response = client.get(&url).send()?;

    if response.status() == reqwest::StatusCode::NOT_FOUND {
        return Err(DownloadError::NotFound(pdb_id.to_string()));
    }

    if !response.status().is_success() {
        return Err(DownloadError::RequestFailed(format!(
            "HTTP {}: {}",
            response.status(),
            response.text().unwrap_or_default()
        )));
    }

    Ok(response.text()?)
}

/// Download and parse a structure from RCSB PDB.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `format` - The desired file format
///
/// # Returns
///
/// A parsed `PdbStructure`.
///
/// # Errors
///
/// Returns a `DownloadError` if the download or parsing fails.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_structure, FileFormat};
///
/// let structure = download_structure("1UBQ", FileFormat::Pdb)?;
/// println!("Downloaded {} atoms", structure.atoms.len());
/// ```
pub fn download_structure(
    pdb_id: &str,
    format: FileFormat,
) -> Result<PdbStructure, DownloadError> {
    let content = download_pdb_string(pdb_id, format)?;

    let structure = match format {
        FileFormat::Pdb => parse_pdb_string(&content)?,
        FileFormat::Cif => parse_mmcif_string(&content)?,
    };

    Ok(structure)
}

/// Download a structure to a file.
///
/// # Arguments
///
/// * `pdb_id` - The 4-character PDB ID
/// * `path` - The path to save the file
/// * `format` - The desired file format
///
/// # Errors
///
/// Returns a `DownloadError` if the download or file writing fails.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_to_file, FileFormat};
///
/// download_to_file("1UBQ", "1UBQ.pdb", FileFormat::Pdb)?;
/// ```
pub fn download_to_file<P: AsRef<Path>>(
    pdb_id: &str,
    path: P,
    format: FileFormat,
) -> Result<(), DownloadError> {
    let content = download_pdb_string(pdb_id, format)?;
    let mut file = File::create(path)?;
    file.write_all(content.as_bytes())?;
    Ok(())
}

/// Download multiple structures.
///
/// # Arguments
///
/// * `pdb_ids` - A slice of PDB IDs to download
/// * `format` - The desired file format
///
/// # Returns
///
/// A vector of results, one for each PDB ID.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_multiple, FileFormat};
///
/// let results = download_multiple(&["1UBQ", "8HM2"], FileFormat::Pdb);
/// for (pdb_id, result) in results {
///     match result {
///         Ok(structure) => println!("{}: {} atoms", pdb_id, structure.atoms.len()),
///         Err(e) => eprintln!("{}: {}", pdb_id, e),
///     }
/// }
/// ```
pub fn download_multiple(
    pdb_ids: &[&str],
    format: FileFormat,
) -> Vec<(String, Result<PdbStructure, DownloadError>)> {
    pdb_ids
        .iter()
        .map(|&pdb_id| {
            let result = download_structure(pdb_id, format);
            (pdb_id.to_string(), result)
        })
        .collect()
}

/// Download multiple structures to files.
///
/// # Arguments
///
/// * `pdb_ids` - A slice of PDB IDs to download
/// * `output_dir` - The directory to save files in
/// * `format` - The desired file format
///
/// # Returns
///
/// A vector of results indicating success or failure for each download.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::rcsb::{download_multiple_to_files, FileFormat};
///
/// let results = download_multiple_to_files(&["1UBQ", "8HM2"], "./structures", FileFormat::Pdb);
/// ```
pub fn download_multiple_to_files<P: AsRef<Path>>(
    pdb_ids: &[&str],
    output_dir: P,
    format: FileFormat,
) -> Vec<(String, Result<(), DownloadError>)> {
    let dir = output_dir.as_ref();

    pdb_ids
        .iter()
        .map(|&pdb_id| {
            let filename = format!("{}.{}", pdb_id.to_uppercase(), format.extension());
            let path = dir.join(filename);
            let result = download_to_file(pdb_id, path, format);
            (pdb_id.to_string(), result)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_format_extension() {
        assert_eq!(FileFormat::Pdb.extension(), "pdb");
        assert_eq!(FileFormat::Cif.extension(), "cif");
    }

    #[test]
    fn test_file_format_compressed_extension() {
        assert_eq!(FileFormat::Pdb.compressed_extension(), "pdb.gz");
        assert_eq!(FileFormat::Cif.compressed_extension(), "cif.gz");
    }

    #[test]
    fn test_file_format_display() {
        assert_eq!(FileFormat::Pdb.to_string(), "PDB");
        assert_eq!(FileFormat::Cif.to_string(), "mmCIF");
    }

    #[test]
    fn test_build_download_url_pdb() {
        let url = build_download_url("1ubq", FileFormat::Pdb);
        assert_eq!(url, "https://files.rcsb.org/download/1UBQ.pdb");
    }

    #[test]
    fn test_build_download_url_cif() {
        let url = build_download_url("8hm2", FileFormat::Cif);
        assert_eq!(url, "https://files.rcsb.org/download/8HM2.cif");
    }

    #[test]
    fn test_download_error_display() {
        let err = DownloadError::NotFound("XXXX".to_string());
        assert_eq!(err.to_string(), "PDB ID not found: XXXX");

        let err = DownloadError::RequestFailed("timeout".to_string());
        assert_eq!(err.to_string(), "Download failed: timeout");
    }
}
