//! RCSB PDB Search and Download functionality.
//!
//! This module provides utilities for searching and downloading structures
//! from the RCSB Protein Data Bank using their REST API.
//!
//! # Examples
//!
//! ## Download a structure by PDB ID
//!
//! ```ignore
//! use pdbrust::rcsb::{download_structure, FileFormat};
//!
//! // Download and parse a structure
//! let structure = download_structure("1UBQ", FileFormat::Pdb)?;
//! println!("Downloaded {} atoms", structure.atoms.len());
//! ```
//!
//! ## Search for structures
//!
//! ```ignore
//! use pdbrust::rcsb::{SearchQuery, rcsb_search};
//!
//! // Search for human insulin structures
//! let query = SearchQuery::new()
//!     .with_text("insulin")
//!     .with_organism("Homo sapiens")
//!     .with_resolution_max(2.0);
//!
//! let pdb_ids = rcsb_search(&query, 10)?;
//! println!("Found {} structures", pdb_ids.len());
//! ```
//!
//! ## Download to file
//!
//! ```ignore
//! use pdbrust::rcsb::{download_to_file, FileFormat};
//!
//! download_to_file("1UBQ", "1UBQ.pdb", FileFormat::Pdb)?;
//! ```
//!
//! # Feature Flag
//!
//! This module requires the `rcsb` feature:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.6", features = ["rcsb"] }
//! ```

mod download;
mod search;

#[cfg(feature = "rcsb-async")]
mod async_download;

pub use download::{
    DownloadError, FileFormat, download_multiple, download_multiple_to_files, download_pdb_string,
    download_structure, download_to_file,
};
pub use search::{
    ExperimentalMethod, PolymerType, SearchError, SearchQuery, SearchResult, rcsb_search,
};

#[cfg(feature = "rcsb-async")]
pub use async_download::{
    AsyncDownloadOptions, download_multiple_async, download_multiple_to_files_async,
    download_pdb_string_async, download_structure_async, download_to_file_async,
};

/// RCSB PDB base URLs
pub const RCSB_DOWNLOAD_URL: &str = "https://files.rcsb.org/download";
pub const RCSB_SEARCH_URL: &str = "https://search.rcsb.org/rcsbsearch/v2/query";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_format_extension() {
        assert_eq!(FileFormat::Pdb.extension(), "pdb");
        assert_eq!(FileFormat::Cif.extension(), "cif");
        assert_eq!(FileFormat::Pdb.compressed_extension(), "pdb.gz");
        assert_eq!(FileFormat::Cif.compressed_extension(), "cif.gz");
    }

    #[test]
    fn test_search_query_builder() {
        let query = SearchQuery::new()
            .with_text("kinase")
            .with_resolution_max(2.5);

        assert!(query.text.is_some());
        assert_eq!(query.text.as_deref(), Some("kinase"));
        assert_eq!(query.resolution_max, Some(2.5));
    }

    #[test]
    fn test_search_query_to_json() {
        let query = SearchQuery::new().with_text("ubiquitin");

        let json = query.to_json();
        assert!(json.contains("ubiquitin"));
    }
}
