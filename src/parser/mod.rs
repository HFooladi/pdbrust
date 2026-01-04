//! Parser module for different molecular structure file formats
//!
//! This module provides parsing functionality for both PDB and mmCIF formats,
//! with automatic format detection capabilities.
//!
//! # Supported Formats
//!
//! - **PDB Format**: Legacy text-based format with fixed-width columns
//! - **mmCIF Format**: Modern dictionary-based format used by PDB
//! - **Gzip-compressed files**: With the `gzip` feature (`.ent.gz`, `.pdb.gz`, `.cif.gz`)
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::parser::{parse_pdb_file, parse_mmcif_file, parse_structure_file};
//!
//! // Parse specific formats
//! let pdb_structure = parse_pdb_file("structure.pdb")?;
//! let mmcif_structure = parse_mmcif_file("structure.cif")?;
//!
//! // Auto-detect format
//! let structure = parse_structure_file("structure.ent")?;
//!
//! // Parse gzip-compressed files (requires "gzip" feature)
//! #[cfg(feature = "gzip")]
//! let structure = pdbrust::parser::gzip::parse_gzip_pdb_file("pdb1ubq.ent.gz")?;
//! ```

/// Parse mmCIF/PDBx files
pub mod mmcif;

/// Parse PDB files
pub mod pdb;

/// Parse gzip-compressed PDB/mmCIF files
#[cfg(feature = "gzip")]
pub mod gzip;

// Re-export all public parsing functions
pub use mmcif::{parse_mmcif_file, parse_mmcif_reader, parse_mmcif_string, parse_structure_file};
pub use pdb::*;

#[cfg(feature = "gzip")]
pub use gzip::{
    parse_gzip_mmcif_file, parse_gzip_mmcif_reader, parse_gzip_pdb_file, parse_gzip_pdb_reader,
    parse_gzip_structure_file,
};
