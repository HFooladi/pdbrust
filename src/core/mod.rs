//! Core module for handling molecular structure file parsing and processing
//!
//! This module provides the core functionality for parsing and processing different
//! molecular structure file formats, including PDB and mmCIF files. It includes
//! structures and traits for handling atoms, chains, residues, and other molecular
//! components.
//!
//! # Module Organization
//!
//! - `pdb`: Handles Protein Data Bank (PDB) format parsing and representation
//!   - Supports standard PDB format (*.pdb files)
//!   - Handles ATOM, HETATM, SEQRES, and other record types
//!
//! - `mmcif`: Handles macromolecular Crystallographic Information File (mmCIF) format
//!   - Supports mmCIF format (*.cif files)
//!   - Provides category-based data access
//!   - Handles both loop_ and non-loop data structures
//!
//! # Examples
//!
//! ```no_run
//! use your_crate_name::core::{PdbParser, MmcifParser};
//!
//! // Parse a PDB file
//! let pdb_structure = PdbParser::parse_file("structure.pdb").unwrap();
//!
//! // Parse an mmCIF file
//! let mut mmcif_parser = MmcifParser::new();
//! mmcif_parser.parse_file("structure.cif").unwrap();
//! ```
//!
//! # Feature Flags
//!
//! - `pdb`: Enabled by default, provides PDB format support
//! - `mmcif`: Enabled by default, provides mmCIF format support

mod pdb;
mod mmcif;

// Re-export the public interfaces
pub use pdb::*;
pub use mmcif::*;

// Common traits and types that are shared between formats could be defined here
/// Represents a molecular structure, regardless of the source format
pub trait Structure {
    /// Get the number of atoms in the structure
    fn atom_count(&self) -> usize;
    
    /// Get the number of residues in the structure
    fn residue_count(&self) -> usize;
    
    /// Get the number of chains in the structure
    fn chain_count(&self) -> usize;
}