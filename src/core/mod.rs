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
//! - `mmcif_converter`: Converts mmCIF parsed data into PdbStructure format
//!   - Unified interface for both PDB and mmCIF data
//!   - Maintains compatibility with existing code
//!
//! # Examples
//!
//! ```no_run
//! use pdbrust::core::{PdbStructure};
//! use pdbrust::parser::{parse_pdb_file, parse_mmcif_file, parse_structure_file};
//!
//! // Parse a PDB file
//! let pdb_structure = parse_pdb_file("structure.pdb").unwrap();
//!
//! // Parse an mmCIF file
//! let mmcif_structure = parse_mmcif_file("structure.cif").unwrap();
//!
//! // Auto-detect format and parse
//! let structure = parse_structure_file("structure.ent").unwrap();
//! ```
//!
//! # Feature Flags
//!
//! - `pdb`: Enabled by default, provides PDB format support
//! - `mmcif`: Enabled by default, provides mmCIF format support

mod mmcif;
mod mmcif_converter;
mod pdb;

// Re-export the public interfaces
pub use mmcif::*;
pub use mmcif_converter::*;
pub use pdb::*;

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