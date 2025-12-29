//! PDBRust: A Rust library for parsing and analyzing Protein Data Bank (PDB) files
//!
//! This library provides a robust and efficient way to work with both PDB and mmCIF files, 
//! supporting various record types including ATOM, SEQRES, CONECT, SSBOND, and more. 
//! It offers functionality for structural analysis, sequence information retrieval, 
//! and connectivity analysis of molecular structures.
//!
//! For detailed documentation, examples, and best practices, see the [guide] module.
//!
//! # Features
//!
//! - Parse both PDB and mmCIF files with comprehensive error handling
//! - Automatic format detection
//! - Support for multiple models in a single structure
//! - Chain and residue analysis
//! - Connectivity information through CONECT records
//! - Sequence information through SEQRES records
//! - Support for disulfide bonds through SSBOND records
//! - Remark handling for additional structural information
//!
//! # Format Support
//!
//! ## PDB Format
//! Traditional fixed-width text format with support for:
//! - ATOM/HETATM records
//! - HEADER, TITLE, REMARK records
//! - SEQRES, CONECT, SSBOND records
//! - MODEL/ENDMDL for multi-model structures
//!
//! ## mmCIF Format
//! Modern dictionary-based format with support for:
//! - _atom_site category (converted to ATOM records)
//! - _entity_poly_seq category (converted to SEQRES records)
//! - _struct_disulfid category (converted to SSBOND records)
//! - Header and metadata information
//!
//! # Examples
//!
//! ## Basic Usage
//!
//! ```rust,no_run
//! use pdbrust::{PdbStructure, parse_structure_file};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Auto-detect format and parse (works with both .pdb and .cif files)
//!     let structure = parse_structure_file("example.cif")?;
//!     
//!     // Get all chain IDs in the structure
//!     let chains = structure.get_chain_ids();
//!     
//!     // Get sequence for a specific chain
//!     if let Some(chain_id) = chains.first() {
//!         let sequence = structure.get_sequence(chain_id);
//!         println!("Sequence for chain {}: {:?}", chain_id, sequence);
//!     }
//!     
//!     Ok(())
//! }
//! ```
//!
//! ## Format-Specific Parsing
//!
//! ```rust,no_run
//! use pdbrust::{parse_pdb_file, parse_mmcif_file};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Parse PDB format explicitly
//!     let pdb_structure = parse_pdb_file("structure.pdb")?;
//!
//!     // Parse mmCIF format explicitly
//!     let mmcif_structure = parse_mmcif_file("structure.cif")?;
//!
//!     // Both produce the same PdbStructure type
//!     assert_eq!(pdb_structure.atoms.len(), mmcif_structure.atoms.len());
//!     Ok(())
//! }
//! ```

// Module exports
pub mod core;
pub mod error;
pub mod guide;
pub mod parser;
pub mod records;
mod utils;
pub mod writer;

// Conditional modules
#[cfg(feature = "filter")]
pub mod filter;

#[cfg(feature = "descriptors")]
pub mod descriptors;

#[cfg(feature = "quality")]
pub mod quality;

#[cfg(feature = "summary")]
pub mod summary;

#[cfg(feature = "rcsb")]
pub mod rcsb;

// Re-exports for convenience
pub use core::PdbStructure;
pub use error::PdbError;
pub use parser::{parse_pdb_file, parse_pdb_string, parse_mmcif_file, parse_mmcif_string, parse_structure_file};
pub use records::{Atom, Conect, Model, Remark, Residue, SSBond, SeqRes};
pub use writer::write_pdb_file;