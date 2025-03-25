//! PDBRust: A Rust library for parsing and analyzing Protein Data Bank (PDB) files
//!
//! This library provides a robust and efficient way to work with PDB files, supporting various record
//! types including ATOM, SEQRES, CONECT, SSBOND, and more. It offers functionality for structural
//! analysis, sequence information retrieval, and connectivity analysis of molecular structures.
//!
//! For detailed documentation, examples, and best practices, see the [guide](guide) module.
//!
//! # Features
//!
//! - Parse PDB files with comprehensive error handling
//! - Support for multiple models in a single structure
//! - Chain and residue analysis
//! - Connectivity information through CONECT records
//! - Sequence information through SEQRES records
//! - Support for disulfide bonds through SSBOND records
//! - Remark handling for additional structural information
//!
//! # Example
//!
//! ```rust
//! use pdbrust::PdbStructure;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let structure = PdbStructure::from_file("example.pdb")?;
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

// Module exports
pub mod core;
pub mod error;
pub mod guide;
pub mod parser;
pub mod records;
mod utils;
pub mod writer;

// Re-exports for convenience
pub use core::PdbStructure;
pub use error::PdbError;
pub use parser::parse_pdb_file;
pub use records::{Atom, Conect, Model, Remark, SSBond, SeqRes};
pub use writer::write_pdb_file;
