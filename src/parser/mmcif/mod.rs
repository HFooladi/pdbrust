//! mmCIF (macromolecular Crystallographic Information File) parsing module
//!
//! This module provides functionality to parse mmCIF format files and convert
//! them to the unified PdbStructure format used throughout the library.
//!
//! mmCIF is the standard format used by the Protein Data Bank (PDB) for
//! distributing structural data. It's more flexible and comprehensive than
//! the legacy PDB format.
//!
//! # Examples
//!
//! ```no_run
//! use pdbrust::parser::mmcif::parse_mmcif_file;
//!
//! let structure = parse_mmcif_file("structure.cif")?;
//! println!("Loaded {} atoms", structure.atoms.len());
//! ```

mod parser;

pub use parser::*;