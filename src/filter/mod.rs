//! Filtering and cleaning operations for PDB structures.
//!
//! This module provides functionality for:
//! - Extracting specific atoms (Cα only, specific chains)
//! - Removing ligands and heteroatoms
//! - Normalizing chain identifiers
//! - Reindexing residues to continuous numbering
//! - Centering structures at the origin
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("example.pdb")?;
//!
//! // Extract Cα coordinates
//! let ca_coords = structure.get_ca_coords(None);
//!
//! // Get a cleaned structure with only protein atoms
//! let protein_only = structure.remove_ligands();
//!
//! // Chain operations with fluent API
//! let mut cleaned = structure.clone();
//! cleaned.normalize_chain_ids()
//!        .reindex_residues()
//!        .center_structure();
//! ```
//!
//! # Feature Flag
//!
//! This module is only available when the `filter` feature is enabled:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.7", features = ["filter"] }
//! ```

mod cleaning;
mod extraction;
pub mod selection;

// The cleaning and extraction modules extend PdbStructure with impl blocks.
// Selection module is public for accessing SelectionError and other types.

// Re-export classification constants and functions from the canonical source.
pub use crate::classify::{
    is_standard_amino_acid, is_standard_nucleotide, is_standard_residue, STANDARD_AMINO_ACIDS,
    STANDARD_NUCLEOTIDES,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_standard_amino_acid() {
        assert!(is_standard_amino_acid("ALA"));
        assert!(is_standard_amino_acid("GLY"));
        assert!(is_standard_amino_acid(" ALA "));
        assert!(!is_standard_amino_acid("HOH"));
        assert!(!is_standard_amino_acid("ATP"));
    }

    #[test]
    fn test_is_standard_nucleotide() {
        assert!(is_standard_nucleotide("A"));
        assert!(is_standard_nucleotide("DA"));
        assert!(!is_standard_nucleotide("ALA"));
    }

    #[test]
    fn test_is_standard_residue() {
        assert!(is_standard_residue("ALA"));
        assert!(is_standard_residue("A"));
        assert!(!is_standard_residue("HOH"));
    }
}
