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
//! pdbrust = { version = "0.1", features = ["filter"] }
//! ```

mod cleaning;
mod extraction;

// The cleaning and extraction modules extend PdbStructure with impl blocks.
// No standalone items need to be re-exported.

/// Standard amino acid residue names (3-letter codes).
pub const STANDARD_AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Standard nucleotide residue names.
pub const STANDARD_NUCLEOTIDES: &[&str] = &[
    "A", "C", "G", "U", // RNA
    "DA", "DC", "DG", "DT", // DNA
];

/// Check if a residue name is a standard amino acid.
#[inline]
pub fn is_standard_amino_acid(residue_name: &str) -> bool {
    STANDARD_AMINO_ACIDS.contains(&residue_name.trim())
}

/// Check if a residue name is a standard nucleotide.
#[inline]
pub fn is_standard_nucleotide(residue_name: &str) -> bool {
    STANDARD_NUCLEOTIDES.contains(&residue_name.trim())
}

/// Check if a residue name is a standard biological residue (amino acid or nucleotide).
#[inline]
pub fn is_standard_residue(residue_name: &str) -> bool {
    is_standard_amino_acid(residue_name) || is_standard_nucleotide(residue_name)
}

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
