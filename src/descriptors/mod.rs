//! Structural descriptors and analysis functions for PDB structures.
//!
//! This module provides functionality for computing various structural
//! descriptors commonly used in structural bioinformatics:
//!
//! - **Composition analysis**: Amino acid frequencies, hydrophobicity ratios
//! - **Geometric descriptors**: Radius of gyration, compactness, density
//! - **Structural completeness**: Missing residue detection
//! - **Secondary structure**: Heuristic estimation based on Cα distances
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Composition analysis
//! let composition = structure.aa_composition();
//! let hydrophobic = structure.hydrophobic_ratio();
//!
//! // Geometric analysis
//! let rg = structure.radius_of_gyration();
//! let compactness = structure.compactness_index();
//!
//! // Get comprehensive summary
//! let summary = structure.structure_descriptors();
//! ```
//!
//! # Feature Flag
//!
//! This module requires the `descriptors` feature:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.1", features = ["descriptors"] }
//! ```

mod composition;
mod distance;
mod geometry;

// Re-export public constants and types
pub use composition::HYDROPHOBIC_RESIDUES;
pub use distance::{DEFAULT_ATOM_CONTACT_THRESHOLD, DEFAULT_CA_CONTACT_THRESHOLD};

use std::collections::HashMap;

/// Comprehensive structure descriptors combining composition and geometry.
///
/// This struct contains all computed descriptors for a protein structure,
/// suitable for batch processing or downstream analysis.
#[derive(Debug, Clone)]
pub struct StructureDescriptors {
    /// Number of residues (based on Cα count)
    pub num_residues: usize,
    /// Total number of atoms
    pub num_atoms: usize,
    /// Amino acid composition as fractions (0.0 to 1.0)
    pub aa_composition: HashMap<String, f64>,
    /// Fraction of glycine residues
    pub glycine_ratio: f64,
    /// Fraction of hydrophobic residues
    pub hydrophobic_ratio: f64,
    /// Radius of gyration (Å)
    pub radius_of_gyration: f64,
    /// Maximum Cα-Cα distance (Å)
    pub max_ca_distance: f64,
    /// Fraction of missing residues based on sequence gaps
    pub missing_residue_ratio: f64,
    /// Heuristic secondary structure content
    pub secondary_structure_ratio: f64,
    /// Compactness index: Rg / n^(1/3)
    pub compactness_index: f64,
    /// Cα density: count / bounding box volume
    pub ca_density: f64,
}

impl Default for StructureDescriptors {
    fn default() -> Self {
        Self {
            num_residues: 0,
            num_atoms: 0,
            aa_composition: HashMap::new(),
            glycine_ratio: 0.0,
            hydrophobic_ratio: 0.0,
            radius_of_gyration: 0.0,
            max_ca_distance: 0.0,
            missing_residue_ratio: 0.0,
            secondary_structure_ratio: 0.0,
            compactness_index: 0.0,
            ca_density: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_structure_descriptors_default() {
        let desc = StructureDescriptors::default();
        assert_eq!(desc.num_residues, 0);
        assert_eq!(desc.num_atoms, 0);
        assert!(desc.aa_composition.is_empty());
    }
}
