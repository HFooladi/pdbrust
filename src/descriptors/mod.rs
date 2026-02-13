//! Structural descriptors and analysis functions for PDB structures.
//!
//! This module provides functionality for computing various structural
//! descriptors commonly used in structural bioinformatics:
//!
//! - **Composition analysis**: Amino acid frequencies, hydrophobicity ratios
//! - **Geometric descriptors**: Radius of gyration, compactness, density
//! - **Structural completeness**: Missing residue detection
//! - **Secondary structure**: Heuristic estimation based on Cα distances
//! - **B-factor analysis**: Temperature factor statistics and flexibility identification
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
//! // B-factor analysis
//! let mean_b = structure.b_factor_mean();
//! let profile = structure.b_factor_profile();
//! let flexible = structure.flexible_residues(50.0);
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
//! pdbrust = { version = "0.7", features = ["descriptors"] }
//! ```

mod bfactor;
mod composition;
mod distance;
mod geometry;

#[cfg(feature = "dssp")]
mod dihedrals;

#[cfg(feature = "dssp")]
mod hbonds;

mod alphafold;
mod interactions;

// Re-export public constants and types
pub use composition::HYDROPHOBIC_RESIDUES;
pub use distance::{DEFAULT_ATOM_CONTACT_THRESHOLD, DEFAULT_CA_CONTACT_THRESHOLD};

// Re-export dihedral types (requires dssp feature)
#[cfg(feature = "dssp")]
pub use dihedrals::{RamachandranRegion, RamachandranStats, ResidueDihedrals, ResidueRef};

// Re-export hydrogen bond types (requires dssp feature)
#[cfg(feature = "dssp")]
pub use hbonds::{HBondStats, HBondType, MainchainHBond, ResidueHBonds};

// Re-export AlphaFold/pLDDT types
pub use alphafold::{ConfidenceCategory, ResiduePlddt};

// Re-export protein-ligand interaction types
pub use interactions::{
    BindingSite, ContactResidue, HydrophobicContact, LigandInteractionProfile, ProteinLigandHBond,
    SaltBridge,
};

use std::collections::HashMap;

/// Per-residue B-factor statistics.
///
/// Contains the mean, min, max B-factors and atom count for a single residue.
/// Used by `PdbStructure::b_factor_profile()` and related methods.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueBFactor {
    /// Chain identifier (e.g., "A")
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
    /// Residue name (e.g., "ALA", "GLY")
    pub residue_name: String,
    /// Mean B-factor across all atoms in the residue (Å²)
    pub b_factor_mean: f64,
    /// Minimum B-factor among atoms in the residue (Å²)
    pub b_factor_min: f64,
    /// Maximum B-factor among atoms in the residue (Å²)
    pub b_factor_max: f64,
    /// Number of atoms in the residue
    pub atom_count: usize,
}

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
    /// Mean B-factor (temperature factor) across all atoms (Å²)
    pub b_factor_mean: f64,
    /// Mean B-factor for Cα atoms only (Å²)
    pub b_factor_mean_ca: f64,
    /// Minimum B-factor in the structure (Å²)
    pub b_factor_min: f64,
    /// Maximum B-factor in the structure (Å²)
    pub b_factor_max: f64,
    /// Standard deviation of B-factors (Å²)
    pub b_factor_std: f64,
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
            b_factor_mean: 0.0,
            b_factor_mean_ca: 0.0,
            b_factor_min: 0.0,
            b_factor_max: 0.0,
            b_factor_std: 0.0,
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
