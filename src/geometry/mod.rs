//! Geometric analysis and structure superposition.
//!
//! This module provides RMSD calculation and structure alignment using
//! the Kabsch algorithm. It requires the `geometry` feature flag.
//!
//! # Feature Flag
//!
//! Enable this module in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.6", features = ["geometry"] }
//! ```
//!
//! # Overview
//!
//! The main functionality includes:
//!
//! - **RMSD Calculation**: Compute root mean square deviation between structures
//! - **Structure Alignment**: Optimal superposition using Kabsch algorithm
//! - **Per-Residue RMSD**: Identify flexible regions in proteins
//!
//! # Quick Start
//!
//! ```rust,ignore
//! use pdbrust::{parse_pdb_file, PdbStructure};
//! use pdbrust::geometry::AtomSelection;
//!
//! // Parse two structures
//! let structure1 = parse_pdb_file("model1.pdb")?;
//! let structure2 = parse_pdb_file("model2.pdb")?;
//!
//! // Calculate RMSD (without alignment)
//! let rmsd = structure1.rmsd_to(&structure2)?;
//! println!("Direct RMSD: {:.2} Angstroms", rmsd);
//!
//! // Align structures and get aligned RMSD
//! let (aligned, result) = structure1.align_to(&structure2)?;
//! println!("Aligned RMSD: {:.4} Angstroms", result.rmsd);
//!
//! // Get per-residue RMSD for flexibility analysis
//! let per_res = structure1.per_residue_rmsd_to(&structure2)?;
//! for r in &per_res {
//!     if r.rmsd > 2.0 {
//!         println!("Flexible: {}{} - {:.2} A", r.residue_id.0, r.residue_id.1, r.rmsd);
//!     }
//! }
//! ```
//!
//! # Atom Selection
//!
//! By default, RMSD and alignment use CA (alpha-carbon) atoms only.
//! You can customize this with [`AtomSelection`]:
//!
//! ```rust,ignore
//! use pdbrust::geometry::AtomSelection;
//!
//! // Use backbone atoms (N, CA, C, O)
//! let rmsd = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
//!
//! // Use all atoms (requires exact correspondence)
//! let rmsd = structure1.rmsd_to_with_selection(&structure2, AtomSelection::AllAtoms)?;
//!
//! // Use custom atom set
//! let selection = AtomSelection::Custom(vec!["CA".into(), "CB".into()]);
//! let rmsd = structure1.rmsd_to_with_selection(&structure2, selection)?;
//! ```

mod rmsd;
mod superposition;
mod transform;

// Re-export public types
pub use rmsd::{calculate_rmsd, calculate_rmsd_chain, rmsd_from_coords};
pub use superposition::{
    AlignmentResult, PerResidueRmsd, align_structures, calculate_alignment, per_residue_rmsd,
};
pub use transform::{
    AtomSelection, apply_transform, apply_transform_to_coords, center_coords, compute_centroid,
    extract_coords_by_selection, extract_coords_with_residue_info, translate_structure,
};

use crate::core::PdbStructure;
use crate::error::PdbError;

/// Extension methods for PdbStructure providing RMSD and alignment functionality.
impl PdbStructure {
    /// Calculate RMSD to another structure using CA atoms.
    ///
    /// This is a convenience method that uses the default CA-only atom selection.
    /// For other atom selections, use [`rmsd_to_with_selection`](Self::rmsd_to_with_selection).
    ///
    /// **Note**: This calculates RMSD without alignment. For aligned RMSD, use
    /// [`align_to`](Self::align_to) which returns both the aligned structure and RMSD.
    ///
    /// # Arguments
    /// * `other` - The structure to compare against
    ///
    /// # Returns
    /// RMSD in Angstroms
    ///
    /// # Errors
    /// - `AtomCountMismatch` if structures have different numbers of CA atoms
    /// - `NoAtomsSelected` if either structure has no CA atoms
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let rmsd = structure1.rmsd_to(&structure2)?;
    /// println!("RMSD: {:.2} Angstroms", rmsd);
    /// ```
    pub fn rmsd_to(&self, other: &PdbStructure) -> Result<f64, PdbError> {
        calculate_rmsd(self, other, AtomSelection::CaOnly)
    }

    /// Calculate RMSD with a custom atom selection.
    ///
    /// # Arguments
    /// * `other` - The structure to compare against
    /// * `selection` - Atom selection criteria
    ///
    /// # Returns
    /// RMSD in Angstroms
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use pdbrust::geometry::AtomSelection;
    ///
    /// // RMSD using backbone atoms
    /// let rmsd = structure1.rmsd_to_with_selection(&structure2, AtomSelection::Backbone)?;
    /// ```
    pub fn rmsd_to_with_selection(
        &self,
        other: &PdbStructure,
        selection: AtomSelection,
    ) -> Result<f64, PdbError> {
        calculate_rmsd(self, other, selection)
    }

    /// Align this structure onto a target structure.
    ///
    /// Uses the Kabsch algorithm to find the optimal rotation and translation
    /// that minimizes RMSD between the structures. Returns both the aligned
    /// structure and an [`AlignmentResult`] containing the RMSD and transformation.
    ///
    /// # Arguments
    /// * `target` - The reference structure to align onto
    ///
    /// # Returns
    /// Tuple of (aligned_structure, AlignmentResult)
    ///
    /// # Errors
    /// - `AtomCountMismatch` if structures have different numbers of CA atoms
    /// - `NoAtomsSelected` if either structure has no CA atoms
    /// - `InsufficientAtoms` if fewer than 3 CA atoms are present
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let (aligned, result) = mobile.align_to(&target)?;
    /// println!("RMSD after alignment: {:.4} Angstroms", result.rmsd);
    /// println!("Aligned {} atoms", result.num_atoms);
    ///
    /// // Save the aligned structure
    /// aligned.to_file("aligned.pdb")?;
    /// ```
    pub fn align_to(&self, target: &PdbStructure) -> Result<(Self, AlignmentResult), PdbError> {
        align_structures(self, target, AtomSelection::CaOnly)
    }

    /// Align with a custom atom selection.
    ///
    /// # Arguments
    /// * `target` - The reference structure to align onto
    /// * `selection` - Atom selection criteria for alignment
    ///
    /// # Returns
    /// Tuple of (aligned_structure, AlignmentResult)
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use pdbrust::geometry::AtomSelection;
    ///
    /// // Align using backbone atoms
    /// let (aligned, result) = mobile.align_to_with_selection(
    ///     &target,
    ///     AtomSelection::Backbone
    /// )?;
    /// ```
    pub fn align_to_with_selection(
        &self,
        target: &PdbStructure,
        selection: AtomSelection,
    ) -> Result<(Self, AlignmentResult), PdbError> {
        align_structures(self, target, selection)
    }

    /// Get per-residue RMSD after alignment to target.
    ///
    /// This method first aligns the structures optimally, then computes
    /// RMSD for each residue individually. Useful for identifying flexible
    /// regions in protein structures.
    ///
    /// # Arguments
    /// * `target` - The reference structure
    ///
    /// # Returns
    /// Vector of [`PerResidueRmsd`] containing RMSD for each residue
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let per_res = mobile.per_residue_rmsd_to(&target)?;
    ///
    /// // Find flexible regions (RMSD > 2.0 Angstroms)
    /// for r in &per_res {
    ///     if r.rmsd > 2.0 {
    ///         println!("{}{} {}: {:.2} A",
    ///             r.residue_id.0,
    ///             r.residue_id.1,
    ///             r.residue_name,
    ///             r.rmsd
    ///         );
    ///     }
    /// }
    /// ```
    pub fn per_residue_rmsd_to(
        &self,
        target: &PdbStructure,
    ) -> Result<Vec<PerResidueRmsd>, PdbError> {
        per_residue_rmsd(self, target, AtomSelection::CaOnly)
    }

    /// Get per-residue RMSD with a custom atom selection.
    ///
    /// # Arguments
    /// * `target` - The reference structure
    /// * `selection` - Atom selection criteria
    ///
    /// # Returns
    /// Vector of [`PerResidueRmsd`] containing RMSD for each residue
    pub fn per_residue_rmsd_to_with_selection(
        &self,
        target: &PdbStructure,
        selection: AtomSelection,
    ) -> Result<Vec<PerResidueRmsd>, PdbError> {
        per_residue_rmsd(self, target, selection)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_atom(x: f64, y: f64, z: f64, residue_seq: i32) -> Atom {
        Atom {
            serial: residue_seq,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
            residue_seq,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        }
    }

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(0.0, 0.0, 0.0, 1),
            create_atom(3.8, 0.0, 0.0, 2),
            create_atom(7.6, 0.0, 0.0, 3),
            create_atom(11.4, 0.0, 0.0, 4),
        ];
        structure
    }

    #[test]
    fn test_pdbstructure_rmsd_to_self() {
        let structure = create_test_structure();
        let rmsd = structure.rmsd_to(&structure).unwrap();
        assert!(rmsd < 1e-10);
    }

    #[test]
    fn test_pdbstructure_align_to_self() {
        let structure = create_test_structure();
        let (aligned, result) = structure.align_to(&structure).unwrap();

        assert!(result.rmsd < 1e-10);
        assert_eq!(result.num_atoms, 4);
        assert_eq!(aligned.atoms.len(), structure.atoms.len());
    }

    #[test]
    fn test_pdbstructure_align_translated() {
        let target = create_test_structure();
        let mut mobile = target.clone();
        for atom in &mut mobile.atoms {
            atom.x += 100.0;
            atom.y += 100.0;
            atom.z += 100.0;
        }

        let (aligned, result) = mobile.align_to(&target).unwrap();

        // RMSD should be ~0 after alignment
        assert!(result.rmsd < 1e-6);

        // Aligned atoms should match target
        for (a, t) in aligned.atoms.iter().zip(target.atoms.iter()) {
            assert!((a.x - t.x).abs() < 1e-6);
            assert!((a.y - t.y).abs() < 1e-6);
            assert!((a.z - t.z).abs() < 1e-6);
        }
    }

    #[test]
    fn test_pdbstructure_per_residue_rmsd() {
        let target = create_test_structure();
        let mut mobile = target.clone();
        for atom in &mut mobile.atoms {
            atom.x += 10.0;
        }

        let per_res = mobile.per_residue_rmsd_to(&target).unwrap();

        assert_eq!(per_res.len(), 4);
        for r in &per_res {
            // After alignment, all per-residue RMSDs should be ~0
            assert!(r.rmsd < 1e-6);
        }
    }
}
