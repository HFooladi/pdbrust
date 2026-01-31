//! Geometric analysis and structure superposition.
//!
//! This module provides RMSD calculation, LDDT scoring, and structure alignment
//! using the Kabsch algorithm. It requires the `geometry` feature flag.
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
//! - **LDDT Calculation**: Compute Local Distance Difference Test (superposition-free)
//! - **Structure Alignment**: Optimal superposition using Kabsch algorithm
//! - **Per-Residue Analysis**: Identify flexible or poorly modeled regions
//!
//! # Quick Start
//!
//! ```rust,ignore
//! use pdbrust::{parse_pdb_file, PdbStructure};
//! use pdbrust::geometry::{AtomSelection, LddtOptions};
//!
//! // Parse two structures
//! let model = parse_pdb_file("model.pdb")?;
//! let reference = parse_pdb_file("reference.pdb")?;
//!
//! // Calculate RMSD (without alignment)
//! let rmsd = model.rmsd_to(&reference)?;
//! println!("Direct RMSD: {:.2} Angstroms", rmsd);
//!
//! // Calculate LDDT (superposition-free)
//! let lddt = model.lddt_to(&reference)?;
//! println!("LDDT: {:.4}", lddt.score);
//!
//! // Align structures and get aligned RMSD
//! let (aligned, result) = model.align_to(&reference)?;
//! println!("Aligned RMSD: {:.4} Angstroms", result.rmsd);
//!
//! // Get per-residue RMSD for flexibility analysis
//! let per_res = model.per_residue_rmsd_to(&reference)?;
//! for r in &per_res {
//!     if r.rmsd > 2.0 {
//!         println!("Flexible: {}{} - {:.2} A", r.residue_id.0, r.residue_id.1, r.rmsd);
//!     }
//! }
//!
//! // Get per-residue LDDT for quality analysis
//! let per_res_lddt = model.per_residue_lddt_to(&reference)?;
//! for r in per_res_lddt.iter().filter(|r| r.score < 0.7) {
//!     println!("Low LDDT: {}{} = {:.2}", r.residue_id.0, r.residue_id.1, r.score);
//! }
//! ```
//!
//! # Atom Selection
//!
//! By default, RMSD, LDDT, and alignment use CA (alpha-carbon) atoms only.
//! You can customize this with [`AtomSelection`]:
//!
//! ```rust,ignore
//! use pdbrust::geometry::AtomSelection;
//!
//! // Use backbone atoms (N, CA, C, O)
//! let rmsd = model.rmsd_to_with_selection(&reference, AtomSelection::Backbone)?;
//!
//! // Use all atoms (requires exact correspondence)
//! let rmsd = model.rmsd_to_with_selection(&reference, AtomSelection::AllAtoms)?;
//!
//! // Use custom atom set
//! let selection = AtomSelection::Custom(vec!["CA".into(), "CB".into()]);
//! let rmsd = model.rmsd_to_with_selection(&reference, selection)?;
//! ```
//!
//! # LDDT vs RMSD
//!
//! LDDT (Local Distance Difference Test) complements RMSD:
//!
//! - **LDDT** is superposition-free and focuses on local structure quality
//! - **RMSD** requires alignment and measures global deviation
//! - LDDT is used in AlphaFold (pLDDT) and CASP evaluations
//! - For comparing predicted structures, LDDT is often preferred

mod lddt;
mod rmsd;
mod superposition;
mod transform;

// Re-export public types
pub use lddt::{LddtOptions, LddtResult, PerResidueLddt, calculate_lddt, per_residue_lddt};
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

    // ==================== LDDT Methods ====================

    /// Calculate LDDT (Local Distance Difference Test) to a reference structure.
    ///
    /// LDDT is a superposition-free metric that measures the fraction of
    /// inter-atomic distances that are preserved within specified thresholds.
    /// It ranges from 0.0 (poor) to 1.0 (perfect).
    ///
    /// This is the same metric used by AlphaFold (pLDDT) and CASP evaluations.
    ///
    /// # Arguments
    /// * `reference` - The reference structure (ground truth)
    ///
    /// # Returns
    /// [`LddtResult`] containing the global score and per-threshold scores
    ///
    /// # Errors
    /// - `AtomCountMismatch` if structures have different numbers of CA atoms
    /// - `NoAtomsSelected` if either structure has no CA atoms
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let result = model.lddt_to(&reference)?;
    /// println!("LDDT: {:.4}", result.score);
    /// println!("Evaluated {} distance pairs", result.num_pairs);
    /// ```
    pub fn lddt_to(&self, reference: &PdbStructure) -> Result<LddtResult, PdbError> {
        calculate_lddt(
            self,
            reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
    }

    /// Calculate LDDT with custom atom selection and options.
    ///
    /// # Arguments
    /// * `reference` - The reference structure (ground truth)
    /// * `selection` - Atom selection criteria (e.g., CA only, backbone)
    /// * `options` - LDDT options (inclusion radius, thresholds)
    ///
    /// # Returns
    /// [`LddtResult`] containing the global score and per-threshold scores
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use pdbrust::geometry::{AtomSelection, LddtOptions};
    ///
    /// let options = LddtOptions::default()
    ///     .with_inclusion_radius(10.0)
    ///     .with_thresholds(vec![0.5, 1.0, 2.0, 4.0]);
    ///
    /// let result = model.lddt_to_with_options(
    ///     &reference,
    ///     AtomSelection::Backbone,
    ///     options
    /// )?;
    /// ```
    pub fn lddt_to_with_options(
        &self,
        reference: &PdbStructure,
        selection: AtomSelection,
        options: LddtOptions,
    ) -> Result<LddtResult, PdbError> {
        calculate_lddt(self, reference, selection, options)
    }

    /// Get per-residue LDDT scores.
    ///
    /// Returns LDDT scores for each residue individually, useful for
    /// identifying poorly modeled regions in predicted structures.
    ///
    /// # Arguments
    /// * `reference` - The reference structure (ground truth)
    ///
    /// # Returns
    /// Vector of [`PerResidueLddt`] containing LDDT for each residue
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let per_res = model.per_residue_lddt_to(&reference)?;
    ///
    /// // Find poorly modeled regions (LDDT < 0.7)
    /// for r in per_res.iter().filter(|r| r.score < 0.7) {
    ///     println!("{}{} {}: LDDT = {:.2}",
    ///         r.residue_id.0,
    ///         r.residue_id.1,
    ///         r.residue_name,
    ///         r.score
    ///     );
    /// }
    /// ```
    pub fn per_residue_lddt_to(
        &self,
        reference: &PdbStructure,
    ) -> Result<Vec<PerResidueLddt>, PdbError> {
        per_residue_lddt(
            self,
            reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
    }

    /// Get per-residue LDDT with custom atom selection and options.
    ///
    /// # Arguments
    /// * `reference` - The reference structure (ground truth)
    /// * `selection` - Atom selection criteria
    /// * `options` - LDDT options (inclusion radius, thresholds)
    ///
    /// # Returns
    /// Vector of [`PerResidueLddt`] containing LDDT for each residue
    pub fn per_residue_lddt_to_with_options(
        &self,
        reference: &PdbStructure,
        selection: AtomSelection,
        options: LddtOptions,
    ) -> Result<Vec<PerResidueLddt>, PdbError> {
        per_residue_lddt(self, reference, selection, options)
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
            is_hetatm: false,
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

    // ==================== LDDT Tests ====================

    #[test]
    fn test_pdbstructure_lddt_to_self() {
        let structure = create_test_structure();
        let result = structure.lddt_to(&structure).unwrap();

        // Self-comparison should have perfect LDDT
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "Self-LDDT should be 1.0, got {}",
            result.score
        );
    }

    #[test]
    fn test_pdbstructure_lddt_translation_invariant() {
        let reference = create_test_structure();
        let mut model = reference.clone();

        // Translate the model
        for atom in &mut model.atoms {
            atom.x += 100.0;
            atom.y += 50.0;
            atom.z += 25.0;
        }

        let result = model.lddt_to(&reference).unwrap();

        // LDDT should be invariant to translation
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "LDDT should be translation invariant, got {}",
            result.score
        );
    }

    #[test]
    fn test_pdbstructure_lddt_rotation_invariant() {
        let reference = create_test_structure();
        let mut model = reference.clone();

        // Rotate 90 degrees around z-axis
        for atom in &mut model.atoms {
            let x = atom.x;
            let y = atom.y;
            atom.x = -y;
            atom.y = x;
        }

        let result = model.lddt_to(&reference).unwrap();

        // LDDT should be invariant to rotation
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "LDDT should be rotation invariant, got {}",
            result.score
        );
    }

    #[test]
    fn test_pdbstructure_lddt_perturbed() {
        let reference = create_test_structure();
        let mut model = reference.clone();

        // Perturb one atom
        model.atoms[1].y += 5.0;

        let result = model.lddt_to(&reference).unwrap();

        // Perturbed structure should have LDDT < 1.0
        assert!(
            result.score < 1.0,
            "Perturbed LDDT should be < 1.0, got {}",
            result.score
        );
        assert!(
            result.score > 0.0,
            "Perturbed LDDT should be > 0.0, got {}",
            result.score
        );
    }

    #[test]
    fn test_pdbstructure_per_residue_lddt() {
        let reference = create_test_structure();
        let model = reference.clone();

        let per_res = model.per_residue_lddt_to(&reference).unwrap();

        assert_eq!(per_res.len(), 4, "Should have 4 residues");
        for r in &per_res {
            assert!(
                (r.score - 1.0).abs() < 1e-10,
                "Self per-residue LDDT should be 1.0"
            );
        }
    }
}
