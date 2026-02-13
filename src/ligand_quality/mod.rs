//! Ligand pose quality assessment module (PoseBusters-style geometry checks).
//!
//! This module provides functionality for validating the geometric quality of
//! protein-ligand complexes using PoseBusters-inspired criteria. It focuses on
//! intermolecular checks that don't require cheminformatics libraries like RDKit.
//!
//! # Features
//!
//! - **Steric clash detection**: Identifies atom pairs where the distance is below
//!   the expected van der Waals contact threshold (0.75 × sum of vdW radii)
//! - **Volume overlap calculation**: Estimates the percentage of ligand volume
//!   that overlaps with protein atoms
//! - **Cofactor clash detection**: Identifies clashes with non-protein heteroatoms
//!
//! # PoseBusters Context
//!
//! PoseBusters is a validation suite for protein-ligand complexes that has become
//! a standard benchmark in computational drug discovery and co-folding papers
//! (AlphaFold3, RoseTTAFold, DiffDock). Over 50% of AI-generated poses fail its
//! validity checks.
//!
//! This module implements the geometry-based subset of PoseBusters checks that
//! can be performed without RDKit:
//!
//! - Minimum protein-ligand distance (>0.75 × sum of vdW radii)
//! - Protein volume overlap (<7.5% at 0.8 vdW scaling)
//! - Cofactor distance checks
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("complex.pdb")?;
//!
//! // Check quality of a specific ligand
//! if let Some(report) = structure.ligand_pose_quality("LIG") {
//!     if report.has_protein_clash {
//!         println!("WARNING: {} clashes detected", report.num_clashes);
//!         for clash in &report.clashes {
//!             println!("  Distance: {:.2}Å < expected {:.2}Å",
//!                      clash.distance, clash.expected_min_distance);
//!         }
//!     }
//!
//!     if report.is_geometry_valid {
//!         println!("Pose passes geometry checks");
//!     }
//! }
//!
//! // Check all ligands at once
//! let reports = structure.all_ligand_pose_quality();
//! for report in &reports {
//!     println!("{}: {} (valid: {})",
//!              report.ligand_name,
//!              if report.is_geometry_valid { "PASS" } else { "FAIL" },
//!              report.is_geometry_valid);
//! }
//! ```
//!
//! # Feature Flag
//!
//! This module requires the `ligand-quality` feature:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.7", features = ["ligand-quality"] }
//! ```
//!
//! # References
//!
//! - [PoseBusters Paper (Chemical Science 2024)](https://pubs.rsc.org/en/content/articlehtml/2024/sc/d3sc04185a)
//! - [PoseBusters GitHub](https://github.com/maabuu/posebusters)

mod clash;
mod overlap;
mod radii;
mod report;

// Re-export public types
pub use clash::AtomClash;
pub use report::LigandPoseReport;

// Re-export radii functions for advanced users
pub use radii::{covalent_radius, vdw_radius};

use crate::core::PdbStructure;

/// Default clash threshold multiplier for van der Waals radii.
/// PoseBusters uses 0.75 × sum of vdW radii as the minimum allowed distance.
pub const CLASH_VDW_MULTIPLIER: f64 = 0.75;

/// Default volume overlap threshold in percent.
/// PoseBusters considers poses with >7.5% protein volume overlap as failures.
pub const MAX_VOLUME_OVERLAP_PCT: f64 = 7.5;

/// Default vdW scaling factor for volume overlap calculation.
/// PoseBusters uses 0.8 × vdW radii for volume overlap.
pub const VOLUME_VDW_SCALING: f64 = 0.8;

/// Grid spacing for volume overlap calculation (in Angstroms).
pub const VOLUME_GRID_SPACING: f64 = 0.5;

/// Standard water residue names to exclude from ligand checks.
pub const WATER_RESIDUES: [&str; 4] = ["HOH", "WAT", "H2O", "DOD"];

impl PdbStructure {
    /// Assess the pose quality of a specific ligand.
    ///
    /// Evaluates the geometric validity of a ligand pose using PoseBusters-style
    /// criteria including steric clash detection and volume overlap calculation.
    ///
    /// # Arguments
    ///
    /// * `ligand_name` - The residue name of the ligand (e.g., "LIG", "ATP", "NAG")
    ///
    /// # Returns
    ///
    /// `Some(LigandPoseReport)` if the ligand exists in the structure, `None` otherwise.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("complex.pdb")?;
    ///
    /// if let Some(report) = structure.ligand_pose_quality("LIG") {
    ///     println!("Min distance: {:.2} Å", report.min_protein_ligand_distance);
    ///     println!("Clashes: {}", report.num_clashes);
    ///     println!("Volume overlap: {:.1}%", report.protein_volume_overlap_pct);
    ///     println!("Valid: {}", report.is_geometry_valid);
    /// }
    /// ```
    pub fn ligand_pose_quality(&self, ligand_name: &str) -> Option<LigandPoseReport> {
        report::compute_ligand_pose_report(self, ligand_name)
    }

    /// Assess the pose quality of all ligands in the structure.
    ///
    /// Evaluates the geometric validity of all HETATM residues (excluding water)
    /// using PoseBusters-style criteria.
    ///
    /// # Returns
    ///
    /// A vector of `LigandPoseReport` for each unique ligand in the structure.
    /// Returns an empty vector if no ligands are present.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("complex.pdb")?;
    ///
    /// let reports = structure.all_ligand_pose_quality();
    /// for report in &reports {
    ///     let status = if report.is_geometry_valid { "PASS" } else { "FAIL" };
    ///     println!("{}: {} ({} clashes, {:.1}% overlap)",
    ///              report.ligand_name, status,
    ///              report.num_clashes, report.protein_volume_overlap_pct);
    /// }
    /// ```
    pub fn all_ligand_pose_quality(&self) -> Vec<LigandPoseReport> {
        report::compute_all_ligand_reports(self)
    }

    /// Get a list of all unique ligand names in the structure.
    ///
    /// Returns the residue names of all HETATM residues excluding water.
    /// This is useful for discovering which ligands are present before
    /// running quality checks.
    ///
    /// # Returns
    ///
    /// A sorted vector of unique ligand residue names.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("complex.pdb")?;
    ///
    /// let ligands = structure.get_ligand_names();
    /// println!("Found {} ligands: {:?}", ligands.len(), ligands);
    /// ```
    pub fn get_ligand_names(&self) -> Vec<String> {
        use std::collections::HashSet;

        let mut ligand_names: HashSet<String> = self
            .atoms
            .iter()
            .filter(|atom| atom.is_hetatm)
            .filter(|atom| !WATER_RESIDUES.contains(&atom.residue_name.as_str()))
            .map(|atom| atom.residue_name.clone())
            .collect();

        let mut names: Vec<String> = ligand_names.drain().collect();
        names.sort();
        names
    }

    /// Check if a residue name corresponds to water.
    ///
    /// # Arguments
    ///
    /// * `residue_name` - The residue name to check
    ///
    /// # Returns
    ///
    /// `true` if the residue name is a known water identifier.
    pub fn is_water_residue(residue_name: &str) -> bool {
        WATER_RESIDUES.contains(&residue_name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    #[allow(clippy::too_many_arguments)]
    fn create_test_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
        element: &str,
        is_hetatm: bool,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    fn create_simple_complex() -> PdbStructure {
        let mut structure = PdbStructure::new();

        // Add protein atoms (chain A)
        structure.atoms.push(create_test_atom(
            1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", false,
        ));
        structure.atoms.push(create_test_atom(
            2, "CB", "ALA", "A", 1, 1.5, 0.0, 0.0, "C", false,
        ));
        structure.atoms.push(create_test_atom(
            3, "N", "ALA", "A", 1, -1.5, 0.0, 0.0, "N", false,
        ));
        structure.atoms.push(create_test_atom(
            4, "O", "ALA", "A", 1, 0.0, 1.5, 0.0, "O", false,
        ));

        // Add ligand atoms (good pose - no clashes)
        structure.atoms.push(create_test_atom(
            5, "C1", "LIG", "A", 100, 5.0, 5.0, 5.0, "C", true,
        ));
        structure.atoms.push(create_test_atom(
            6, "O1", "LIG", "A", 100, 6.0, 5.0, 5.0, "O", true,
        ));
        structure.atoms.push(create_test_atom(
            7, "N1", "LIG", "A", 100, 5.0, 6.0, 5.0, "N", true,
        ));

        // Add water (should be excluded)
        structure.atoms.push(create_test_atom(
            8, "O", "HOH", "A", 200, 10.0, 10.0, 10.0, "O", true,
        ));

        structure
    }

    #[test]
    fn test_get_ligand_names() {
        let structure = create_simple_complex();
        let ligands = structure.get_ligand_names();

        assert_eq!(ligands.len(), 1);
        assert!(ligands.contains(&"LIG".to_string()));
        assert!(!ligands.contains(&"HOH".to_string()));
    }

    #[test]
    fn test_is_water_residue() {
        assert!(PdbStructure::is_water_residue("HOH"));
        assert!(PdbStructure::is_water_residue("WAT"));
        assert!(PdbStructure::is_water_residue("H2O"));
        assert!(PdbStructure::is_water_residue("DOD"));
        assert!(!PdbStructure::is_water_residue("LIG"));
        assert!(!PdbStructure::is_water_residue("ATP"));
    }

    #[test]
    fn test_ligand_pose_quality_good_pose() {
        let structure = create_simple_complex();
        let report = structure.ligand_pose_quality("LIG");

        assert!(report.is_some());
        let report = report.unwrap();

        assert_eq!(report.ligand_name, "LIG");
        assert!(!report.has_protein_clash);
        assert_eq!(report.num_clashes, 0);
        assert!(report.is_geometry_valid);
    }

    #[test]
    fn test_ligand_pose_quality_nonexistent_ligand() {
        let structure = create_simple_complex();
        let report = structure.ligand_pose_quality("XYZ");

        assert!(report.is_none());
    }

    #[test]
    fn test_ligand_pose_quality_with_clash() {
        let mut structure = create_simple_complex();

        // Add a clashing ligand atom very close to the protein
        structure.atoms.push(create_test_atom(
            9, "C2", "BAD", "A", 101, 0.5, 0.0, 0.0, "C", true,
        ));

        let report = structure.ligand_pose_quality("BAD");

        assert!(report.is_some());
        let report = report.unwrap();

        assert!(report.has_protein_clash);
        assert!(report.num_clashes > 0);
        assert!(!report.is_geometry_valid);
    }

    #[test]
    fn test_all_ligand_pose_quality() {
        let structure = create_simple_complex();
        let reports = structure.all_ligand_pose_quality();

        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].ligand_name, "LIG");
    }

    #[test]
    fn test_empty_structure() {
        let structure = PdbStructure::new();
        let ligands = structure.get_ligand_names();
        let reports = structure.all_ligand_pose_quality();

        assert!(ligands.is_empty());
        assert!(reports.is_empty());
    }
}
