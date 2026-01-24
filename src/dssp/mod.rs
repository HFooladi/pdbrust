//! DSSP 4-like secondary structure assignment.
//!
//! This module implements a secondary structure assignment algorithm based on the
//! DSSP 4 method (Kabsch & Sander with updates from Hekkelman et al., 2025).
//!
//! The algorithm computes secondary structure from atomic coordinates using:
//! - Hydrogen bond detection via electrostatic energy calculation
//! - Pattern recognition for helices, sheets, turns, and bends
//! - Backbone dihedral analysis for PPII (κ-helix) detection
//!
//! # Secondary Structure Types
//!
//! The following 9-state classification is used:
//!
//! | Code | Name | Description |
//! |------|------|-------------|
//! | H | α-helix | i → i+4 hydrogen bond pattern |
//! | G | 3₁₀-helix | i → i+3 hydrogen bond pattern |
//! | I | π-helix | i → i+5 hydrogen bond pattern |
//! | P | κ-helix | PPII helix (polyproline II), dihedral-based |
//! | E | Extended strand | Part of β-sheet |
//! | B | Beta bridge | Isolated β-bridge |
//! | T | Turn | Hydrogen-bonded turn |
//! | S | Bend | High backbone curvature |
//! | C | Coil | None of the above |
//!
//! # Usage
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Get complete assignment
//! let ss = structure.assign_secondary_structure();
//! println!("Secondary structure: {}", ss);
//!
//! // Get as string of codes (e.g., "HHHHEEEECCCC")
//! let ss_string = structure.secondary_structure_string();
//!
//! // Get composition fractions
//! let (helix, sheet, coil) = structure.secondary_structure_composition();
//! println!("Helix: {:.1}%, Sheet: {:.1}%, Coil: {:.1}%",
//!     helix * 100.0, sheet * 100.0, coil * 100.0);
//! ```
//!
//! # Algorithm Details
//!
//! ## Hydrogen Bond Detection
//!
//! H-bonds are detected using the Kabsch-Sander electrostatic model:
//!
//! ```text
//! E = 27.888 × (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
//! ```
//!
//! A hydrogen bond is assigned when E < -0.5 kcal/mol.
//!
//! ## Helix Detection
//!
//! Helices are detected based on repeating H-bond patterns:
//! - α-helix: consecutive i → i+4 H-bonds
//! - 3₁₀-helix: consecutive i → i+3 H-bonds
//! - π-helix: consecutive i → i+5 H-bonds
//!
//! ## β-Sheet Detection
//!
//! β-sheets are detected from H-bond patterns between non-consecutive residues:
//! - Parallel: alternating H-bond directions
//! - Antiparallel: symmetric H-bond pairs
//!
//! ## PPII/κ-Helix Detection (DSSP 4)
//!
//! PPII helices are detected in coil regions based on backbone dihedrals:
//! - φ = -75° ± 29°
//! - ψ = +145° ± 29°
//! - Minimum 2 consecutive residues
//!
//! # Feature Flag
//!
//! This module requires the `dssp` feature:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.5", features = ["dssp"] }
//! ```
//!
//! # References
//!
//! - Kabsch W, Sander C. Dictionary of protein secondary structure:
//!   Pattern recognition of hydrogen-bonded and geometrical features.
//!   Biopolymers. 1983;22(12):2577-2637.
//! - Hekkelman ML et al. DSSP 4: FAIR annotation of protein secondary structure.
//!   Protein Sci. 2025;34(8):e70208.

mod assignment;
mod dihedral;
mod hbond;
mod patterns;
mod types;

// Re-export public types
pub use types::{ResidueSSAssignment, SecondaryStructure, SecondaryStructureAssignment};

// Re-export constants that users might need
pub use dihedral::{PPII_ANGLE_TOLERANCE, PPII_PHI_CENTER, PPII_PSI_CENTER};
pub use hbond::HBOND_THRESHOLD;
pub use patterns::BEND_ANGLE_THRESHOLD;

// Re-export dihedral calculation types and functions for the dihedrals API
pub use dihedral::{
    BackboneDihedrals, PPII_MIN_CONSECUTIVE, calculate_all_dihedrals, calculate_dihedral,
    calculate_omega, calculate_phi, calculate_psi,
};

// Re-export hydrogen bond types and functions for the hbonds API
pub use hbond::{
    BackboneAtoms, HydrogenBond, calculate_hbond_energy, compute_all_virtual_hydrogens,
    detect_hydrogen_bonds, extract_backbone_atoms,
};

use crate::PdbStructure;

/// Extension trait for secondary structure assignment on PdbStructure.
impl PdbStructure {
    /// Computes DSSP-like secondary structure assignment.
    ///
    /// This method analyzes the protein structure and assigns secondary
    /// structure to each residue using the DSSP algorithm.
    ///
    /// # Returns
    ///
    /// A `SecondaryStructureAssignment` containing:
    /// - Per-residue assignments
    /// - Statistics (helix/sheet/coil counts and fractions)
    /// - Any warnings generated during assignment
    ///
    /// # Example
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let ss = structure.assign_secondary_structure();
    ///
    /// println!("Total residues: {}", ss.len());
    /// println!("Helix fraction: {:.1}%", ss.helix_fraction * 100.0);
    /// println!("Sheet fraction: {:.1}%", ss.sheet_fraction * 100.0);
    ///
    /// // Iterate over assignments
    /// for res in &ss.residue_assignments {
    ///     println!("{}{}: {}", res.chain_id, res.residue_seq, res.ss.code());
    /// }
    /// ```
    pub fn assign_secondary_structure(&self) -> SecondaryStructureAssignment {
        assignment::assign_secondary_structure(&self.atoms)
    }

    /// Returns the secondary structure as a string of single-character codes.
    ///
    /// The string contains one character per residue in sequence order,
    /// using standard DSSP codes (H, G, I, P, E, B, T, S, C).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let ss_string = structure.secondary_structure_string();
    /// println!("Secondary structure: {}", ss_string);
    /// // Output: "CCCHHHHHHHHHCCCEEEEECCC"
    /// ```
    pub fn secondary_structure_string(&self) -> String {
        self.assign_secondary_structure().to_string()
    }

    /// Returns the secondary structure composition as fractions.
    ///
    /// # Returns
    ///
    /// A tuple of (helix_fraction, sheet_fraction, coil_fraction) where:
    /// - helix includes H, G, I, and P (κ-helix/PPII)
    /// - sheet includes E and B
    /// - coil includes T, S, and C
    ///
    /// All fractions sum to 1.0 (within floating-point precision).
    ///
    /// # Example
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let (helix, sheet, coil) = structure.secondary_structure_composition();
    ///
    /// println!("Helix: {:.1}%", helix * 100.0);
    /// println!("Sheet: {:.1}%", sheet * 100.0);
    /// println!("Coil:  {:.1}%", coil * 100.0);
    /// ```
    pub fn secondary_structure_composition(&self) -> (f64, f64, f64) {
        self.assign_secondary_structure().composition()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_minimal_residue(chain_id: &str, seq: i32, res_name: &str) -> Vec<Atom> {
        // Create a minimal residue with backbone atoms
        // This is not geometrically accurate but tests the API
        vec![
            Atom {
                serial: seq * 4,
                name: "N".to_string(),
                alt_loc: None,
                residue_name: res_name.to_string(),
                chain_id: chain_id.to_string(),
                residue_seq: seq,
                ins_code: None,
                x: (seq as f64) * 3.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "N".to_string(),
            },
            Atom {
                serial: seq * 4 + 1,
                name: "CA".to_string(),
                alt_loc: None,
                residue_name: res_name.to_string(),
                chain_id: chain_id.to_string(),
                residue_seq: seq,
                ins_code: None,
                x: (seq as f64) * 3.0 + 1.5,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            Atom {
                serial: seq * 4 + 2,
                name: "C".to_string(),
                alt_loc: None,
                residue_name: res_name.to_string(),
                chain_id: chain_id.to_string(),
                residue_seq: seq,
                ins_code: None,
                x: (seq as f64) * 3.0 + 2.5,
                y: 1.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            Atom {
                serial: seq * 4 + 3,
                name: "O".to_string(),
                alt_loc: None,
                residue_name: res_name.to_string(),
                chain_id: chain_id.to_string(),
                residue_seq: seq,
                ins_code: None,
                x: (seq as f64) * 3.0 + 2.5,
                y: 2.2,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "O".to_string(),
            },
        ]
    }

    #[test]
    fn test_empty_structure() {
        let structure = PdbStructure::new();
        let ss = structure.assign_secondary_structure();
        assert!(ss.is_empty());
    }

    #[test]
    fn test_secondary_structure_string() {
        let mut structure = PdbStructure::new();
        for i in 1..=5 {
            structure
                .atoms
                .extend(create_minimal_residue("A", i, "ALA"));
        }

        let ss_string = structure.secondary_structure_string();
        assert_eq!(ss_string.len(), 5);
        // All characters should be valid DSSP codes
        for c in ss_string.chars() {
            assert!("HGIPEBTSC".contains(c), "Invalid SS code: {}", c);
        }
    }

    #[test]
    fn test_secondary_structure_composition() {
        let mut structure = PdbStructure::new();
        for i in 1..=10 {
            structure
                .atoms
                .extend(create_minimal_residue("A", i, "ALA"));
        }

        let (helix, sheet, coil) = structure.secondary_structure_composition();

        // Fractions should sum to 1.0
        assert!((helix + sheet + coil - 1.0).abs() < 0.01);
        // All fractions should be between 0 and 1
        assert!((0.0..=1.0).contains(&helix));
        assert!((0.0..=1.0).contains(&sheet));
        assert!((0.0..=1.0).contains(&coil));
    }

    #[test]
    fn test_secondary_structure_codes() {
        // Test that all codes are recognized
        assert_eq!(SecondaryStructure::AlphaHelix.code(), 'H');
        assert_eq!(SecondaryStructure::Helix310.code(), 'G');
        assert_eq!(SecondaryStructure::PiHelix.code(), 'I');
        assert_eq!(SecondaryStructure::KappaHelix.code(), 'P');
        assert_eq!(SecondaryStructure::ExtendedStrand.code(), 'E');
        assert_eq!(SecondaryStructure::BetaBridge.code(), 'B');
        assert_eq!(SecondaryStructure::Turn.code(), 'T');
        assert_eq!(SecondaryStructure::Bend.code(), 'S');
        assert_eq!(SecondaryStructure::Coil.code(), 'C');
    }
}
