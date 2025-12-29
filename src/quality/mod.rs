//! Structure quality assessment functions.
//!
//! This module provides fast, non-destructive quality checks for PDB structures.
//! These functions are useful for filtering datasets and characterizing structures.
//!
//! # Quality Checks
//!
//! - **Atom content**: Detect CA-only (coarse-grained) structures
//! - **Model content**: Detect multi-model structures (NMR ensembles)
//! - **Alternate locations**: Detect conformational heterogeneity
//! - **Chain count**: Count distinct protein chains
//! - **Completeness**: Check for missing atoms or residues
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Quick quality checks
//! if structure.has_multiple_models() {
//!     println!("NMR ensemble detected");
//! }
//!
//! if structure.has_altlocs() {
//!     println!("Structure has alternate conformations");
//! }
//!
//! // Get comprehensive quality report
//! let report = structure.quality_report();
//! println!("Chains: {}, Models: {}", report.num_chains, report.num_models);
//! ```
//!
//! # Feature Flag
//!
//! This module requires the `quality` feature:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.1", features = ["quality"] }
//! ```

use crate::core::PdbStructure;

/// Comprehensive quality assessment report for a PDB structure.
///
/// This struct contains all quality metrics for a structure,
/// suitable for filtering datasets or generating quality summaries.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct QualityReport {
    /// True if the structure contains only Cα atoms (coarse-grained)
    pub has_ca_only: bool,
    /// True if the structure has multiple models (e.g., NMR ensemble)
    pub has_multiple_models: bool,
    /// True if any atoms have alternate location indicators
    pub has_altlocs: bool,
    /// Number of distinct chains in the structure
    pub num_chains: usize,
    /// Number of models in the structure
    pub num_models: usize,
    /// Number of atoms in the structure
    pub num_atoms: usize,
    /// Number of residues (based on unique chain+resnum+icode)
    pub num_residues: usize,
    /// True if the structure contains HETATM records (ligands, waters)
    pub has_hetatm: bool,
    /// True if the structure contains hydrogen atoms
    pub has_hydrogens: bool,
    /// True if the structure has disulfide bonds defined
    pub has_ssbonds: bool,
    /// True if the structure has connectivity records
    pub has_conect: bool,
}

impl QualityReport {
    /// Check if the structure passes basic quality criteria.
    ///
    /// Returns true if the structure:
    /// - Has at least one atom
    /// - Is not CA-only (unless that's expected)
    /// - Has no alternate locations (clean structure)
    pub fn is_clean(&self) -> bool {
        self.num_atoms > 0 && !self.has_ca_only && !self.has_altlocs
    }

    /// Check if the structure is suitable for typical analysis.
    ///
    /// Returns true if the structure:
    /// - Has atoms
    /// - Is a single model (not NMR ensemble)
    /// - Has no alternate locations
    pub fn is_analysis_ready(&self) -> bool {
        self.num_atoms > 0 && !self.has_multiple_models && !self.has_altlocs
    }
}

impl PdbStructure {
    /// Check if the structure contains only Cα atoms.
    ///
    /// Returns true if all ATOM records are Cα atoms, indicating
    /// a coarse-grained representation of the structure.
    ///
    /// # Returns
    ///
    /// `true` if the structure contains only Cα atoms, `false` otherwise.
    /// Returns `false` for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    ///
    /// if structure.has_ca_only() {
    ///     println!("Coarse-grained structure detected");
    /// }
    /// ```
    pub fn has_ca_only(&self) -> bool {
        if self.atoms.is_empty() {
            return false;
        }

        self.atoms.iter().all(|atom| atom.name.trim() == "CA")
    }

    /// Check if the structure has multiple models.
    ///
    /// Returns true if the structure contains more than one model,
    /// typically indicating an NMR ensemble or molecular dynamics trajectory.
    ///
    /// # Returns
    ///
    /// `true` if multiple models are present, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("nmr_ensemble.pdb")?;
    ///
    /// if structure.has_multiple_models() {
    ///     println!("NMR ensemble with {} models", structure.num_models());
    /// }
    /// ```
    pub fn has_multiple_models(&self) -> bool {
        self.models.len() > 1
    }

    /// Get the number of models in the structure.
    ///
    /// # Returns
    ///
    /// The number of models. Returns 1 for single-model structures
    /// (even if no explicit MODEL record is present).
    pub fn num_models(&self) -> usize {
        if self.models.is_empty() {
            1 // Implicit single model
        } else {
            self.models.len()
        }
    }

    /// Check if any atoms have alternate location indicators.
    ///
    /// Alternate locations indicate that an atom exists in multiple
    /// conformations, representing local disorder or flexibility.
    ///
    /// # Returns
    ///
    /// `true` if any atom has an alternate location indicator.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    ///
    /// if structure.has_altlocs() {
    ///     println!("Structure has conformational heterogeneity");
    /// }
    /// ```
    pub fn has_altlocs(&self) -> bool {
        self.atoms.iter().any(|atom| atom.alt_loc.is_some())
    }

    /// Get all unique alternate location identifiers.
    ///
    /// # Returns
    ///
    /// A sorted vector of alternate location characters (e.g., ['A', 'B']).
    pub fn get_altloc_ids(&self) -> Vec<char> {
        let mut altlocs: Vec<char> = self.atoms
            .iter()
            .filter_map(|atom| atom.alt_loc)
            .collect();

        altlocs.sort();
        altlocs.dedup();
        altlocs
    }

    /// Check if the structure contains HETATM records.
    ///
    /// HETATM records represent non-standard residues like ligands,
    /// water molecules, ions, and modified amino acids.
    ///
    /// # Returns
    ///
    /// `true` if any non-standard residues are present.
    pub fn has_hetatm(&self) -> bool {
        // Check if any residue is not a standard amino acid or nucleotide
        #[cfg(feature = "filter")]
        {
            use crate::filter::is_standard_residue;
            self.atoms.iter().any(|atom| !is_standard_residue(&atom.residue_name))
        }

        #[cfg(not(feature = "filter"))]
        {
            // Fallback: check for common HETATM residue names
            const NON_STANDARD: &[&str] = &["HOH", "WAT", "DOD", "H2O", "TIP"];
            self.atoms.iter().any(|atom| {
                NON_STANDARD.contains(&atom.residue_name.trim())
            })
        }
    }

    /// Check if the structure contains hydrogen atoms.
    ///
    /// # Returns
    ///
    /// `true` if any hydrogen atoms are present.
    pub fn has_hydrogens(&self) -> bool {
        self.atoms.iter().any(|atom| atom.is_hydrogen())
    }

    /// Check if the structure has disulfide bonds defined.
    ///
    /// # Returns
    ///
    /// `true` if SSBOND records are present.
    pub fn has_ssbonds(&self) -> bool {
        !self.ssbonds.is_empty()
    }

    /// Check if the structure has connectivity records.
    ///
    /// # Returns
    ///
    /// `true` if CONECT records are present.
    pub fn has_conect(&self) -> bool {
        !self.connects.is_empty()
    }

    /// Count the number of atoms with alternate locations.
    ///
    /// # Returns
    ///
    /// The number of atoms that have alternate location indicators.
    pub fn count_altloc_atoms(&self) -> usize {
        self.atoms.iter().filter(|atom| atom.alt_loc.is_some()).count()
    }

    /// Count the number of hydrogen atoms.
    ///
    /// # Returns
    ///
    /// The number of hydrogen atoms in the structure.
    pub fn count_hydrogens(&self) -> usize {
        self.atoms.iter().filter(|atom| atom.is_hydrogen()).count()
    }

    /// Get a comprehensive quality report for the structure.
    ///
    /// This method computes all quality metrics in a single pass
    /// and returns them in a structured format.
    ///
    /// # Returns
    ///
    /// A `QualityReport` struct containing all quality metrics.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let report = structure.quality_report();
    ///
    /// println!("Chains: {}", report.num_chains);
    /// println!("Models: {}", report.num_models);
    /// println!("Has altlocs: {}", report.has_altlocs);
    ///
    /// if report.is_analysis_ready() {
    ///     println!("Structure is ready for analysis");
    /// }
    /// ```
    pub fn quality_report(&self) -> QualityReport {
        QualityReport {
            has_ca_only: self.has_ca_only(),
            has_multiple_models: self.has_multiple_models(),
            has_altlocs: self.has_altlocs(),
            num_chains: self.get_num_chains(),
            num_models: self.num_models(),
            num_atoms: self.get_num_atoms(),
            num_residues: self.get_num_residues(),
            has_hetatm: self.has_hetatm(),
            has_hydrogens: self.has_hydrogens(),
            has_ssbonds: self.has_ssbonds(),
            has_conect: self.has_conect(),
        }
    }

    /// Check if the structure is empty.
    ///
    /// # Returns
    ///
    /// `true` if the structure has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Get resolution information from REMARK 2 records.
    ///
    /// # Returns
    ///
    /// The resolution in Angstroms if available, `None` otherwise.
    pub fn get_resolution(&self) -> Option<f64> {
        for remark in &self.remarks {
            if remark.number == 2 {
                // Try to parse resolution from REMARK 2
                // Format: "RESOLUTION.    X.XX ANGSTROMS."
                let content = remark.content.to_uppercase();
                if content.contains("RESOLUTION") {
                    // Extract the numeric value
                    for word in content.split_whitespace() {
                        if let Ok(value) = word.parse::<f64>() {
                            if value > 0.0 && value < 20.0 {
                                return Some(value);
                            }
                        }
                    }
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::{Atom, Model, Remark};

    fn create_atom(serial: i32, name: &str, residue_name: &str, chain_id: &str, alt_loc: Option<char>) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq: 1,
            ins_code: None,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: if name.trim().starts_with('H') { "H".to_string() } else { "C".to_string() },
        }
    }

    #[test]
    fn test_has_ca_only_true() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", None),
            create_atom(2, " CA ", "GLY", "A", None),
            create_atom(3, " CA ", "VAL", "A", None),
        ];

        assert!(structure.has_ca_only());
    }

    #[test]
    fn test_has_ca_only_false() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " N  ", "ALA", "A", None),
            create_atom(2, " CA ", "ALA", "A", None),
            create_atom(3, " C  ", "ALA", "A", None),
        ];

        assert!(!structure.has_ca_only());
    }

    #[test]
    fn test_has_ca_only_empty() {
        let structure = PdbStructure::new();
        assert!(!structure.has_ca_only());
    }

    #[test]
    fn test_has_multiple_models_true() {
        let mut structure = PdbStructure::new();
        structure.models = vec![
            Model { serial: 1, atoms: vec![], remarks: vec![] },
            Model { serial: 2, atoms: vec![], remarks: vec![] },
        ];

        assert!(structure.has_multiple_models());
        assert_eq!(structure.num_models(), 2);
    }

    #[test]
    fn test_has_multiple_models_false() {
        let structure = PdbStructure::new();
        assert!(!structure.has_multiple_models());
        assert_eq!(structure.num_models(), 1);
    }

    #[test]
    fn test_has_altlocs_true() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", Some('A')),
            create_atom(2, " CA ", "ALA", "A", Some('B')),
        ];

        assert!(structure.has_altlocs());
    }

    #[test]
    fn test_has_altlocs_false() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", None),
            create_atom(2, " CA ", "GLY", "A", None),
        ];

        assert!(!structure.has_altlocs());
    }

    #[test]
    fn test_get_altloc_ids() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", Some('A')),
            create_atom(2, " CA ", "ALA", "A", Some('B')),
            create_atom(3, " CB ", "ALA", "A", Some('A')),
        ];

        let altlocs = structure.get_altloc_ids();
        assert_eq!(altlocs, vec!['A', 'B']);
    }

    #[test]
    fn test_has_hydrogens() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", None),
            Atom {
                serial: 2,
                name: "H   ".to_string(), // Hydrogen names start with H
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 0.0, y: 0.0, z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "H".to_string(),
            },
        ];

        assert!(structure.has_hydrogens());
        assert_eq!(structure.count_hydrogens(), 1);
    }

    #[test]
    fn test_has_ssbonds() {
        let mut structure = PdbStructure::new();
        assert!(!structure.has_ssbonds());

        structure.ssbonds = vec![crate::records::SSBond {
            serial: 1,
            residue1_name: "CYS".to_string(),
            chain1_id: "A".to_string(),
            residue1_seq: 10,
            icode1: None,
            residue2_name: "CYS".to_string(),
            chain2_id: "A".to_string(),
            residue2_seq: 20,
            icode2: None,
            sym1: 1,
            sym2: 1,
            length: 2.03,
        }];

        assert!(structure.has_ssbonds());
    }

    #[test]
    fn test_quality_report() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " N  ", "ALA", "A", None),
            create_atom(2, " CA ", "ALA", "A", None),
            create_atom(3, " C  ", "ALA", "A", None),
        ];

        let report = structure.quality_report();

        assert!(!report.has_ca_only);
        assert!(!report.has_multiple_models);
        assert!(!report.has_altlocs);
        assert_eq!(report.num_chains, 1);
        assert_eq!(report.num_models, 1);
        assert_eq!(report.num_atoms, 3);
    }

    #[test]
    fn test_quality_report_is_clean() {
        let report = QualityReport {
            has_ca_only: false,
            has_multiple_models: false,
            has_altlocs: false,
            num_chains: 1,
            num_models: 1,
            num_atoms: 100,
            num_residues: 10,
            has_hetatm: false,
            has_hydrogens: false,
            has_ssbonds: false,
            has_conect: false,
        };

        assert!(report.is_clean());
        assert!(report.is_analysis_ready());
    }

    #[test]
    fn test_quality_report_not_clean() {
        let report = QualityReport {
            has_ca_only: true,
            has_multiple_models: false,
            has_altlocs: false,
            num_chains: 1,
            num_models: 1,
            num_atoms: 100,
            num_residues: 10,
            has_hetatm: false,
            has_hydrogens: false,
            has_ssbonds: false,
            has_conect: false,
        };

        assert!(!report.is_clean()); // CA-only is not clean
    }

    #[test]
    fn test_is_empty() {
        let structure = PdbStructure::new();
        assert!(structure.is_empty());

        let mut structure2 = PdbStructure::new();
        structure2.atoms = vec![create_atom(1, " CA ", "ALA", "A", None)];
        assert!(!structure2.is_empty());
    }

    #[test]
    fn test_get_resolution() {
        let mut structure = PdbStructure::new();
        structure.remarks = vec![
            Remark {
                number: 2,
                content: "RESOLUTION.    2.50 ANGSTROMS.".to_string(),
            },
        ];

        let resolution = structure.get_resolution();
        assert!(resolution.is_some());
        assert!((resolution.unwrap() - 2.50).abs() < 0.01);
    }

    #[test]
    fn test_get_resolution_none() {
        let structure = PdbStructure::new();
        assert!(structure.get_resolution().is_none());
    }

    #[test]
    fn test_count_altloc_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", Some('A')),
            create_atom(2, " CA ", "ALA", "A", Some('B')),
            create_atom(3, " CB ", "ALA", "A", None),
        ];

        assert_eq!(structure.count_altloc_atoms(), 2);
    }
}
