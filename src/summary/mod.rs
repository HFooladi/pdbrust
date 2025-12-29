//! Unified structure summary functionality.
//!
//! This module provides a single-call interface for generating comprehensive
//! structural summaries that combine quality assessment and structural descriptors.
//! This is useful for batch processing, dataset characterization, and generating
//! tabular outputs for downstream analysis.
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Get comprehensive summary
//! let summary = structure.summary();
//!
//! println!("Atoms: {}", summary.num_atoms);
//! println!("Residues: {}", summary.num_residues);
//! println!("Rg: {:.2} Å", summary.radius_of_gyration);
//! println!("Hydrophobic: {:.1}%", summary.hydrophobic_ratio * 100.0);
//!
//! // Check if suitable for analysis
//! if summary.is_analysis_ready() {
//!     println!("Structure is ready for analysis");
//! }
//! ```
//!
//! # Feature Flag
//!
//! This module requires the `summary` feature, which automatically enables
//! the `descriptors` and `quality` features:
//!
//! ```toml
//! [dependencies]
//! pdbrust = { version = "0.1", features = ["summary"] }
//! ```

use crate::core::PdbStructure;
use crate::descriptors::StructureDescriptors;
use crate::quality::QualityReport;
use std::collections::HashMap;

/// Comprehensive structure summary combining quality and descriptors.
///
/// This struct provides a unified view of all computed metrics for a
/// protein structure, suitable for:
/// - Batch processing of multiple structures
/// - Generating CSV/tabular outputs
/// - Filtering datasets based on quality criteria
/// - Statistical analysis and visualization
#[derive(Debug, Clone)]
pub struct StructureSummary {
    // ========== Quality Indicators ==========
    /// True if the structure contains only Cα atoms
    pub has_ca_only: bool,
    /// True if multiple models are present (NMR ensemble)
    pub has_multiple_models: bool,
    /// True if alternate locations are present
    pub has_altlocs: bool,
    /// Number of distinct chains
    pub num_chains: usize,
    /// Number of models
    pub num_models: usize,
    /// True if HETATM records are present
    pub has_hetatm: bool,
    /// True if hydrogen atoms are present
    pub has_hydrogens: bool,
    /// True if disulfide bonds are defined
    pub has_ssbonds: bool,

    // ========== Size Descriptors ==========
    /// Number of residues (based on Cα count)
    pub num_residues: usize,
    /// Total number of atoms
    pub num_atoms: usize,
    /// Fraction of missing residues (0.0 = complete)
    pub missing_residue_ratio: f64,

    // ========== Composition ==========
    /// Fraction of glycine residues
    pub glycine_ratio: f64,
    /// Fraction of hydrophobic residues
    pub hydrophobic_ratio: f64,
    /// Full amino acid composition
    pub aa_composition: HashMap<String, f64>,

    // ========== Geometric Descriptors ==========
    /// Radius of gyration (Å)
    pub radius_of_gyration: f64,
    /// Maximum Cα-Cα distance (Å)
    pub max_ca_distance: f64,
    /// Heuristic secondary structure content
    pub secondary_structure_ratio: f64,
    /// Compactness index: Rg / n^(1/3)
    pub compactness_index: f64,
    /// Cα density: count / bounding box volume
    pub ca_density: f64,
}

impl Default for StructureSummary {
    fn default() -> Self {
        Self {
            has_ca_only: false,
            has_multiple_models: false,
            has_altlocs: false,
            num_chains: 0,
            num_models: 0,
            has_hetatm: false,
            has_hydrogens: false,
            has_ssbonds: false,
            num_residues: 0,
            num_atoms: 0,
            missing_residue_ratio: 0.0,
            glycine_ratio: 0.0,
            hydrophobic_ratio: 0.0,
            aa_composition: HashMap::new(),
            radius_of_gyration: 0.0,
            max_ca_distance: 0.0,
            secondary_structure_ratio: 0.0,
            compactness_index: 0.0,
            ca_density: 0.0,
        }
    }
}

impl StructureSummary {
    /// Check if the structure is suitable for typical analysis.
    ///
    /// Returns true if the structure:
    /// - Has atoms
    /// - Is a single model (not NMR ensemble)
    /// - Has no alternate locations
    /// - Is not CA-only
    pub fn is_analysis_ready(&self) -> bool {
        self.num_atoms > 0
            && !self.has_multiple_models
            && !self.has_altlocs
            && !self.has_ca_only
    }

    /// Check if the structure passes basic quality criteria.
    pub fn is_clean(&self) -> bool {
        self.num_atoms > 0 && !self.has_ca_only && !self.has_altlocs
    }

    /// Get a vector of field names for CSV header generation.
    pub fn field_names() -> Vec<&'static str> {
        vec![
            "has_ca_only",
            "has_multiple_models",
            "has_altlocs",
            "num_chains",
            "num_models",
            "has_hetatm",
            "has_hydrogens",
            "has_ssbonds",
            "num_residues",
            "num_atoms",
            "missing_residue_ratio",
            "glycine_ratio",
            "hydrophobic_ratio",
            "radius_of_gyration",
            "max_ca_distance",
            "secondary_structure_ratio",
            "compactness_index",
            "ca_density",
        ]
    }

    /// Convert summary to a vector of string values for CSV output.
    ///
    /// Values are ordered to match `field_names()`.
    pub fn to_csv_values(&self) -> Vec<String> {
        vec![
            self.has_ca_only.to_string(),
            self.has_multiple_models.to_string(),
            self.has_altlocs.to_string(),
            self.num_chains.to_string(),
            self.num_models.to_string(),
            self.has_hetatm.to_string(),
            self.has_hydrogens.to_string(),
            self.has_ssbonds.to_string(),
            self.num_residues.to_string(),
            self.num_atoms.to_string(),
            format!("{:.6}", self.missing_residue_ratio),
            format!("{:.6}", self.glycine_ratio),
            format!("{:.6}", self.hydrophobic_ratio),
            format!("{:.4}", self.radius_of_gyration),
            format!("{:.4}", self.max_ca_distance),
            format!("{:.6}", self.secondary_structure_ratio),
            format!("{:.6}", self.compactness_index),
            format!("{:.8}", self.ca_density),
        ]
    }

    /// Create a summary from quality report and descriptors.
    pub fn from_parts(quality: QualityReport, descriptors: StructureDescriptors) -> Self {
        Self {
            has_ca_only: quality.has_ca_only,
            has_multiple_models: quality.has_multiple_models,
            has_altlocs: quality.has_altlocs,
            num_chains: quality.num_chains,
            num_models: quality.num_models,
            has_hetatm: quality.has_hetatm,
            has_hydrogens: quality.has_hydrogens,
            has_ssbonds: quality.has_ssbonds,
            num_residues: descriptors.num_residues,
            num_atoms: descriptors.num_atoms,
            missing_residue_ratio: descriptors.missing_residue_ratio,
            glycine_ratio: descriptors.glycine_ratio,
            hydrophobic_ratio: descriptors.hydrophobic_ratio,
            aa_composition: descriptors.aa_composition,
            radius_of_gyration: descriptors.radius_of_gyration,
            max_ca_distance: descriptors.max_ca_distance,
            secondary_structure_ratio: descriptors.secondary_structure_ratio,
            compactness_index: descriptors.compactness_index,
            ca_density: descriptors.ca_density,
        }
    }
}

impl PdbStructure {
    /// Generate a comprehensive structure summary.
    ///
    /// This method combines quality assessment and structural descriptors
    /// into a single unified summary object. It's more efficient than
    /// calling `quality_report()` and `structure_descriptors()` separately
    /// when you need both.
    ///
    /// # Returns
    ///
    /// A `StructureSummary` containing all quality metrics and descriptors.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let summary = structure.summary();
    ///
    /// // Access quality metrics
    /// println!("Chains: {}", summary.num_chains);
    /// println!("Has altlocs: {}", summary.has_altlocs);
    ///
    /// // Access structural descriptors
    /// println!("Rg: {:.2} Å", summary.radius_of_gyration);
    /// println!("Compactness: {:.2}", summary.compactness_index);
    ///
    /// // Check suitability for analysis
    /// if summary.is_analysis_ready() {
    ///     println!("Ready for analysis!");
    /// }
    /// ```
    pub fn summary(&self) -> StructureSummary {
        let quality = self.quality_report();
        let descriptors = self.structure_descriptors();

        StructureSummary::from_parts(quality, descriptors)
    }
}

/// Generate summaries for multiple structures.
///
/// This is a convenience function for batch processing.
///
/// # Arguments
///
/// * `structures` - A slice of PdbStructure references
///
/// # Returns
///
/// A vector of StructureSummary, one for each input structure.
///
/// # Examples
///
/// ```ignore
/// use pdbrust::{PdbStructure, parse_pdb_file};
/// use pdbrust::summary::batch_summarize;
///
/// let structures: Vec<PdbStructure> = vec![
///     parse_pdb_file("protein1.pdb")?,
///     parse_pdb_file("protein2.pdb")?,
/// ];
///
/// let summaries = batch_summarize(&structures);
/// ```
pub fn batch_summarize(structures: &[PdbStructure]) -> Vec<StructureSummary> {
    structures.iter().map(|s| s.summary()).collect()
}

/// Generate a CSV string from multiple structure summaries.
///
/// # Arguments
///
/// * `summaries` - A slice of StructureSummary
/// * `include_header` - Whether to include a header row
///
/// # Returns
///
/// A CSV-formatted string.
pub fn summaries_to_csv(summaries: &[StructureSummary], include_header: bool) -> String {
    let mut output = String::new();

    if include_header {
        output.push_str(&StructureSummary::field_names().join(","));
        output.push('\n');
    }

    for summary in summaries {
        output.push_str(&summary.to_csv_values().join(","));
        output.push('\n');
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            create_atom(1, " N  ", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_atom(2, " CA ", "ALA", "A", 1, 1.5, 0.0, 0.0),
            create_atom(3, " C  ", "ALA", "A", 1, 3.0, 0.0, 0.0),
            create_atom(4, " CA ", "GLY", "A", 2, 6.8, 0.0, 0.0),
            create_atom(5, " CA ", "VAL", "A", 3, 10.6, 0.0, 0.0),
        ];

        structure
    }

    #[allow(clippy::too_many_arguments)]
    fn create_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_summary() {
        let structure = create_test_structure();
        let summary = structure.summary();

        assert_eq!(summary.num_atoms, 5);
        assert_eq!(summary.num_chains, 1);
        assert!(!summary.has_ca_only);
        assert!(!summary.has_multiple_models);
        assert!(!summary.has_altlocs);
    }

    #[test]
    fn test_summary_descriptors() {
        let structure = create_test_structure();
        let summary = structure.summary();

        // Should have 3 CA atoms (residues)
        assert_eq!(summary.num_residues, 3);
        assert!(summary.radius_of_gyration > 0.0);
        assert!(summary.max_ca_distance > 0.0);
    }

    #[test]
    fn test_summary_is_analysis_ready() {
        let structure = create_test_structure();
        let summary = structure.summary();

        // Simple structure should be analysis-ready
        assert!(summary.is_analysis_ready());
        assert!(summary.is_clean());
    }

    #[test]
    fn test_summary_empty_structure() {
        let structure = PdbStructure::new();
        let summary = structure.summary();

        assert_eq!(summary.num_atoms, 0);
        assert_eq!(summary.num_residues, 0);
        assert!(!summary.is_analysis_ready());
    }

    #[test]
    fn test_summary_default() {
        let summary = StructureSummary::default();

        assert_eq!(summary.num_atoms, 0);
        assert!(!summary.has_ca_only);
        assert!(summary.aa_composition.is_empty());
    }

    #[test]
    fn test_field_names() {
        let names = StructureSummary::field_names();

        assert!(names.contains(&"num_atoms"));
        assert!(names.contains(&"radius_of_gyration"));
        assert!(names.contains(&"has_altlocs"));
    }

    #[test]
    fn test_to_csv_values() {
        let structure = create_test_structure();
        let summary = structure.summary();
        let values = summary.to_csv_values();

        // Should have same length as field names
        assert_eq!(values.len(), StructureSummary::field_names().len());
    }

    #[test]
    fn test_batch_summarize() {
        let structures = vec![
            create_test_structure(),
            create_test_structure(),
        ];

        let summaries = batch_summarize(&structures);

        assert_eq!(summaries.len(), 2);
        assert_eq!(summaries[0].num_atoms, summaries[1].num_atoms);
    }

    #[test]
    fn test_summaries_to_csv() {
        let structures = vec![create_test_structure()];
        let summaries = batch_summarize(&structures);

        let csv_with_header = summaries_to_csv(&summaries, true);
        let csv_without_header = summaries_to_csv(&summaries, false);

        // With header should have more content
        assert!(csv_with_header.len() > csv_without_header.len());

        // Should contain field names in header
        assert!(csv_with_header.contains("num_atoms"));
        assert!(csv_with_header.contains("radius_of_gyration"));
    }

    #[test]
    fn test_from_parts() {
        let quality = QualityReport {
            has_ca_only: false,
            has_multiple_models: true,
            has_altlocs: false,
            num_chains: 2,
            num_models: 5,
            num_atoms: 100,
            num_residues: 10,
            has_hetatm: true,
            has_hydrogens: false,
            has_ssbonds: true,
            has_conect: false,
        };

        let descriptors = StructureDescriptors {
            num_residues: 50,
            num_atoms: 500,
            aa_composition: HashMap::new(),
            glycine_ratio: 0.1,
            hydrophobic_ratio: 0.4,
            radius_of_gyration: 15.0,
            max_ca_distance: 40.0,
            missing_residue_ratio: 0.05,
            secondary_structure_ratio: 0.8,
            compactness_index: 2.5,
            ca_density: 0.005,
        };

        let summary = StructureSummary::from_parts(quality, descriptors);

        // Quality fields
        assert!(summary.has_multiple_models);
        assert_eq!(summary.num_models, 5);
        assert!(summary.has_ssbonds);

        // Descriptor fields
        assert_eq!(summary.num_residues, 50);
        assert!((summary.radius_of_gyration - 15.0).abs() < 1e-10);
    }
}
