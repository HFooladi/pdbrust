//! Python bindings for unified structure summaries

use pdbrust::summary::StructureSummary;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Unified structure summary combining quality and descriptors
///
/// Provides a comprehensive overview of a structure's characteristics
/// in a single object.
#[pyclass(name = "StructureSummary")]
#[derive(Clone)]
pub struct PyStructureSummary {
    inner: StructureSummary,
}

#[pymethods]
impl PyStructureSummary {
    // ==================== Quality Fields ====================

    /// True if structure appears to be CA-only
    #[getter]
    fn has_ca_only(&self) -> bool {
        self.inner.has_ca_only
    }

    /// True if structure has multiple models
    #[getter]
    fn has_multiple_models(&self) -> bool {
        self.inner.has_multiple_models
    }

    /// True if any atoms have alternate locations
    #[getter]
    fn has_altlocs(&self) -> bool {
        self.inner.has_altlocs
    }

    /// Number of chains
    #[getter]
    fn num_chains(&self) -> usize {
        self.inner.num_chains
    }

    /// Number of models
    #[getter]
    fn num_models(&self) -> usize {
        self.inner.num_models
    }

    /// True if HETATM records present
    #[getter]
    fn has_hetatm(&self) -> bool {
        self.inner.has_hetatm
    }

    /// True if hydrogens present
    #[getter]
    fn has_hydrogens(&self) -> bool {
        self.inner.has_hydrogens
    }

    /// True if disulfide bonds defined
    #[getter]
    fn has_ssbonds(&self) -> bool {
        self.inner.has_ssbonds
    }

    // ==================== Size Fields ====================

    /// Number of residues
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues
    }

    /// Number of atoms
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    /// Missing residue ratio
    #[getter]
    fn missing_residue_ratio(&self) -> f64 {
        self.inner.missing_residue_ratio
    }

    // ==================== Composition Fields ====================

    /// Glycine ratio
    #[getter]
    fn glycine_ratio(&self) -> f64 {
        self.inner.glycine_ratio
    }

    /// Hydrophobic residue ratio
    #[getter]
    fn hydrophobic_ratio(&self) -> f64 {
        self.inner.hydrophobic_ratio
    }

    /// Amino acid composition dictionary
    #[getter]
    fn aa_composition(&self) -> HashMap<String, f64> {
        self.inner.aa_composition.clone()
    }

    // ==================== Geometric Fields ====================

    /// Radius of gyration in Angstroms
    #[getter]
    fn radius_of_gyration(&self) -> f64 {
        self.inner.radius_of_gyration
    }

    /// Maximum CA-CA distance in Angstroms
    #[getter]
    fn max_ca_distance(&self) -> f64 {
        self.inner.max_ca_distance
    }

    /// Secondary structure ratio (heuristic)
    #[getter]
    fn secondary_structure_ratio(&self) -> f64 {
        self.inner.secondary_structure_ratio
    }

    /// Compactness index
    #[getter]
    fn compactness_index(&self) -> f64 {
        self.inner.compactness_index
    }

    /// CA density
    #[getter]
    fn ca_density(&self) -> f64 {
        self.inner.ca_density
    }

    // ==================== Methods ====================

    /// Check if structure is analysis-ready
    fn is_analysis_ready(&self) -> bool {
        self.inner.is_analysis_ready()
    }

    /// Check if structure is clean
    fn is_clean(&self) -> bool {
        self.inner.is_clean()
    }

    /// Get CSV-formatted values
    fn to_csv_values(&self) -> Vec<String> {
        self.inner.to_csv_values()
    }

    /// Get field names for CSV header
    #[staticmethod]
    fn field_names() -> Vec<&'static str> {
        StructureSummary::field_names()
    }

    fn __repr__(&self) -> String {
        format!(
            "StructureSummary(atoms={}, residues={}, Rg={:.2}, analysis_ready={})",
            self.inner.num_atoms,
            self.inner.num_residues,
            self.inner.radius_of_gyration,
            self.inner.is_analysis_ready()
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, PyObject> {
        Python::with_gil(|py| {
            let mut dict = HashMap::new();
            // Quality
            dict.insert("has_ca_only".to_string(), self.inner.has_ca_only.into_py(py));
            dict.insert("has_multiple_models".to_string(), self.inner.has_multiple_models.into_py(py));
            dict.insert("has_altlocs".to_string(), self.inner.has_altlocs.into_py(py));
            dict.insert("num_chains".to_string(), self.inner.num_chains.into_py(py));
            dict.insert("num_models".to_string(), self.inner.num_models.into_py(py));
            dict.insert("has_hetatm".to_string(), self.inner.has_hetatm.into_py(py));
            dict.insert("has_hydrogens".to_string(), self.inner.has_hydrogens.into_py(py));
            dict.insert("has_ssbonds".to_string(), self.inner.has_ssbonds.into_py(py));
            // Size
            dict.insert("num_residues".to_string(), self.inner.num_residues.into_py(py));
            dict.insert("num_atoms".to_string(), self.inner.num_atoms.into_py(py));
            dict.insert("missing_residue_ratio".to_string(), self.inner.missing_residue_ratio.into_py(py));
            // Composition
            dict.insert("glycine_ratio".to_string(), self.inner.glycine_ratio.into_py(py));
            dict.insert("hydrophobic_ratio".to_string(), self.inner.hydrophobic_ratio.into_py(py));
            // Geometry
            dict.insert("radius_of_gyration".to_string(), self.inner.radius_of_gyration.into_py(py));
            dict.insert("max_ca_distance".to_string(), self.inner.max_ca_distance.into_py(py));
            dict.insert("secondary_structure_ratio".to_string(), self.inner.secondary_structure_ratio.into_py(py));
            dict.insert("compactness_index".to_string(), self.inner.compactness_index.into_py(py));
            dict.insert("ca_density".to_string(), self.inner.ca_density.into_py(py));
            dict
        })
    }
}

impl From<StructureSummary> for PyStructureSummary {
    fn from(summary: StructureSummary) -> Self {
        PyStructureSummary { inner: summary }
    }
}
