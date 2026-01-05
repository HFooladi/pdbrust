//! Python bindings for structure descriptors

use pdbrust::descriptors::StructureDescriptors;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Structural descriptors computed from a PDB structure
///
/// Contains geometric, compositional, and size metrics.
#[pyclass(name = "StructureDescriptors")]
#[derive(Clone)]
pub struct PyStructureDescriptors {
    inner: StructureDescriptors,
}

#[pymethods]
impl PyStructureDescriptors {
    /// Number of residues (based on CA count)
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues
    }

    /// Number of atoms
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

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

    /// Fraction of hydrophobic residues (0.0 to 1.0)
    #[getter]
    fn hydrophobic_ratio(&self) -> f64 {
        self.inner.hydrophobic_ratio
    }

    /// Fraction of glycine residues (0.0 to 1.0)
    #[getter]
    fn glycine_ratio(&self) -> f64 {
        self.inner.glycine_ratio
    }

    /// Missing residue ratio based on sequence gaps
    #[getter]
    fn missing_residue_ratio(&self) -> f64 {
        self.inner.missing_residue_ratio
    }

    /// Secondary structure ratio (heuristic estimate)
    #[getter]
    fn secondary_structure_ratio(&self) -> f64 {
        self.inner.secondary_structure_ratio
    }

    /// Compactness index (Rg / n^(1/3))
    #[getter]
    fn compactness_index(&self) -> f64 {
        self.inner.compactness_index
    }

    /// CA density (atoms per cubic Angstrom)
    #[getter]
    fn ca_density(&self) -> f64 {
        self.inner.ca_density
    }

    /// Amino acid composition as dictionary
    #[getter]
    fn aa_composition(&self) -> HashMap<String, f64> {
        self.inner.aa_composition.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "StructureDescriptors(residues={}, Rg={:.2}, max_dist={:.2}, hydrophobic={:.2})",
            self.inner.num_residues,
            self.inner.radius_of_gyration,
            self.inner.max_ca_distance,
            self.inner.hydrophobic_ratio
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, PyObject> {
        Python::with_gil(|py| {
            let mut dict = HashMap::new();
            dict.insert("num_residues".to_string(), self.inner.num_residues.into_py(py));
            dict.insert("num_atoms".to_string(), self.inner.num_atoms.into_py(py));
            dict.insert("radius_of_gyration".to_string(), self.inner.radius_of_gyration.into_py(py));
            dict.insert("max_ca_distance".to_string(), self.inner.max_ca_distance.into_py(py));
            dict.insert("hydrophobic_ratio".to_string(), self.inner.hydrophobic_ratio.into_py(py));
            dict.insert("glycine_ratio".to_string(), self.inner.glycine_ratio.into_py(py));
            dict.insert("missing_residue_ratio".to_string(), self.inner.missing_residue_ratio.into_py(py));
            dict.insert("secondary_structure_ratio".to_string(), self.inner.secondary_structure_ratio.into_py(py));
            dict.insert("compactness_index".to_string(), self.inner.compactness_index.into_py(py));
            dict.insert("ca_density".to_string(), self.inner.ca_density.into_py(py));
            dict
        })
    }
}

impl From<StructureDescriptors> for PyStructureDescriptors {
    fn from(desc: StructureDescriptors) -> Self {
        PyStructureDescriptors { inner: desc }
    }
}
