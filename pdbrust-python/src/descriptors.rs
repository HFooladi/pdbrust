//! Python bindings for structure descriptors

use pdbrust::descriptors::{ResidueBFactor, StructureDescriptors};
use pyo3::prelude::*;
use std::collections::HashMap;

/// Per-residue B-factor statistics
///
/// Contains B-factor mean, min, max and atom count for a single residue.
#[pyclass(name = "ResidueBFactor")]
#[derive(Clone)]
pub struct PyResidueBFactor {
    pub(crate) inner: ResidueBFactor,
}

#[pymethods]
impl PyResidueBFactor {
    /// Chain identifier (e.g., "A")
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    /// Insertion code (if any)
    #[getter]
    fn ins_code(&self) -> Option<char> {
        self.inner.ins_code
    }

    /// Residue name (e.g., "ALA", "GLY")
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// Mean B-factor across all atoms in the residue (Å²)
    #[getter]
    fn b_factor_mean(&self) -> f64 {
        self.inner.b_factor_mean
    }

    /// Minimum B-factor among atoms in the residue (Å²)
    #[getter]
    fn b_factor_min(&self) -> f64 {
        self.inner.b_factor_min
    }

    /// Maximum B-factor among atoms in the residue (Å²)
    #[getter]
    fn b_factor_max(&self) -> f64 {
        self.inner.b_factor_max
    }

    /// Number of atoms in the residue
    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count
    }

    fn __repr__(&self) -> String {
        format!(
            "ResidueBFactor(chain={}, seq={}, name={}, mean={:.2})",
            self.inner.chain_id,
            self.inner.residue_seq,
            self.inner.residue_name,
            self.inner.b_factor_mean
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, PyObject> {
        Python::with_gil(|py| {
            let mut dict = HashMap::new();
            dict.insert(
                "chain_id".to_string(),
                self.inner.chain_id.clone().into_py(py),
            );
            dict.insert(
                "residue_seq".to_string(),
                self.inner.residue_seq.into_py(py),
            );
            dict.insert("ins_code".to_string(), self.inner.ins_code.into_py(py));
            dict.insert(
                "residue_name".to_string(),
                self.inner.residue_name.clone().into_py(py),
            );
            dict.insert(
                "b_factor_mean".to_string(),
                self.inner.b_factor_mean.into_py(py),
            );
            dict.insert(
                "b_factor_min".to_string(),
                self.inner.b_factor_min.into_py(py),
            );
            dict.insert(
                "b_factor_max".to_string(),
                self.inner.b_factor_max.into_py(py),
            );
            dict.insert("atom_count".to_string(), self.inner.atom_count.into_py(py));
            dict
        })
    }
}

impl From<ResidueBFactor> for PyResidueBFactor {
    fn from(res: ResidueBFactor) -> Self {
        PyResidueBFactor { inner: res }
    }
}

impl From<&ResidueBFactor> for PyResidueBFactor {
    fn from(res: &ResidueBFactor) -> Self {
        PyResidueBFactor { inner: res.clone() }
    }
}

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

    /// Mean B-factor (temperature factor) across all atoms (Å²)
    #[getter]
    fn b_factor_mean(&self) -> f64 {
        self.inner.b_factor_mean
    }

    /// Mean B-factor for Cα atoms only (Å²)
    #[getter]
    fn b_factor_mean_ca(&self) -> f64 {
        self.inner.b_factor_mean_ca
    }

    /// Minimum B-factor in the structure (Å²)
    #[getter]
    fn b_factor_min(&self) -> f64 {
        self.inner.b_factor_min
    }

    /// Maximum B-factor in the structure (Å²)
    #[getter]
    fn b_factor_max(&self) -> f64 {
        self.inner.b_factor_max
    }

    /// Standard deviation of B-factors (Å²)
    #[getter]
    fn b_factor_std(&self) -> f64 {
        self.inner.b_factor_std
    }

    fn __repr__(&self) -> String {
        format!(
            "StructureDescriptors(residues={}, Rg={:.2}, max_dist={:.2}, hydrophobic={:.2}, b_factor_mean={:.2})",
            self.inner.num_residues,
            self.inner.radius_of_gyration,
            self.inner.max_ca_distance,
            self.inner.hydrophobic_ratio,
            self.inner.b_factor_mean
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, PyObject> {
        Python::with_gil(|py| {
            let mut dict = HashMap::new();
            dict.insert(
                "num_residues".to_string(),
                self.inner.num_residues.into_py(py),
            );
            dict.insert("num_atoms".to_string(), self.inner.num_atoms.into_py(py));
            dict.insert(
                "radius_of_gyration".to_string(),
                self.inner.radius_of_gyration.into_py(py),
            );
            dict.insert(
                "max_ca_distance".to_string(),
                self.inner.max_ca_distance.into_py(py),
            );
            dict.insert(
                "hydrophobic_ratio".to_string(),
                self.inner.hydrophobic_ratio.into_py(py),
            );
            dict.insert(
                "glycine_ratio".to_string(),
                self.inner.glycine_ratio.into_py(py),
            );
            dict.insert(
                "missing_residue_ratio".to_string(),
                self.inner.missing_residue_ratio.into_py(py),
            );
            dict.insert(
                "secondary_structure_ratio".to_string(),
                self.inner.secondary_structure_ratio.into_py(py),
            );
            dict.insert(
                "compactness_index".to_string(),
                self.inner.compactness_index.into_py(py),
            );
            dict.insert("ca_density".to_string(), self.inner.ca_density.into_py(py));
            dict.insert(
                "b_factor_mean".to_string(),
                self.inner.b_factor_mean.into_py(py),
            );
            dict.insert(
                "b_factor_mean_ca".to_string(),
                self.inner.b_factor_mean_ca.into_py(py),
            );
            dict.insert(
                "b_factor_min".to_string(),
                self.inner.b_factor_min.into_py(py),
            );
            dict.insert(
                "b_factor_max".to_string(),
                self.inner.b_factor_max.into_py(py),
            );
            dict.insert(
                "b_factor_std".to_string(),
                self.inner.b_factor_std.into_py(py),
            );
            dict
        })
    }
}

impl From<StructureDescriptors> for PyStructureDescriptors {
    fn from(desc: StructureDescriptors) -> Self {
        PyStructureDescriptors { inner: desc }
    }
}
