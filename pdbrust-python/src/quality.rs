//! Python bindings for structure quality assessment

use pdbrust::quality::QualityReport;
use pyo3::conversion::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::Py;
use std::collections::HashMap;

/// Quality assessment report for a PDB structure
///
/// Contains various quality metrics and flags indicating
/// structure characteristics.
#[pyclass(name = "QualityReport")]
#[derive(Clone)]
pub struct PyQualityReport {
    inner: QualityReport,
}

#[pymethods]
impl PyQualityReport {
    /// True if structure appears to be CA-only (coarse-grained)
    #[getter]
    fn has_ca_only(&self) -> bool {
        self.inner.has_ca_only
    }

    /// True if structure has multiple models (NMR ensemble)
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

    /// Number of atoms
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    /// Number of residues
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues
    }

    /// True if structure contains HETATM records (ligands, waters)
    #[getter]
    fn has_hetatm(&self) -> bool {
        self.inner.has_hetatm
    }

    /// True if structure contains hydrogen atoms
    #[getter]
    fn has_hydrogens(&self) -> bool {
        self.inner.has_hydrogens
    }

    /// True if structure has disulfide bonds defined
    #[getter]
    fn has_ssbonds(&self) -> bool {
        self.inner.has_ssbonds
    }

    /// True if structure has CONECT records
    #[getter]
    fn has_conect(&self) -> bool {
        self.inner.has_conect
    }

    /// Check if structure is suitable for typical analysis
    ///
    /// Returns True if:
    /// - Single model (not NMR ensemble)
    /// - No alternate locations
    /// - Not CA-only
    ///
    /// Returns:
    ///     True if structure is analysis-ready
    fn is_analysis_ready(&self) -> bool {
        self.inner.is_analysis_ready()
    }

    /// Check if structure is "clean"
    ///
    /// Returns True if:
    /// - No alternate locations
    /// - Not CA-only
    ///
    /// Returns:
    ///     True if structure is clean
    fn is_clean(&self) -> bool {
        self.inner.is_clean()
    }

    fn __repr__(&self) -> String {
        format!(
            "QualityReport(atoms={}, chains={}, models={}, analysis_ready={})",
            self.inner.num_atoms,
            self.inner.num_chains,
            self.inner.num_models,
            self.inner.is_analysis_ready()
        )
    }

    /// Convert to dictionary
    fn to_dict(&self) -> HashMap<String, Py<PyAny>> {
        Python::attach(|py| {
            let mut dict = HashMap::new();
            dict.insert(
                "has_ca_only".to_string(),
                self.inner.has_ca_only.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_multiple_models".to_string(),
                self.inner.has_multiple_models.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_altlocs".to_string(),
                self.inner.has_altlocs.into_py_any(py).unwrap(),
            );
            dict.insert(
                "num_chains".to_string(),
                self.inner.num_chains.into_py_any(py).unwrap(),
            );
            dict.insert(
                "num_models".to_string(),
                self.inner.num_models.into_py_any(py).unwrap(),
            );
            dict.insert(
                "num_atoms".to_string(),
                self.inner.num_atoms.into_py_any(py).unwrap(),
            );
            dict.insert(
                "num_residues".to_string(),
                self.inner.num_residues.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_hetatm".to_string(),
                self.inner.has_hetatm.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_hydrogens".to_string(),
                self.inner.has_hydrogens.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_ssbonds".to_string(),
                self.inner.has_ssbonds.into_py_any(py).unwrap(),
            );
            dict.insert(
                "has_conect".to_string(),
                self.inner.has_conect.into_py_any(py).unwrap(),
            );
            dict.insert(
                "is_analysis_ready".to_string(),
                self.inner.is_analysis_ready().into_py_any(py).unwrap(),
            );
            dict.insert(
                "is_clean".to_string(),
                self.inner.is_clean().into_py_any(py).unwrap(),
            );
            dict
        })
    }
}

impl From<QualityReport> for PyQualityReport {
    fn from(report: QualityReport) -> Self {
        PyQualityReport { inner: report }
    }
}
