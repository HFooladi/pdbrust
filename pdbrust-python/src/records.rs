//! Python bindings for PDB record types (SSBond, SeqRes, Conect, Remark, Model)

use pdbrust::{Conect, Model, Remark, SSBond, SeqRes};
use pyo3::prelude::*;

use crate::atom::PyAtom;

/// Represents a disulfide bond between two cysteine residues
#[pyclass(name = "SSBond")]
#[derive(Clone)]
pub struct PySSBond {
    inner: SSBond,
}

#[pymethods]
impl PySSBond {
    /// Serial number of the SSBond record
    #[getter]
    fn serial(&self) -> i32 {
        self.inner.serial
    }

    /// Chain ID of the first cysteine
    #[getter]
    fn chain1_id(&self) -> &str {
        &self.inner.chain1_id
    }

    /// Residue sequence number of the first cysteine
    #[getter]
    fn residue1_seq(&self) -> i32 {
        self.inner.residue1_seq
    }

    /// Insertion code of the first cysteine
    #[getter]
    fn icode1(&self) -> Option<char> {
        self.inner.icode1
    }

    /// Chain ID of the second cysteine
    #[getter]
    fn chain2_id(&self) -> &str {
        &self.inner.chain2_id
    }

    /// Residue sequence number of the second cysteine
    #[getter]
    fn residue2_seq(&self) -> i32 {
        self.inner.residue2_seq
    }

    /// Insertion code of the second cysteine
    #[getter]
    fn icode2(&self) -> Option<char> {
        self.inner.icode2
    }

    /// Bond length in Angstroms
    #[getter]
    fn length(&self) -> f64 {
        self.inner.length
    }

    fn __repr__(&self) -> String {
        format!(
            "SSBond({}:{} - {}:{}, length={:.2})",
            self.inner.chain1_id,
            self.inner.residue1_seq,
            self.inner.chain2_id,
            self.inner.residue2_seq,
            self.inner.length
        )
    }
}

impl From<SSBond> for PySSBond {
    fn from(ssbond: SSBond) -> Self {
        PySSBond { inner: ssbond }
    }
}

impl From<&SSBond> for PySSBond {
    fn from(ssbond: &SSBond) -> Self {
        PySSBond {
            inner: ssbond.clone(),
        }
    }
}

/// Represents a SEQRES record containing sequence information
#[pyclass(name = "SeqRes")]
#[derive(Clone)]
pub struct PySeqRes {
    inner: SeqRes,
}

#[pymethods]
impl PySeqRes {
    /// Chain identifier
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Number of residues in the chain
    #[getter]
    fn num_residues(&self) -> i32 {
        self.inner.num_residues
    }

    /// List of residue names
    #[getter]
    fn residues(&self) -> Vec<String> {
        self.inner.residues.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "SeqRes(chain='{}', num_residues={}, residues={})",
            self.inner.chain_id,
            self.inner.num_residues,
            self.inner.residues.len()
        )
    }
}

impl From<SeqRes> for PySeqRes {
    fn from(seqres: SeqRes) -> Self {
        PySeqRes { inner: seqres }
    }
}

impl From<&SeqRes> for PySeqRes {
    fn from(seqres: &SeqRes) -> Self {
        PySeqRes {
            inner: seqres.clone(),
        }
    }
}

/// Represents a CONECT record containing atom connectivity information
#[pyclass(name = "Conect")]
#[derive(Clone)]
pub struct PyConect {
    inner: Conect,
}

#[pymethods]
impl PyConect {
    /// Serial number of atom 1
    #[getter]
    fn atom1(&self) -> i32 {
        self.inner.atom1
    }

    /// Serial number of atom 2
    #[getter]
    fn atom2(&self) -> i32 {
        self.inner.atom2
    }

    /// Serial number of atom 3 (optional)
    #[getter]
    fn atom3(&self) -> Option<i32> {
        self.inner.atom3
    }

    /// Serial number of atom 4 (optional)
    #[getter]
    fn atom4(&self) -> Option<i32> {
        self.inner.atom4
    }

    /// Get all connected atom serials as a list
    fn bonded_atoms(&self) -> Vec<i32> {
        let mut atoms = vec![self.inner.atom1, self.inner.atom2];
        if let Some(a3) = self.inner.atom3 {
            atoms.push(a3);
        }
        if let Some(a4) = self.inner.atom4 {
            atoms.push(a4);
        }
        atoms
    }

    fn __repr__(&self) -> String {
        format!(
            "Conect(atom1={}, atom2={}, atom3={:?}, atom4={:?})",
            self.inner.atom1, self.inner.atom2, self.inner.atom3, self.inner.atom4
        )
    }
}

impl From<Conect> for PyConect {
    fn from(conect: Conect) -> Self {
        PyConect { inner: conect }
    }
}

impl From<&Conect> for PyConect {
    fn from(conect: &Conect) -> Self {
        PyConect {
            inner: conect.clone(),
        }
    }
}

/// Represents a REMARK record
#[pyclass(name = "Remark")]
#[derive(Clone)]
pub struct PyRemark {
    inner: Remark,
}

#[pymethods]
impl PyRemark {
    /// Remark number
    #[getter]
    fn number(&self) -> i32 {
        self.inner.number
    }

    /// Remark content
    #[getter]
    fn content(&self) -> &str {
        &self.inner.content
    }

    fn __repr__(&self) -> String {
        let content_preview = if self.inner.content.len() > 40 {
            format!("{}...", &self.inner.content[..40])
        } else {
            self.inner.content.clone()
        };
        format!("Remark(number={}, content='{}')", self.inner.number, content_preview)
    }
}

impl From<Remark> for PyRemark {
    fn from(remark: Remark) -> Self {
        PyRemark { inner: remark }
    }
}

impl From<&Remark> for PyRemark {
    fn from(remark: &Remark) -> Self {
        PyRemark {
            inner: remark.clone(),
        }
    }
}

/// Represents a MODEL record for multi-model structures (e.g., NMR ensembles)
#[pyclass(name = "Model")]
#[derive(Clone)]
pub struct PyModel {
    inner: Model,
}

#[pymethods]
impl PyModel {
    /// Model serial number
    #[getter]
    fn serial(&self) -> i32 {
        self.inner.serial
    }

    /// Atoms in this model
    #[getter]
    fn atoms(&self) -> Vec<PyAtom> {
        self.inner.atoms.iter().map(PyAtom::from).collect()
    }

    /// Number of atoms in this model
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.atoms.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Model(serial={}, atoms={})",
            self.inner.serial,
            self.inner.atoms.len()
        )
    }

    fn __len__(&self) -> usize {
        self.inner.atoms.len()
    }
}

impl From<Model> for PyModel {
    fn from(model: Model) -> Self {
        PyModel { inner: model }
    }
}

impl From<&Model> for PyModel {
    fn from(model: &Model) -> Self {
        PyModel {
            inner: model.clone(),
        }
    }
}
