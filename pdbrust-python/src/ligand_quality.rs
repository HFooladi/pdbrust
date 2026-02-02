//! Python bindings for ligand pose quality assessment.
//!
//! Provides PoseBusters-style geometry checks for protein-ligand complexes.

use pdbrust::ligand_quality::{AtomClash, LigandPoseReport};
use pyo3::prelude::*;

/// Represents a steric clash between two atoms.
///
/// A clash occurs when two non-bonded atoms are closer than the expected
/// minimum distance based on their van der Waals radii.
#[pyclass(name = "AtomClash")]
#[derive(Clone)]
pub struct PyAtomClash {
    inner: AtomClash,
}

#[pymethods]
impl PyAtomClash {
    /// Serial number of the protein/receptor atom.
    #[getter]
    fn protein_atom_serial(&self) -> i32 {
        self.inner.protein_atom_serial
    }

    /// Serial number of the ligand atom.
    #[getter]
    fn ligand_atom_serial(&self) -> i32 {
        self.inner.ligand_atom_serial
    }

    /// Chain ID of the protein atom.
    #[getter]
    fn protein_chain_id(&self) -> String {
        self.inner.protein_chain_id.clone()
    }

    /// Residue name of the protein atom.
    #[getter]
    fn protein_residue_name(&self) -> String {
        self.inner.protein_residue_name.clone()
    }

    /// Residue sequence number of the protein atom.
    #[getter]
    fn protein_residue_seq(&self) -> i32 {
        self.inner.protein_residue_seq
    }

    /// Atom name of the protein atom.
    #[getter]
    fn protein_atom_name(&self) -> String {
        self.inner.protein_atom_name.clone()
    }

    /// Element of the protein atom.
    #[getter]
    fn protein_element(&self) -> String {
        self.inner.protein_element.clone()
    }

    /// Atom name of the ligand atom.
    #[getter]
    fn ligand_atom_name(&self) -> String {
        self.inner.ligand_atom_name.clone()
    }

    /// Element of the ligand atom.
    #[getter]
    fn ligand_element(&self) -> String {
        self.inner.ligand_element.clone()
    }

    /// Actual distance between atoms in Angstroms.
    #[getter]
    fn distance(&self) -> f64 {
        self.inner.distance
    }

    /// Expected minimum distance (0.75 × sum of vdW radii) in Angstroms.
    #[getter]
    fn expected_min_distance(&self) -> f64 {
        self.inner.expected_min_distance
    }

    /// Severity ratio: expected_min_distance / distance.
    /// Higher values indicate more severe clashes.
    #[getter]
    fn severity(&self) -> f64 {
        self.inner.severity
    }

    fn __repr__(&self) -> String {
        format!(
            "AtomClash(protein={} {} {}, ligand={}, dist={:.2}Å, expected={:.2}Å, severity={:.2}x)",
            self.inner.protein_chain_id,
            self.inner.protein_residue_name,
            self.inner.protein_residue_seq,
            self.inner.ligand_atom_name.trim(),
            self.inner.distance,
            self.inner.expected_min_distance,
            self.inner.severity
        )
    }
}

impl From<AtomClash> for PyAtomClash {
    fn from(inner: AtomClash) -> Self {
        Self { inner }
    }
}

/// Comprehensive report on ligand pose quality.
///
/// Contains results from PoseBusters-style geometry checks including
/// steric clash detection and volume overlap calculation.
#[pyclass(name = "LigandPoseReport")]
#[derive(Clone)]
pub struct PyLigandPoseReport {
    inner: LigandPoseReport,
}

#[pymethods]
impl PyLigandPoseReport {
    /// Residue name of the ligand (e.g., "LIG", "ATP").
    #[getter]
    fn ligand_name(&self) -> String {
        self.inner.ligand_name.clone()
    }

    /// Chain ID where the ligand is located.
    #[getter]
    fn ligand_chain_id(&self) -> String {
        self.inner.ligand_chain_id.clone()
    }

    /// Residue sequence number of the ligand.
    #[getter]
    fn ligand_residue_seq(&self) -> i32 {
        self.inner.ligand_residue_seq
    }

    /// Number of atoms in the ligand.
    #[getter]
    fn ligand_atom_count(&self) -> usize {
        self.inner.ligand_atom_count
    }

    /// Minimum distance between any ligand atom and any protein atom (Å).
    #[getter]
    fn min_protein_ligand_distance(&self) -> f64 {
        self.inner.min_protein_ligand_distance
    }

    /// List of detected steric clashes with protein atoms.
    #[getter]
    fn clashes(&self) -> Vec<PyAtomClash> {
        self.inner
            .clashes
            .iter()
            .cloned()
            .map(PyAtomClash::from)
            .collect()
    }

    /// Whether any protein clashes were detected.
    #[getter]
    fn has_protein_clash(&self) -> bool {
        self.inner.has_protein_clash
    }

    /// Total number of protein-ligand clashes.
    #[getter]
    fn num_clashes(&self) -> usize {
        self.inner.num_clashes
    }

    /// Severity of the worst clash (expected_min / actual distance).
    #[getter]
    fn worst_clash_severity(&self) -> f64 {
        self.inner.worst_clash_severity
    }

    /// Percentage of ligand volume overlapping with protein (0-100%).
    #[getter]
    fn protein_volume_overlap_pct(&self) -> f64 {
        self.inner.protein_volume_overlap_pct
    }

    /// List of detected clashes with other HETATM atoms (cofactors).
    #[getter]
    fn cofactor_clashes(&self) -> Vec<PyAtomClash> {
        self.inner
            .cofactor_clashes
            .iter()
            .cloned()
            .map(PyAtomClash::from)
            .collect()
    }

    /// Whether any cofactor clashes were detected.
    #[getter]
    fn has_cofactor_clash(&self) -> bool {
        self.inner.has_cofactor_clash
    }

    /// Number of clashes with cofactors.
    #[getter]
    fn num_cofactor_clashes(&self) -> usize {
        self.inner.num_cofactor_clashes
    }

    /// Whether the pose passes the minimum distance check.
    #[getter]
    fn passes_distance_check(&self) -> bool {
        self.inner.passes_distance_check
    }

    /// Whether the pose passes the volume overlap check (<7.5%).
    #[getter]
    fn passes_overlap_check(&self) -> bool {
        self.inner.passes_overlap_check
    }

    /// Whether the pose passes all geometry checks.
    #[getter]
    fn is_geometry_valid(&self) -> bool {
        self.inner.is_geometry_valid
    }

    /// Get a summary string of the report.
    fn summary(&self) -> String {
        self.inner.summary()
    }

    fn __repr__(&self) -> String {
        let status = if self.inner.is_geometry_valid {
            "PASS"
        } else {
            "FAIL"
        };
        format!(
            "LigandPoseReport(ligand='{}', chain='{}', seq={}, status={}, clashes={}, overlap={:.1}%)",
            self.inner.ligand_name,
            self.inner.ligand_chain_id,
            self.inner.ligand_residue_seq,
            status,
            self.inner.num_clashes,
            self.inner.protein_volume_overlap_pct
        )
    }
}

impl From<LigandPoseReport> for PyLigandPoseReport {
    fn from(inner: LigandPoseReport) -> Self {
        Self { inner }
    }
}

/// Get the van der Waals radius for an element.
///
/// Args:
///     element: Element symbol (case-insensitive)
///
/// Returns:
///     The van der Waals radius in Angstroms.
#[pyfunction]
pub fn vdw_radius(element: &str) -> f64 {
    pdbrust::ligand_quality::vdw_radius(element)
}

/// Get the covalent radius for an element.
///
/// Args:
///     element: Element symbol (case-insensitive)
///
/// Returns:
///     The covalent radius in Angstroms.
#[pyfunction]
pub fn covalent_radius(element: &str) -> f64 {
    pdbrust::ligand_quality::covalent_radius(element)
}
