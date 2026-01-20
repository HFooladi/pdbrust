//! Python bindings for DSSP secondary structure assignment

use pdbrust::dssp::{ResidueSSAssignment, SecondaryStructure, SecondaryStructureAssignment};
use pyo3::prelude::*;

/// Secondary structure classification (DSSP codes).
///
/// The DSSP algorithm classifies residues into 9 states:
/// - H: α-helix (i → i+4 H-bond pattern)
/// - G: 3₁₀-helix (i → i+3 H-bond pattern)
/// - I: π-helix (i → i+5 H-bond pattern)
/// - P: κ-helix / PPII helix (polyproline II, dihedral-based)
/// - E: Extended strand in β-sheet
/// - B: Isolated β-bridge
/// - T: Hydrogen-bonded turn
/// - S: Bend (high backbone curvature)
/// - C: Coil (none of the above)
#[pyclass(name = "SecondaryStructure")]
#[derive(Clone)]
pub struct PySecondaryStructure {
    pub(crate) inner: SecondaryStructure,
}

#[pymethods]
impl PySecondaryStructure {
    /// α-helix: i → i+4 hydrogen bond pattern
    #[classattr]
    fn ALPHA_HELIX() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::AlphaHelix,
        }
    }

    /// 3₁₀-helix: i → i+3 hydrogen bond pattern
    #[classattr]
    fn HELIX_310() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::Helix310,
        }
    }

    /// π-helix: i → i+5 hydrogen bond pattern
    #[classattr]
    fn PI_HELIX() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::PiHelix,
        }
    }

    /// κ-helix / PPII helix: polyproline II helix (dihedral-based)
    #[classattr]
    fn KAPPA_HELIX() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::KappaHelix,
        }
    }

    /// Extended strand: part of β-sheet
    #[classattr]
    fn EXTENDED_STRAND() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::ExtendedStrand,
        }
    }

    /// Isolated β-bridge: single H-bond pair
    #[classattr]
    fn BETA_BRIDGE() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::BetaBridge,
        }
    }

    /// Turn: H-bonded turn, not part of helix
    #[classattr]
    fn TURN() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::Turn,
        }
    }

    /// Bend: high backbone curvature
    #[classattr]
    fn BEND() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::Bend,
        }
    }

    /// Coil: none of the above
    #[classattr]
    fn COIL() -> Self {
        PySecondaryStructure {
            inner: SecondaryStructure::Coil,
        }
    }

    /// Get the single-character DSSP code.
    ///
    /// Returns:
    ///     Single character code (H, G, I, P, E, B, T, S, or C)
    fn code(&self) -> String {
        self.inner.code().to_string()
    }

    /// Check if this is a helical secondary structure.
    ///
    /// Returns:
    ///     True for H, G, I, or P
    fn is_helix(&self) -> bool {
        self.inner.is_helix()
    }

    /// Check if this is a β-structure.
    ///
    /// Returns:
    ///     True for E or B
    fn is_sheet(&self) -> bool {
        self.inner.is_sheet()
    }

    /// Check if this is coil/loop.
    ///
    /// Returns:
    ///     True for C, T, or S
    fn is_coil(&self) -> bool {
        self.inner.is_coil()
    }

    fn __repr__(&self) -> String {
        format!("SecondaryStructure('{}')", self.inner.code())
    }

    fn __str__(&self) -> String {
        self.inner.code().to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut hasher = DefaultHasher::new();
        self.inner.code().hash(&mut hasher);
        hasher.finish()
    }
}

impl From<SecondaryStructure> for PySecondaryStructure {
    fn from(ss: SecondaryStructure) -> Self {
        PySecondaryStructure { inner: ss }
    }
}

/// Secondary structure assignment for a single residue.
#[pyclass(name = "ResidueSSAssignment")]
#[derive(Clone)]
pub struct PyResidueSSAssignment {
    pub(crate) inner: ResidueSSAssignment,
}

#[pymethods]
impl PyResidueSSAssignment {
    /// Chain identifier
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    /// Residue name (3-letter code)
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// Insertion code (if any)
    #[getter]
    fn ins_code(&self) -> Option<char> {
        self.inner.ins_code
    }

    /// Assigned secondary structure
    #[getter]
    fn ss(&self) -> PySecondaryStructure {
        PySecondaryStructure::from(self.inner.ss)
    }

    /// Get the DSSP code for this residue
    fn code(&self) -> String {
        self.inner.ss.code().to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "ResidueSSAssignment(chain='{}', resid={}, name='{}', ss='{}')",
            self.inner.chain_id,
            self.inner.residue_seq,
            self.inner.residue_name,
            self.inner.ss.code()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}{}:{} ({})",
            self.inner.chain_id,
            self.inner.residue_seq,
            self.inner.residue_name,
            self.inner.ss.code()
        )
    }
}

impl From<ResidueSSAssignment> for PyResidueSSAssignment {
    fn from(assignment: ResidueSSAssignment) -> Self {
        PyResidueSSAssignment { inner: assignment }
    }
}

/// Complete secondary structure assignment for a protein structure.
///
/// Contains per-residue assignments and summary statistics.
#[pyclass(name = "SecondaryStructureAssignment")]
#[derive(Clone)]
pub struct PySecondaryStructureAssignment {
    pub(crate) inner: SecondaryStructureAssignment,
}

#[pymethods]
impl PySecondaryStructureAssignment {
    /// Number of residues with assignments
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Check if there are no assignments
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// Get per-residue assignments
    #[getter]
    fn residue_assignments(&self) -> Vec<PyResidueSSAssignment> {
        self.inner
            .residue_assignments
            .iter()
            .cloned()
            .map(PyResidueSSAssignment::from)
            .collect()
    }

    /// Number of residues assigned to helix (H, G, I, P)
    #[getter]
    fn helix_count(&self) -> usize {
        self.inner.helix_count
    }

    /// Number of residues assigned to sheet (E, B)
    #[getter]
    fn sheet_count(&self) -> usize {
        self.inner.sheet_count
    }

    /// Number of residues assigned to coil (C, T, S)
    #[getter]
    fn coil_count(&self) -> usize {
        self.inner.coil_count
    }

    /// Fraction of residues in helical conformation
    #[getter]
    fn helix_fraction(&self) -> f64 {
        self.inner.helix_fraction
    }

    /// Fraction of residues in sheet conformation
    #[getter]
    fn sheet_fraction(&self) -> f64 {
        self.inner.sheet_fraction
    }

    /// Fraction of residues in coil conformation
    #[getter]
    fn coil_fraction(&self) -> f64 {
        self.inner.coil_fraction
    }

    /// Warnings generated during assignment
    #[getter]
    fn warnings(&self) -> Vec<String> {
        self.inner.warnings.clone()
    }

    /// Get the secondary structure as a string of single-character codes.
    ///
    /// Returns:
    ///     String with one character per residue (e.g., "HHHHEEEECCCC")
    fn to_string(&self) -> String {
        self.inner.as_codes()
    }

    /// Get the secondary structure composition.
    ///
    /// Returns:
    ///     Tuple of (helix_fraction, sheet_fraction, coil_fraction)
    fn composition(&self) -> (f64, f64, f64) {
        self.inner.composition()
    }

    fn __repr__(&self) -> String {
        format!(
            "SecondaryStructureAssignment(residues={}, helix={:.1}%, sheet={:.1}%, coil={:.1}%)",
            self.inner.len(),
            self.inner.helix_fraction * 100.0,
            self.inner.sheet_fraction * 100.0,
            self.inner.coil_fraction * 100.0
        )
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Iterate over residue assignments
    fn __iter__(slf: PyRef<'_, Self>) -> PyResidueSSAssignmentIterator {
        PyResidueSSAssignmentIterator {
            inner: slf.inner.residue_assignments.clone(),
            index: 0,
        }
    }

    /// Get assignment by index
    fn __getitem__(&self, index: isize) -> PyResult<PyResidueSSAssignment> {
        let len = self.inner.len() as isize;
        let normalized_index = if index < 0 { len + index } else { index };

        if normalized_index < 0 || normalized_index >= len {
            return Err(pyo3::exceptions::PyIndexError::new_err(
                "index out of range",
            ));
        }

        Ok(PyResidueSSAssignment::from(
            self.inner.residue_assignments[normalized_index as usize].clone(),
        ))
    }
}

impl From<SecondaryStructureAssignment> for PySecondaryStructureAssignment {
    fn from(assignment: SecondaryStructureAssignment) -> Self {
        PySecondaryStructureAssignment { inner: assignment }
    }
}

/// Iterator for residue assignments
#[pyclass]
pub struct PyResidueSSAssignmentIterator {
    inner: Vec<ResidueSSAssignment>,
    index: usize,
}

#[pymethods]
impl PyResidueSSAssignmentIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyResidueSSAssignment> {
        if slf.index < slf.inner.len() {
            let item = slf.inner[slf.index].clone();
            slf.index += 1;
            Some(PyResidueSSAssignment::from(item))
        } else {
            None
        }
    }
}
