//! Python bindings for geometry/RMSD/LDDT functionality

use pdbrust::geometry::{
    AlignmentResult, AtomSelection, LddtOptions, LddtResult, PerResidueLddt, PerResidueRmsd,
};
use pyo3::prelude::*;

/// Atom selection criteria for RMSD and alignment calculations.
///
/// Use the class methods to create selection instances:
///
///     selection = AtomSelection.ca_only()    # Default: CA atoms only
///     selection = AtomSelection.backbone()   # N, CA, C, O atoms
///     selection = AtomSelection.all_atoms()  # All atoms
///     selection = AtomSelection.custom(["CA", "CB"])  # Custom list
#[pyclass(name = "AtomSelection")]
#[derive(Clone)]
pub struct PyAtomSelection {
    pub(crate) inner: AtomSelection,
}

#[pymethods]
impl PyAtomSelection {
    /// Create a selection for CA (alpha-carbon) atoms only.
    ///
    /// This is the default and most common choice for protein structure
    /// comparison. CA atoms provide a reliable backbone representation.
    ///
    /// Returns:
    ///     AtomSelection: A selection for CA atoms only
    #[staticmethod]
    fn ca_only() -> Self {
        PyAtomSelection {
            inner: AtomSelection::CaOnly,
        }
    }

    /// Create a selection for backbone atoms (N, CA, C, O).
    ///
    /// This provides more detailed backbone comparison than CA-only,
    /// capturing backbone geometry more completely.
    ///
    /// Returns:
    ///     AtomSelection: A selection for backbone atoms
    #[staticmethod]
    fn backbone() -> Self {
        PyAtomSelection {
            inner: AtomSelection::Backbone,
        }
    }

    /// Create a selection for all atoms.
    ///
    /// Requires exact atom correspondence between structures.
    /// Use with caution as minor differences in atom naming can cause issues.
    ///
    /// Returns:
    ///     AtomSelection: A selection for all atoms
    #[staticmethod]
    fn all_atoms() -> Self {
        PyAtomSelection {
            inner: AtomSelection::AllAtoms,
        }
    }

    /// Create a custom selection from a list of atom names.
    ///
    /// Args:
    ///     atom_names: List of atom names to include (e.g., ["CA", "CB"])
    ///
    /// Returns:
    ///     AtomSelection: A custom selection
    #[staticmethod]
    fn custom(atom_names: Vec<String>) -> Self {
        PyAtomSelection {
            inner: AtomSelection::Custom(atom_names),
        }
    }

    fn __repr__(&self) -> String {
        match &self.inner {
            AtomSelection::CaOnly => "AtomSelection.ca_only()".to_string(),
            AtomSelection::Backbone => "AtomSelection.backbone()".to_string(),
            AtomSelection::AllAtoms => "AtomSelection.all_atoms()".to_string(),
            AtomSelection::Custom(names) => format!("AtomSelection.custom({:?})", names),
        }
    }
}

/// Result of structure alignment containing RMSD and transformation parameters.
///
/// Attributes:
///     rmsd (float): Root mean square deviation after alignment (Angstroms)
///     rotation (list): 3x3 rotation matrix as nested lists
///     translation (list): Translation vector [tx, ty, tz]
///     num_atoms (int): Number of atoms used in the alignment
#[pyclass(name = "AlignmentResult")]
#[derive(Clone)]
pub struct PyAlignmentResult {
    inner: AlignmentResult,
}

#[pymethods]
impl PyAlignmentResult {
    /// RMSD after optimal alignment (Angstroms)
    #[getter]
    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    /// 3x3 rotation matrix as nested lists
    #[getter]
    fn rotation(&self) -> [[f64; 3]; 3] {
        self.inner.rotation
    }

    /// Translation vector [tx, ty, tz]
    #[getter]
    fn translation(&self) -> [f64; 3] {
        self.inner.translation
    }

    /// Number of atoms used in the alignment
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentResult(rmsd={:.4}, num_atoms={})",
            self.inner.rmsd, self.inner.num_atoms
        )
    }

    fn __str__(&self) -> String {
        format!(
            "RMSD: {:.4} Angstroms ({} atoms)",
            self.inner.rmsd, self.inner.num_atoms
        )
    }
}

impl From<AlignmentResult> for PyAlignmentResult {
    fn from(result: AlignmentResult) -> Self {
        PyAlignmentResult { inner: result }
    }
}

/// Per-residue RMSD information for flexibility analysis.
///
/// Attributes:
///     chain_id (str): Chain identifier
///     residue_seq (int): Residue sequence number
///     residue_name (str): Residue name (e.g., "ALA", "GLY")
///     rmsd (float): RMSD for this residue (Angstroms)
///     num_atoms (int): Number of atoms used for this residue
#[pyclass(name = "PerResidueRmsd")]
#[derive(Clone)]
pub struct PyPerResidueRmsd {
    inner: PerResidueRmsd,
}

#[pymethods]
impl PyPerResidueRmsd {
    /// Chain identifier
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.residue_id.0
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_id.1
    }

    /// Residue name (e.g., "ALA", "GLY")
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// RMSD for this residue (Angstroms)
    #[getter]
    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    /// Number of atoms used for this residue
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    fn __repr__(&self) -> String {
        format!(
            "PerResidueRmsd({}{} {}: {:.3} A)",
            self.inner.residue_id.0,
            self.inner.residue_id.1,
            self.inner.residue_name,
            self.inner.rmsd
        )
    }
}

impl From<PerResidueRmsd> for PyPerResidueRmsd {
    fn from(result: PerResidueRmsd) -> Self {
        PyPerResidueRmsd { inner: result }
    }
}

// ============================================================================
// LDDT Types
// ============================================================================

/// Options for LDDT (Local Distance Difference Test) calculation.
///
/// LDDT is a superposition-free metric for comparing protein structures,
/// widely used in AlphaFold (pLDDT) and CASP evaluations.
///
/// Example:
///     >>> options = LddtOptions()  # Default options
///     >>> options = LddtOptions(inclusion_radius=10.0)  # Custom radius
///     >>> options = LddtOptions(thresholds=[0.5, 1.0, 2.0])  # Custom thresholds
///
/// Attributes:
///     inclusion_radius (float): Maximum distance to consider (default: 15.0 Angstroms)
///     thresholds (list[float]): Distance difference thresholds (default: [0.5, 1.0, 2.0, 4.0])
#[pyclass(name = "LddtOptions")]
#[derive(Clone)]
pub struct PyLddtOptions {
    pub(crate) inner: LddtOptions,
}

#[pymethods]
impl PyLddtOptions {
    /// Create LDDT options.
    ///
    /// Args:
    ///     inclusion_radius: Maximum distance to consider for local comparisons.
    ///                       Pairs of atoms with reference distance > inclusion_radius
    ///                       are ignored. Default: 15.0 Angstroms.
    ///     thresholds: Distance difference thresholds. For each threshold, count
    ///                 how many distance pairs are preserved (|d_ref - d_model| < threshold).
    ///                 Default: [0.5, 1.0, 2.0, 4.0] Angstroms.
    ///
    /// Returns:
    ///     LddtOptions: Options for LDDT calculation
    #[new]
    #[pyo3(signature = (inclusion_radius=15.0, thresholds=None))]
    fn new(inclusion_radius: f64, thresholds: Option<Vec<f64>>) -> Self {
        let thresholds = thresholds.unwrap_or_else(|| vec![0.5, 1.0, 2.0, 4.0]);
        PyLddtOptions {
            inner: LddtOptions {
                inclusion_radius,
                thresholds,
            },
        }
    }

    /// Maximum distance to consider for local comparisons (Angstroms)
    #[getter]
    fn inclusion_radius(&self) -> f64 {
        self.inner.inclusion_radius
    }

    /// Distance difference thresholds (Angstroms)
    #[getter]
    fn thresholds(&self) -> Vec<f64> {
        self.inner.thresholds.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "LddtOptions(inclusion_radius={:.1}, thresholds={:?})",
            self.inner.inclusion_radius, self.inner.thresholds
        )
    }
}

impl From<LddtOptions> for PyLddtOptions {
    fn from(options: LddtOptions) -> Self {
        PyLddtOptions { inner: options }
    }
}

/// Result of LDDT (Local Distance Difference Test) calculation.
///
/// LDDT scores range from 0.0 (poor) to 1.0 (perfect).
///
/// Attributes:
///     score (float): Global LDDT score (0.0 to 1.0)
///     num_pairs (int): Number of distance pairs evaluated
///     per_threshold_scores (list[float]): Score for each threshold
///     num_residues (int): Number of residues evaluated
#[pyclass(name = "LddtResult")]
#[derive(Clone)]
pub struct PyLddtResult {
    inner: LddtResult,
}

#[pymethods]
impl PyLddtResult {
    /// Global LDDT score (0.0 to 1.0)
    #[getter]
    fn score(&self) -> f64 {
        self.inner.score
    }

    /// Number of distance pairs evaluated
    #[getter]
    fn num_pairs(&self) -> usize {
        self.inner.num_pairs
    }

    /// Score for each threshold
    #[getter]
    fn per_threshold_scores(&self) -> Vec<f64> {
        self.inner.per_threshold_scores.clone()
    }

    /// Number of residues evaluated
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.num_residues
    }

    fn __repr__(&self) -> String {
        format!(
            "LddtResult(score={:.4}, num_pairs={}, num_residues={})",
            self.inner.score, self.inner.num_pairs, self.inner.num_residues
        )
    }

    fn __str__(&self) -> String {
        format!(
            "LDDT: {:.4} ({} pairs, {} residues)",
            self.inner.score, self.inner.num_pairs, self.inner.num_residues
        )
    }
}

impl From<LddtResult> for PyLddtResult {
    fn from(result: LddtResult) -> Self {
        PyLddtResult { inner: result }
    }
}

/// Per-residue LDDT information for quality analysis.
///
/// Identifies poorly modeled regions in predicted structures.
///
/// Attributes:
///     chain_id (str): Chain identifier
///     residue_seq (int): Residue sequence number
///     residue_name (str): Residue name (e.g., "ALA", "GLY")
///     score (float): LDDT score for this residue (0.0 to 1.0)
///     num_pairs (int): Number of distance pairs involving this residue
#[pyclass(name = "PerResidueLddt")]
#[derive(Clone)]
pub struct PyPerResidueLddt {
    inner: PerResidueLddt,
}

#[pymethods]
impl PyPerResidueLddt {
    /// Chain identifier
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.residue_id.0
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_id.1
    }

    /// Residue name (e.g., "ALA", "GLY")
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// LDDT score for this residue (0.0 to 1.0)
    #[getter]
    fn score(&self) -> f64 {
        self.inner.score
    }

    /// Number of distance pairs involving this residue
    #[getter]
    fn num_pairs(&self) -> usize {
        self.inner.num_pairs
    }

    fn __repr__(&self) -> String {
        format!(
            "PerResidueLddt({}{} {}: {:.3})",
            self.inner.residue_id.0,
            self.inner.residue_id.1,
            self.inner.residue_name,
            self.inner.score
        )
    }
}

impl From<PerResidueLddt> for PyPerResidueLddt {
    fn from(result: PerResidueLddt) -> Self {
        PyPerResidueLddt { inner: result }
    }
}
