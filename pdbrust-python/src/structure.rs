//! Python bindings for PdbStructure

use pdbrust::PdbStructure;
use pyo3::prelude::*;
use std::collections::HashMap;

use crate::atom::PyAtom;
use crate::error::convert_error;
use crate::records::{PyConect, PyModel, PyRemark, PySSBond, PySeqRes};

/// Represents a parsed PDB/mmCIF structure.
///
/// This is the main container for all structure data including atoms,
/// chains, residues, and metadata.
#[pyclass(name = "PdbStructure")]
pub struct PyPdbStructure {
    pub(crate) inner: PdbStructure,
}

#[pymethods]
impl PyPdbStructure {
    /// Create a new empty structure
    #[new]
    fn new() -> Self {
        PyPdbStructure {
            inner: PdbStructure::new(),
        }
    }

    // ==================== Basic Properties ====================

    /// Number of atoms in the structure
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.get_num_atoms()
    }

    /// Number of chains in the structure
    #[getter]
    fn num_chains(&self) -> usize {
        self.inner.get_num_chains()
    }

    /// Number of residues in the structure
    #[getter]
    fn num_residues(&self) -> usize {
        self.inner.get_num_residues()
    }

    /// Header string (if present)
    #[getter]
    fn header(&self) -> Option<&str> {
        self.inner.header.as_deref()
    }

    /// Title string (if present)
    #[getter]
    fn title(&self) -> Option<&str> {
        self.inner.title.as_deref()
    }

    // ==================== Atom Access ====================

    /// Get all atoms as a list
    #[getter]
    fn atoms(&self) -> Vec<PyAtom> {
        self.inner.atoms.iter().map(PyAtom::from).collect()
    }

    // ==================== Record Access ====================

    /// Get all disulfide bonds
    #[getter]
    fn ssbonds(&self) -> Vec<PySSBond> {
        self.inner.ssbonds.iter().map(PySSBond::from).collect()
    }

    /// Get all SEQRES records
    #[getter]
    fn seqres(&self) -> Vec<PySeqRes> {
        self.inner.seqres.iter().map(PySeqRes::from).collect()
    }

    /// Get all CONECT records
    #[getter]
    fn connects(&self) -> Vec<PyConect> {
        self.inner.connects.iter().map(PyConect::from).collect()
    }

    /// Get all REMARK records
    #[getter]
    fn remarks(&self) -> Vec<PyRemark> {
        self.inner.remarks.iter().map(PyRemark::from).collect()
    }

    /// Get all models (for multi-model structures like NMR ensembles)
    #[getter]
    fn models(&self) -> Vec<PyModel> {
        self.inner.models.iter().map(PyModel::from).collect()
    }

    // ==================== Query Methods ====================

    /// Get sorted list of unique chain identifiers
    ///
    /// Returns:
    ///     List of chain IDs (e.g., ["A", "B", "C"])
    fn get_chain_ids(&self) -> Vec<String> {
        self.inner.get_chain_ids()
    }

    /// Get all residues as list of (seq_num, residue_name) tuples
    ///
    /// Returns:
    ///     List of (residue_seq, residue_name) tuples
    fn get_residues(&self) -> Vec<(i32, String)> {
        self.inner.get_residues()
    }

    /// Get residues for a specific chain
    ///
    /// Args:
    ///     chain_id: Chain identifier (e.g., "A")
    ///
    /// Returns:
    ///     List of (residue_seq, residue_name) tuples for that chain
    fn get_residues_for_chain(&self, chain_id: &str) -> Vec<(i32, String)> {
        self.inner.get_residues_for_chain(chain_id)
    }

    /// Get sequence for a chain from SEQRES records
    ///
    /// Args:
    ///     chain_id: Chain identifier
    ///
    /// Returns:
    ///     List of residue names from SEQRES
    fn get_sequence(&self, chain_id: &str) -> Vec<String> {
        self.inner.get_sequence(chain_id)
    }

    /// Get remarks by number
    ///
    /// Args:
    ///     number: Remark number (e.g., 2 for resolution)
    ///
    /// Returns:
    ///     List of Remark objects with that number
    fn get_remarks_by_number(&self, number: i32) -> Vec<PyRemark> {
        self.inner
            .get_remarks_by_number(number)
            .into_iter()
            .map(PyRemark::from)
            .collect()
    }

    /// Get atoms connected to a specific atom (from CONECT records)
    ///
    /// Args:
    ///     atom_serial: Serial number of the atom
    ///
    /// Returns:
    ///     List of connected Atom objects
    fn get_connected_atoms(&self, atom_serial: i32) -> Vec<PyAtom> {
        self.inner
            .get_connected_atoms(atom_serial)
            .into_iter()
            .map(PyAtom::from)
            .collect()
    }

    // ==================== Modification Methods ====================

    /// Translate all atoms by (dx, dy, dz)
    ///
    /// Args:
    ///     dx: Translation in X
    ///     dy: Translation in Y
    ///     dz: Translation in Z
    fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        self.inner.translate(dx, dy, dz);
    }

    /// Write structure to a PDB file
    ///
    /// Args:
    ///     path: Output file path
    ///
    /// Raises:
    ///     IOError: If file cannot be written
    fn to_file(&self, path: &str) -> PyResult<()> {
        self.inner.to_file(path).map_err(convert_error)
    }

    // ==================== Filter Methods (feature-gated) ====================

    /// Get CA coordinates, optionally filtered by chain
    ///
    /// Args:
    ///     chain_id: Optional chain identifier to filter by
    ///
    /// Returns:
    ///     List of (x, y, z) coordinate tuples
    #[cfg(feature = "filter")]
    #[pyo3(signature = (chain_id=None))]
    fn get_ca_coords(&self, chain_id: Option<&str>) -> Vec<(f64, f64, f64)> {
        self.inner.get_ca_coords(chain_id)
    }

    /// Remove ligands, waters, and non-standard residues
    ///
    /// Returns:
    ///     New PdbStructure with only standard amino acids/nucleotides
    #[cfg(feature = "filter")]
    fn remove_ligands(&self) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.remove_ligands(),
        }
    }

    /// Keep only atoms from the specified chain
    ///
    /// Args:
    ///     chain_id: Chain identifier to keep
    ///
    /// Returns:
    ///     New PdbStructure with only that chain
    #[cfg(feature = "filter")]
    fn keep_only_chain(&self, chain_id: &str) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.keep_only_chain(chain_id),
        }
    }

    /// Keep only CA (alpha carbon) atoms
    ///
    /// Returns:
    ///     New PdbStructure with only CA atoms
    #[cfg(feature = "filter")]
    fn keep_only_ca(&self) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.keep_only_ca(),
        }
    }

    /// Keep only backbone atoms (N, CA, C, O)
    ///
    /// Returns:
    ///     New PdbStructure with only backbone atoms
    #[cfg(feature = "filter")]
    fn keep_only_backbone(&self) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.keep_only_backbone(),
        }
    }

    /// Remove hydrogen atoms
    ///
    /// Returns:
    ///     New PdbStructure without hydrogen atoms
    #[cfg(feature = "filter")]
    fn remove_hydrogens(&self) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.remove_hydrogens(),
        }
    }

    /// Normalize chain IDs to A, B, C, ...
    #[cfg(feature = "filter")]
    fn normalize_chain_ids(&mut self) {
        self.inner.normalize_chain_ids();
    }

    /// Reindex residues starting from 1
    #[cfg(feature = "filter")]
    fn reindex_residues(&mut self) {
        self.inner.reindex_residues();
    }

    /// Renumber atoms sequentially starting from 1
    #[cfg(feature = "filter")]
    fn renumber_atoms(&mut self) {
        self.inner.renumber_atoms();
    }

    /// Center structure at the origin
    #[cfg(feature = "filter")]
    fn center_structure(&mut self) {
        self.inner.center_structure();
    }

    /// Get the centroid (center of mass) of all atoms
    ///
    /// Returns:
    ///     (x, y, z) coordinates of the centroid
    #[cfg(feature = "filter")]
    fn get_centroid(&self) -> (f64, f64, f64) {
        self.inner.get_centroid()
    }

    /// Get the centroid of CA atoms only
    ///
    /// Returns:
    ///     (x, y, z) coordinates of the CA centroid
    #[cfg(feature = "filter")]
    fn get_ca_centroid(&self) -> (f64, f64, f64) {
        self.inner.get_ca_centroid()
    }

    // ==================== Descriptor Methods (feature-gated) ====================

    /// Calculate radius of gyration
    ///
    /// Returns:
    ///     Radius of gyration in Angstroms
    #[cfg(feature = "descriptors")]
    fn radius_of_gyration(&self) -> f64 {
        self.inner.radius_of_gyration()
    }

    /// Calculate maximum CA-CA distance
    ///
    /// Returns:
    ///     Maximum distance between any two CA atoms in Angstroms
    #[cfg(feature = "descriptors")]
    fn max_ca_distance(&self) -> f64 {
        self.inner.max_ca_distance()
    }

    /// Get amino acid composition
    ///
    /// Returns:
    ///     Dictionary mapping residue names to their fractions (sum to 1.0)
    #[cfg(feature = "descriptors")]
    fn aa_composition(&self) -> HashMap<String, f64> {
        self.inner.aa_composition()
    }

    /// Get glycine ratio
    ///
    /// Returns:
    ///     Fraction of residues that are glycine (0.0 to 1.0)
    #[cfg(feature = "descriptors")]
    fn glycine_ratio(&self) -> f64 {
        self.inner.glycine_ratio()
    }

    /// Get hydrophobic residue ratio
    ///
    /// Returns:
    ///     Fraction of residues that are hydrophobic (0.0 to 1.0)
    #[cfg(feature = "descriptors")]
    fn hydrophobic_ratio(&self) -> f64 {
        self.inner.hydrophobic_ratio()
    }

    /// Get polar residue ratio
    ///
    /// Returns:
    ///     Fraction of residues that are polar (0.0 to 1.0)
    #[cfg(feature = "descriptors")]
    fn polar_ratio(&self) -> f64 {
        self.inner.polar_ratio()
    }

    /// Get charged residue ratio
    ///
    /// Returns:
    ///     Fraction of residues that are charged (0.0 to 1.0)
    #[cfg(feature = "descriptors")]
    fn charged_ratio(&self) -> f64 {
        self.inner.charged_ratio()
    }

    /// Count CA residues (proxy for total residues)
    ///
    /// Returns:
    ///     Number of CA atoms
    #[cfg(feature = "descriptors")]
    fn count_ca_residues(&self) -> usize {
        self.inner.count_ca_residues()
    }

    /// Get CA density (atoms per cubic Angstrom)
    ///
    /// Returns:
    ///     CA atom density
    #[cfg(feature = "descriptors")]
    fn ca_density(&self) -> f64 {
        self.inner.ca_density()
    }

    /// Get compactness index (Rg / n^(1/3))
    ///
    /// Returns:
    ///     Compactness index value
    #[cfg(feature = "descriptors")]
    fn compactness_index(&self) -> f64 {
        self.inner.compactness_index()
    }

    /// Get all structure descriptors at once
    ///
    /// Returns:
    ///     StructureDescriptors object with all metrics
    #[cfg(feature = "descriptors")]
    fn structure_descriptors(&self) -> crate::descriptors::PyStructureDescriptors {
        crate::descriptors::PyStructureDescriptors::from(self.inner.structure_descriptors())
    }

    // ==================== Quality Methods (feature-gated) ====================

    /// Check if structure has only CA atoms (coarse-grained)
    ///
    /// Returns:
    ///     True if structure appears to be CA-only
    #[cfg(feature = "quality")]
    fn has_ca_only(&self) -> bool {
        self.inner.has_ca_only()
    }

    /// Check if structure has multiple models (NMR ensemble)
    ///
    /// Returns:
    ///     True if multiple models present
    #[cfg(feature = "quality")]
    fn has_multiple_models(&self) -> bool {
        self.inner.has_multiple_models()
    }

    /// Check if structure has alternate locations
    ///
    /// Returns:
    ///     True if any atoms have alternate locations
    #[cfg(feature = "quality")]
    fn has_altlocs(&self) -> bool {
        self.inner.has_altlocs()
    }

    /// Check if structure has HETATM records
    ///
    /// Returns:
    ///     True if ligands/waters present
    #[cfg(feature = "quality")]
    fn has_hetatm(&self) -> bool {
        self.inner.has_hetatm()
    }

    /// Check if structure has hydrogen atoms
    ///
    /// Returns:
    ///     True if hydrogens present
    #[cfg(feature = "quality")]
    fn has_hydrogens(&self) -> bool {
        self.inner.has_hydrogens()
    }

    /// Get resolution from REMARK 2 (if available)
    ///
    /// Returns:
    ///     Resolution in Angstroms, or None if not available
    #[cfg(feature = "quality")]
    fn get_resolution(&self) -> Option<f64> {
        self.inner.get_resolution()
    }

    /// Get comprehensive quality report
    ///
    /// Returns:
    ///     QualityReport object with all quality metrics
    #[cfg(feature = "quality")]
    fn quality_report(&self) -> crate::quality::PyQualityReport {
        crate::quality::PyQualityReport::from(self.inner.quality_report())
    }

    // ==================== Summary Methods (feature-gated) ====================

    /// Get unified structure summary (quality + descriptors)
    ///
    /// Returns:
    ///     StructureSummary object with all metrics
    #[cfg(feature = "summary")]
    fn summary(&self) -> crate::summary::PyStructureSummary {
        crate::summary::PyStructureSummary::from(self.inner.summary())
    }

    // ==================== Numpy Methods (feature-gated) ====================

    /// Get all atom coordinates as a numpy array (N x 3)
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, 3) with atom coordinates
    #[cfg(feature = "numpy")]
    fn get_coords_array<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        crate::numpy_support::get_coords_array(py, &self.inner)
    }

    /// Get CA atom coordinates as a numpy array (N x 3)
    ///
    /// Args:
    ///     chain_id: Optional chain identifier to filter by
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, 3) with CA coordinates
    #[cfg(all(feature = "numpy", feature = "filter"))]
    #[pyo3(signature = (chain_id=None))]
    fn get_ca_coords_array<'py>(
        &self,
        py: Python<'py>,
        chain_id: Option<&str>,
    ) -> Bound<'py, numpy::PyArray2<f64>> {
        crate::numpy_support::get_ca_coords_array(py, &self.inner, chain_id)
    }

    /// Get backbone atom coordinates as a numpy array (N x 3)
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, 3) with backbone atom coordinates
    #[cfg(all(feature = "numpy", feature = "filter"))]
    fn get_backbone_coords_array<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        crate::numpy_support::get_backbone_coords_array(py, &self.inner)
    }

    // ==================== Magic Methods ====================

    fn __repr__(&self) -> String {
        format!(
            "PdbStructure(atoms={}, chains={}, residues={})",
            self.inner.get_num_atoms(),
            self.inner.get_num_chains(),
            self.inner.get_num_residues()
        )
    }

    fn __str__(&self) -> String {
        let header = self.inner.header.as_deref().unwrap_or("Unknown");
        format!(
            "PdbStructure: {} ({} atoms, {} chains)",
            header,
            self.inner.get_num_atoms(),
            self.inner.get_num_chains()
        )
    }

    fn __len__(&self) -> usize {
        self.inner.get_num_atoms()
    }
}

impl From<PdbStructure> for PyPdbStructure {
    fn from(structure: PdbStructure) -> Self {
        PyPdbStructure { inner: structure }
    }
}
