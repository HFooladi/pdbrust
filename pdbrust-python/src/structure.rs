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

    /// Select atoms using a PyMOL/VMD-style selection language.
    ///
    /// This method provides a powerful and flexible way to filter atoms
    /// using a query language inspired by molecular visualization tools.
    ///
    /// Basic Selectors:
    ///     - chain A: Select atoms from chain A
    ///     - name CA: Select atoms named CA (alpha carbons)
    ///     - resname ALA: Select alanine residues
    ///     - resid 50: Select residue number 50
    ///     - resid 1:100: Select residues 1 through 100
    ///     - element C: Select carbon atoms
    ///
    /// Keywords:
    ///     - backbone: N, CA, C, O atoms
    ///     - protein: Standard amino acids
    ///     - nucleic: Standard nucleotides
    ///     - water: Water molecules (HOH, WAT)
    ///     - hetero: HETATM records
    ///     - hydrogen: Hydrogen atoms
    ///     - all or *: All atoms
    ///
    /// Numeric Comparisons:
    ///     - bfactor < 30.0: B-factor less than 30
    ///     - occupancy >= 0.5: Occupancy at least 0.5
    ///
    /// Boolean Operations:
    ///     - and: Both conditions must match
    ///     - or: Either condition matches
    ///     - not: Negation
    ///     - (): Grouping with parentheses
    ///
    /// Args:
    ///     selection: A selection string (e.g., "chain A and name CA")
    ///
    /// Returns:
    ///     New PdbStructure containing only the selected atoms
    ///
    /// Raises:
    ///     ValueError: If the selection string is invalid
    ///
    /// Example:
    ///     >>> # Simple selections
    ///     >>> chain_a = structure.select("chain A")
    ///     >>> ca_atoms = structure.select("name CA")
    ///     >>>
    ///     >>> # Combined selections
    ///     >>> chain_a_ca = structure.select("chain A and name CA")
    ///     >>>
    ///     >>> # Complex selections
    ///     >>> selected = structure.select("(chain A or chain B) and backbone")
    ///     >>>
    ///     >>> # Numeric comparisons
    ///     >>> low_bfactor = structure.select("protein and bfactor < 30.0")
    #[cfg(feature = "filter")]
    fn select(&self, selection: &str) -> PyResult<PyPdbStructure> {
        self.inner
            .select(selection)
            .map(|s| PyPdbStructure { inner: s })
            .map_err(crate::error::convert_selection_error)
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

    // ==================== B-factor Methods (feature-gated) ====================

    /// Calculate mean B-factor across all atoms
    ///
    /// Returns:
    ///     Mean B-factor (temperature factor) in Å²
    #[cfg(feature = "descriptors")]
    fn b_factor_mean(&self) -> f64 {
        self.inner.b_factor_mean()
    }

    /// Calculate mean B-factor for CA atoms only
    ///
    /// Returns:
    ///     Mean B-factor for CA atoms in Å²
    #[cfg(feature = "descriptors")]
    fn b_factor_mean_ca(&self) -> f64 {
        self.inner.b_factor_mean_ca()
    }

    /// Get minimum B-factor value
    ///
    /// Returns:
    ///     Minimum B-factor in Å²
    #[cfg(feature = "descriptors")]
    fn b_factor_min(&self) -> f64 {
        self.inner.b_factor_min()
    }

    /// Get maximum B-factor value
    ///
    /// Returns:
    ///     Maximum B-factor in Å²
    #[cfg(feature = "descriptors")]
    fn b_factor_max(&self) -> f64 {
        self.inner.b_factor_max()
    }

    /// Calculate standard deviation of B-factors
    ///
    /// Returns:
    ///     B-factor standard deviation in Å²
    #[cfg(feature = "descriptors")]
    fn b_factor_std(&self) -> f64 {
        self.inner.b_factor_std()
    }

    /// Get per-residue B-factor statistics
    ///
    /// For each residue, computes mean, min, and max B-factors
    /// across all atoms in that residue.
    ///
    /// Returns:
    ///     List of ResidueBFactor objects, sorted by chain and residue number
    ///
    /// Example:
    ///     >>> profile = structure.b_factor_profile()
    ///     >>> for res in profile:
    ///     ...     print(f"{res.chain_id}{res.residue_seq}: {res.b_factor_mean:.2f}")
    #[cfg(feature = "descriptors")]
    fn b_factor_profile(&self) -> Vec<crate::descriptors::PyResidueBFactor> {
        self.inner
            .b_factor_profile()
            .into_iter()
            .map(crate::descriptors::PyResidueBFactor::from)
            .collect()
    }

    /// Identify flexible residues (high B-factors)
    ///
    /// Returns residues where mean B-factor exceeds the threshold.
    /// High B-factors indicate mobile or disordered regions.
    ///
    /// Args:
    ///     threshold: B-factor cutoff in Å². Residues with mean B > threshold
    ///                are returned. Common values: 40-60 Å².
    ///
    /// Returns:
    ///     List of ResidueBFactor for flexible residues
    ///
    /// Example:
    ///     >>> flexible = structure.flexible_residues(50.0)
    ///     >>> print(f"Found {len(flexible)} flexible residues")
    #[cfg(feature = "descriptors")]
    fn flexible_residues(&self, threshold: f64) -> Vec<crate::descriptors::PyResidueBFactor> {
        self.inner
            .flexible_residues(threshold)
            .into_iter()
            .map(crate::descriptors::PyResidueBFactor::from)
            .collect()
    }

    /// Identify rigid residues (low B-factors)
    ///
    /// Returns residues where mean B-factor is below the threshold.
    /// Low B-factors indicate well-ordered, rigid regions.
    ///
    /// Args:
    ///     threshold: B-factor cutoff in Å². Residues with mean B < threshold
    ///                are returned. Common values: 15-25 Å².
    ///
    /// Returns:
    ///     List of ResidueBFactor for rigid residues
    ///
    /// Example:
    ///     >>> rigid = structure.rigid_residues(20.0)
    ///     >>> print(f"Found {len(rigid)} rigid residues")
    #[cfg(feature = "descriptors")]
    fn rigid_residues(&self, threshold: f64) -> Vec<crate::descriptors::PyResidueBFactor> {
        self.inner
            .rigid_residues(threshold)
            .into_iter()
            .map(crate::descriptors::PyResidueBFactor::from)
            .collect()
    }

    /// Create a new structure with Z-score normalized B-factors
    ///
    /// Transforms B-factors to have mean 0 and standard deviation 1,
    /// enabling comparison between different structures.
    ///
    /// Returns:
    ///     New PdbStructure with normalized B-factors
    ///
    /// Example:
    ///     >>> normalized = structure.normalize_b_factors()
    ///     >>> print(f"Normalized mean: {normalized.b_factor_mean():.4f}")
    ///     >>> print(f"Normalized std:  {normalized.b_factor_std():.4f}")
    #[cfg(feature = "descriptors")]
    fn normalize_b_factors(&self) -> PyPdbStructure {
        PyPdbStructure {
            inner: self.inner.normalize_b_factors(),
        }
    }

    /// Get the percentile rank of an atom's B-factor
    ///
    /// Args:
    ///     atom_serial: Serial number of the atom to query
    ///
    /// Returns:
    ///     Percentile (0-100) of the atom's B-factor, or None if not found
    ///
    /// Example:
    ///     >>> percentile = structure.b_factor_percentile(42)
    ///     >>> if percentile is not None:
    ///     ...     print(f"Atom 42 is at the {percentile:.1f}th percentile")
    #[cfg(feature = "descriptors")]
    fn b_factor_percentile(&self, atom_serial: i32) -> Option<f64> {
        self.inner.b_factor_percentile(atom_serial)
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

    // ==================== Distance Matrix / Contact Map Methods ====================

    /// Get the CA distance matrix as a numpy array (N x N)
    ///
    /// Computes pairwise Euclidean distances between all CA atoms.
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, N) with distances in Angstroms
    #[cfg(all(feature = "numpy", feature = "descriptors"))]
    fn distance_matrix_ca<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        let matrix = self.inner.distance_matrix_ca();
        crate::numpy_support::vec2d_to_array2(py, &matrix)
    }

    /// Get the CA contact map as a numpy boolean array (N x N)
    ///
    /// Args:
    ///     threshold: Distance threshold in Angstroms (default: 8.0)
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, N) with True where distance <= threshold
    #[cfg(all(feature = "numpy", feature = "descriptors"))]
    #[pyo3(signature = (threshold=8.0))]
    fn contact_map_ca<'py>(
        &self,
        py: Python<'py>,
        threshold: f64,
    ) -> Bound<'py, numpy::PyArray2<bool>> {
        let matrix = self.inner.contact_map_ca(threshold);
        crate::numpy_support::vec2d_bool_to_array2(py, &matrix)
    }

    /// Get the all-atom distance matrix as a numpy array (N x N)
    ///
    /// Computes pairwise Euclidean distances between all atoms.
    /// Warning: This can be memory-intensive for large structures.
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, N) with distances in Angstroms
    #[cfg(all(feature = "numpy", feature = "descriptors"))]
    fn distance_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        let matrix = self.inner.distance_matrix();
        crate::numpy_support::vec2d_to_array2(py, &matrix)
    }

    /// Get the all-atom contact map as a numpy boolean array (N x N)
    ///
    /// Args:
    ///     threshold: Distance threshold in Angstroms (default: 4.5)
    ///
    /// Returns:
    ///     numpy.ndarray of shape (N, N) with True where distance <= threshold
    #[cfg(all(feature = "numpy", feature = "descriptors"))]
    #[pyo3(signature = (threshold=4.5))]
    fn contact_map<'py>(
        &self,
        py: Python<'py>,
        threshold: f64,
    ) -> Bound<'py, numpy::PyArray2<bool>> {
        let matrix = self.inner.contact_map(threshold);
        crate::numpy_support::vec2d_bool_to_array2(py, &matrix)
    }

    // ==================== Geometry Methods (feature-gated) ====================

    /// Calculate RMSD to another structure.
    ///
    /// This calculates RMSD without alignment. For aligned RMSD, use align_to()
    /// which returns both the aligned structure and RMSD.
    ///
    /// Args:
    ///     other: Target structure to compare against
    ///     selection: Atom selection (default: CA atoms only).
    ///                Use AtomSelection.ca_only(), backbone(), all_atoms(), or custom()
    ///
    /// Returns:
    ///     RMSD in Angstroms
    ///
    /// Raises:
    ///     ValueError: If structures have different numbers of selected atoms
    ///
    /// Example:
    ///     >>> rmsd = structure1.rmsd_to(structure2)
    ///     >>> rmsd = structure1.rmsd_to(structure2, AtomSelection.backbone())
    #[cfg(feature = "geometry")]
    #[pyo3(signature = (other, selection=None))]
    fn rmsd_to(
        &self,
        other: &PyPdbStructure,
        selection: Option<&crate::geometry::PyAtomSelection>,
    ) -> PyResult<f64> {
        let sel = selection
            .map(|s| s.inner.clone())
            .unwrap_or(pdbrust::geometry::AtomSelection::CaOnly);
        self.inner
            .rmsd_to_with_selection(&other.inner, sel)
            .map_err(convert_error)
    }

    /// Align this structure onto a target structure.
    ///
    /// Uses the Kabsch algorithm to find optimal rotation and translation
    /// that minimizes RMSD between the structures.
    ///
    /// Args:
    ///     target: Reference structure to align onto
    ///     selection: Atom selection for alignment (default: CA atoms only)
    ///
    /// Returns:
    ///     Tuple of (aligned_structure, AlignmentResult)
    ///
    /// Raises:
    ///     ValueError: If structures have different numbers of selected atoms,
    ///                 or fewer than 3 atoms are selected
    ///
    /// Example:
    ///     >>> aligned, result = mobile.align_to(target)
    ///     >>> print(f"RMSD: {result.rmsd:.4f} Angstroms")
    ///     >>> aligned.to_file("aligned.pdb")
    #[cfg(feature = "geometry")]
    #[pyo3(signature = (target, selection=None))]
    fn align_to(
        &self,
        target: &PyPdbStructure,
        selection: Option<&crate::geometry::PyAtomSelection>,
    ) -> PyResult<(PyPdbStructure, crate::geometry::PyAlignmentResult)> {
        let sel = selection
            .map(|s| s.inner.clone())
            .unwrap_or(pdbrust::geometry::AtomSelection::CaOnly);
        let (aligned, result) = self
            .inner
            .align_to_with_selection(&target.inner, sel)
            .map_err(convert_error)?;
        Ok((
            PyPdbStructure { inner: aligned },
            crate::geometry::PyAlignmentResult::from(result),
        ))
    }

    /// Get per-residue RMSD after alignment to target.
    ///
    /// This method first aligns the structures optimally, then computes
    /// RMSD for each residue individually. Useful for identifying flexible
    /// regions in protein structures.
    ///
    /// Args:
    ///     target: Reference structure
    ///     selection: Atom selection (default: CA atoms only)
    ///
    /// Returns:
    ///     List of PerResidueRmsd objects
    ///
    /// Example:
    ///     >>> per_res = mobile.per_residue_rmsd_to(target)
    ///     >>> for r in per_res:
    ///     ...     if r.rmsd > 2.0:
    ///     ...         print(f"{r.chain_id}{r.residue_seq}: {r.rmsd:.2f} A")
    #[cfg(feature = "geometry")]
    #[pyo3(signature = (target, selection=None))]
    fn per_residue_rmsd_to(
        &self,
        target: &PyPdbStructure,
        selection: Option<&crate::geometry::PyAtomSelection>,
    ) -> PyResult<Vec<crate::geometry::PyPerResidueRmsd>> {
        let sel = selection
            .map(|s| s.inner.clone())
            .unwrap_or(pdbrust::geometry::AtomSelection::CaOnly);
        self.inner
            .per_residue_rmsd_to_with_selection(&target.inner, sel)
            .map(|v| {
                v.into_iter()
                    .map(crate::geometry::PyPerResidueRmsd::from)
                    .collect()
            })
            .map_err(convert_error)
    }

    /// Calculate LDDT (Local Distance Difference Test) to a reference structure.
    ///
    /// LDDT is a superposition-free metric that measures the fraction of
    /// inter-atomic distances that are preserved within specified thresholds.
    /// It ranges from 0.0 (poor) to 1.0 (perfect).
    ///
    /// This is the same metric used by AlphaFold (pLDDT) and CASP evaluations.
    ///
    /// Args:
    ///     reference: The reference structure (ground truth)
    ///     selection: Atom selection (default: CA atoms only).
    ///                Use AtomSelection.ca_only(), backbone(), all_atoms(), or custom()
    ///     options: LDDT options (default: inclusion_radius=15.0, thresholds=[0.5, 1.0, 2.0, 4.0])
    ///
    /// Returns:
    ///     LddtResult containing:
    ///     - score: Global LDDT score (0.0 to 1.0)
    ///     - num_pairs: Number of distance pairs evaluated
    ///     - per_threshold_scores: Score for each threshold
    ///     - num_residues: Number of residues evaluated
    ///
    /// Raises:
    ///     ValueError: If structures have different numbers of selected atoms
    ///
    /// Example:
    ///     >>> result = model.lddt_to(reference)
    ///     >>> print(f"LDDT: {result.score:.4f}")
    ///     >>> result = model.lddt_to(reference, AtomSelection.backbone())
    ///     >>> result = model.lddt_to(reference, options=LddtOptions(inclusion_radius=10.0))
    #[cfg(feature = "geometry")]
    #[pyo3(signature = (reference, selection=None, options=None))]
    fn lddt_to(
        &self,
        reference: &PyPdbStructure,
        selection: Option<&crate::geometry::PyAtomSelection>,
        options: Option<&crate::geometry::PyLddtOptions>,
    ) -> PyResult<crate::geometry::PyLddtResult> {
        let sel = selection
            .map(|s| s.inner.clone())
            .unwrap_or(pdbrust::geometry::AtomSelection::CaOnly);
        let opts = options
            .map(|o| o.inner.clone())
            .unwrap_or_else(pdbrust::geometry::LddtOptions::default);
        self.inner
            .lddt_to_with_options(&reference.inner, sel, opts)
            .map(crate::geometry::PyLddtResult::from)
            .map_err(convert_error)
    }

    /// Get per-residue LDDT scores.
    ///
    /// Returns LDDT scores for each residue individually, useful for
    /// identifying poorly modeled regions in predicted structures.
    ///
    /// Args:
    ///     reference: The reference structure (ground truth)
    ///     selection: Atom selection (default: CA atoms only)
    ///     options: LDDT options (default: inclusion_radius=15.0, thresholds=[0.5, 1.0, 2.0, 4.0])
    ///
    /// Returns:
    ///     List of PerResidueLddt containing:
    ///     - chain_id: Chain identifier
    ///     - residue_seq: Residue sequence number
    ///     - residue_name: Residue name
    ///     - score: LDDT score for this residue
    ///     - num_pairs: Number of distance pairs involving this residue
    ///
    /// Example:
    ///     >>> per_res = model.per_residue_lddt_to(reference)
    ///     >>> for r in per_res:
    ///     ...     if r.score < 0.7:
    ///     ...         print(f"{r.chain_id}{r.residue_seq}: LDDT = {r.score:.2f}")
    #[cfg(feature = "geometry")]
    #[pyo3(signature = (reference, selection=None, options=None))]
    fn per_residue_lddt_to(
        &self,
        reference: &PyPdbStructure,
        selection: Option<&crate::geometry::PyAtomSelection>,
        options: Option<&crate::geometry::PyLddtOptions>,
    ) -> PyResult<Vec<crate::geometry::PyPerResidueLddt>> {
        let sel = selection
            .map(|s| s.inner.clone())
            .unwrap_or(pdbrust::geometry::AtomSelection::CaOnly);
        let opts = options
            .map(|o| o.inner.clone())
            .unwrap_or_else(pdbrust::geometry::LddtOptions::default);
        self.inner
            .per_residue_lddt_to_with_options(&reference.inner, sel, opts)
            .map(|v| {
                v.into_iter()
                    .map(crate::geometry::PyPerResidueLddt::from)
                    .collect()
            })
            .map_err(convert_error)
    }

    // ==================== DSSP Methods (feature-gated) ====================

    /// Compute DSSP-like secondary structure assignment.
    ///
    /// This method analyzes the protein structure and assigns secondary
    /// structure to each residue using the DSSP algorithm.
    ///
    /// Returns:
    ///     SecondaryStructureAssignment containing:
    ///     - Per-residue assignments
    ///     - Statistics (helix/sheet/coil counts and fractions)
    ///     - Any warnings generated during assignment
    ///
    /// Example:
    ///     >>> ss = structure.assign_secondary_structure()
    ///     >>> print(f"Helix: {ss.helix_fraction * 100:.1f}%")
    ///     >>> print(f"Sheet: {ss.sheet_fraction * 100:.1f}%")
    ///     >>> for res in ss:
    ///     ...     print(f"{res.chain_id}{res.residue_seq}: {res.code()}")
    #[cfg(feature = "dssp")]
    fn assign_secondary_structure(&self) -> crate::dssp::PySecondaryStructureAssignment {
        crate::dssp::PySecondaryStructureAssignment::from(self.inner.assign_secondary_structure())
    }

    /// Get the secondary structure as a string of single-character codes.
    ///
    /// The string contains one character per residue in sequence order,
    /// using standard DSSP codes (H, G, I, P, E, B, T, S, C).
    ///
    /// Returns:
    ///     String like "CCCHHHHHHHHHCCCEEEEECCC"
    ///
    /// Example:
    ///     >>> ss_string = structure.secondary_structure_string()
    ///     >>> print(f"Secondary structure: {ss_string}")
    #[cfg(feature = "dssp")]
    fn secondary_structure_string(&self) -> String {
        self.inner.secondary_structure_string()
    }

    /// Get the secondary structure composition as fractions.
    ///
    /// Returns:
    ///     Tuple of (helix_fraction, sheet_fraction, coil_fraction) where:
    ///     - helix includes H, G, I, and P (κ-helix/PPII)
    ///     - sheet includes E and B
    ///     - coil includes T, S, and C
    ///     All fractions sum to 1.0.
    ///
    /// Example:
    ///     >>> helix, sheet, coil = structure.secondary_structure_composition()
    ///     >>> print(f"Helix: {helix * 100:.1f}%")
    ///     >>> print(f"Sheet: {sheet * 100:.1f}%")
    ///     >>> print(f"Coil:  {coil * 100:.1f}%")
    #[cfg(feature = "dssp")]
    fn secondary_structure_composition(&self) -> (f64, f64, f64) {
        self.inner.secondary_structure_composition()
    }

    // ==================== AlphaFold/pLDDT Methods (descriptors feature) ====================

    /// Check if this structure appears to be from AlphaFold/ESMFold.
    ///
    /// Detection is based on header/title keywords and B-factor distribution
    /// (which stores pLDDT scores in predicted structures).
    ///
    /// Returns:
    ///     True if structure appears to be AI-predicted
    ///
    /// Example:
    ///     >>> if structure.is_predicted():
    ///     ...     print("This is a predicted structure")
    #[cfg(feature = "descriptors")]
    fn is_predicted(&self) -> bool {
        self.inner.is_predicted()
    }

    /// Get mean pLDDT score across all atoms.
    ///
    /// Only meaningful for predicted structures (AlphaFold, ESMFold, etc.)
    /// where B-factors represent pLDDT confidence scores (0-100).
    ///
    /// Returns:
    ///     Mean pLDDT score
    ///
    /// Example:
    ///     >>> if structure.is_predicted():
    ///     ...     print(f"Mean pLDDT: {structure.plddt_mean():.1f}")
    #[cfg(feature = "descriptors")]
    fn plddt_mean(&self) -> f64 {
        self.inner.plddt_mean()
    }

    /// Get per-residue pLDDT scores.
    ///
    /// Returns:
    ///     List of ResiduePlddt objects with confidence scores
    ///
    /// Example:
    ///     >>> for res in structure.per_residue_plddt():
    ///     ...     print(f"{res.chain_id}{res.residue_seq}: {res.plddt:.1f}")
    #[cfg(feature = "descriptors")]
    fn per_residue_plddt(&self) -> Vec<crate::descriptors::PyResiduePlddt> {
        self.inner
            .per_residue_plddt()
            .into_iter()
            .map(crate::descriptors::PyResiduePlddt::from)
            .collect()
    }

    /// Get residues with low pLDDT confidence.
    ///
    /// Args:
    ///     threshold: pLDDT threshold (default: 70.0)
    ///
    /// Returns:
    ///     List of ResiduePlddt for residues below threshold
    ///
    /// Example:
    ///     >>> disordered = structure.low_confidence_regions(50.0)
    ///     >>> print(f"{len(disordered)} likely disordered residues")
    #[cfg(feature = "descriptors")]
    fn low_confidence_regions(&self, threshold: f64) -> Vec<crate::descriptors::PyResiduePlddt> {
        self.inner
            .low_confidence_regions(threshold)
            .into_iter()
            .map(crate::descriptors::PyResiduePlddt::from)
            .collect()
    }

    /// Get residues with high pLDDT confidence.
    ///
    /// Args:
    ///     threshold: pLDDT threshold (default: 70.0)
    ///
    /// Returns:
    ///     List of ResiduePlddt for residues at or above threshold
    #[cfg(feature = "descriptors")]
    fn high_confidence_regions(&self, threshold: f64) -> Vec<crate::descriptors::PyResiduePlddt> {
        self.inner
            .high_confidence_regions(threshold)
            .into_iter()
            .map(crate::descriptors::PyResiduePlddt::from)
            .collect()
    }

    /// Get fraction of residues in each pLDDT confidence category.
    ///
    /// Returns:
    ///     Tuple of (very_high, confident, low, very_low) fractions
    ///
    /// Example:
    ///     >>> vh, conf, low, vlow = structure.plddt_distribution()
    ///     >>> print(f"Very high (>90): {vh*100:.1f}%")
    ///     >>> print(f"Confident (70-90): {conf*100:.1f}%")
    #[cfg(feature = "descriptors")]
    fn plddt_distribution(&self) -> (f64, f64, f64, f64) {
        self.inner.plddt_distribution()
    }

    // ==================== Dihedral/Ramachandran Methods (dssp + descriptors) ====================

    /// Get backbone dihedral angles (phi, psi, omega) for all residues.
    ///
    /// Returns:
    ///     List of ResidueDihedrals with angles and Ramachandran classification
    ///
    /// Example:
    ///     >>> for d in structure.phi_psi_angles():
    ///     ...     if d.phi is not None:
    ///     ...         print(f"{d.residue_name}: phi={d.phi:.1f}, psi={d.psi:.1f}")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn phi_psi_angles(&self) -> Vec<crate::descriptors::PyResidueDihedrals> {
        self.inner
            .phi_psi_angles()
            .into_iter()
            .map(crate::descriptors::PyResidueDihedrals::from)
            .collect()
    }

    /// Get Ramachandran outliers.
    ///
    /// Returns:
    ///     List of ResidueDihedrals for residues in outlier regions
    ///
    /// Example:
    ///     >>> outliers = structure.ramachandran_outliers()
    ///     >>> print(f"Found {len(outliers)} Ramachandran outliers")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn ramachandran_outliers(&self) -> Vec<crate::descriptors::PyResidueDihedrals> {
        self.inner
            .ramachandran_outliers()
            .into_iter()
            .map(crate::descriptors::PyResidueDihedrals::from)
            .collect()
    }

    /// Detect cis peptide bonds in the structure.
    ///
    /// Returns:
    ///     List of (residue_i-1, residue_i) tuples as (chain, resid, resname) pairs
    ///
    /// Example:
    ///     >>> cis = structure.cis_peptide_bonds()
    ///     >>> for (r1, r2) in cis:
    ///     ...     print(f"Cis: {r1[2]}-{r2[2]}")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn cis_peptide_bonds(&self) -> Vec<((String, i32, String), (String, i32, String))> {
        self.inner
            .cis_peptide_bonds()
            .into_iter()
            .map(|(r1, r2)| {
                (
                    (r1.chain_id, r1.residue_seq, r1.residue_name),
                    (r2.chain_id, r2.residue_seq, r2.residue_name),
                )
            })
            .collect()
    }

    /// Get Ramachandran plot statistics for structure validation.
    ///
    /// Returns:
    ///     RamachandranStats with favored/allowed/outlier fractions
    ///
    /// Example:
    ///     >>> stats = structure.ramachandran_statistics()
    ///     >>> print(f"Favored: {stats.favored_fraction*100:.1f}%")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn ramachandran_statistics(&self) -> crate::descriptors::PyRamachandranStats {
        crate::descriptors::PyRamachandranStats::from(self.inner.ramachandran_statistics())
    }

    // ==================== H-Bond Methods (dssp + descriptors) ====================

    /// Get all mainchain (backbone) hydrogen bonds.
    ///
    /// Returns:
    ///     List of MainchainHBond with donor/acceptor info and energy
    ///
    /// Example:
    ///     >>> for hb in structure.mainchain_hbonds():
    ///     ...     print(f"{hb.donor_resid}->{hb.acceptor_resid}: {hb.energy:.2f} kcal/mol")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn mainchain_hbonds(&self) -> Vec<crate::descriptors::PyMainchainHBond> {
        self.inner
            .mainchain_hbonds()
            .into_iter()
            .map(crate::descriptors::PyMainchainHBond::from)
            .collect()
    }

    /// Get H-bonds for a specific residue.
    ///
    /// Args:
    ///     chain: Chain identifier
    ///     resid: Residue sequence number
    ///
    /// Returns:
    ///     ResidueHBonds with donated and accepted H-bonds
    ///
    /// Example:
    ///     >>> hbonds = structure.hbonds_for_residue("A", 42)
    ///     >>> print(f"Donates: {len(hbonds.donated)}, Accepts: {len(hbonds.accepted)}")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn hbonds_for_residue(&self, chain: &str, resid: i32) -> crate::descriptors::PyResidueHBonds {
        crate::descriptors::PyResidueHBonds::from(self.inner.hbonds_for_residue(chain, resid))
    }

    /// Get H-bond network statistics.
    ///
    /// Returns:
    ///     HBondStats with counts by type and energy statistics
    ///
    /// Example:
    ///     >>> stats = structure.hbond_statistics()
    ///     >>> print(f"Total: {stats.total_hbonds}, Helical: {stats.intra_helical}")
    #[cfg(all(feature = "descriptors", feature = "dssp"))]
    fn hbond_statistics(&self) -> crate::descriptors::PyHBondStats {
        crate::descriptors::PyHBondStats::from(self.inner.hbond_statistics())
    }

    // ==================== Protein-Ligand Interaction Methods (descriptors) ====================

    /// Get the binding site around a ligand.
    ///
    /// Args:
    ///     ligand_name: 3-letter code of the ligand (e.g., "ATP", "HEM")
    ///     distance_cutoff: Maximum distance in Angstroms (typically 4.0-6.0)
    ///
    /// Returns:
    ///     BindingSite with contact residues, or None if ligand not found
    ///
    /// Example:
    ///     >>> site = structure.binding_site("ATP", 5.0)
    ///     >>> if site:
    ///     ...     print(f"Found {site.num_residues()} residues in binding site")
    #[cfg(feature = "descriptors")]
    fn binding_site(
        &self,
        ligand_name: &str,
        distance_cutoff: f64,
    ) -> Option<crate::descriptors::PyBindingSite> {
        self.inner
            .binding_site(ligand_name, distance_cutoff)
            .map(crate::descriptors::PyBindingSite::from)
    }

    /// Analyze interactions between protein and a specific ligand.
    ///
    /// Detects H-bonds, salt bridges, and hydrophobic contacts.
    ///
    /// Args:
    ///     ligand_name: 3-letter code of the ligand
    ///
    /// Returns:
    ///     LigandInteractionProfile with all interactions, or None if not found
    ///
    /// Example:
    ///     >>> profile = structure.ligand_interactions("ATP")
    ///     >>> if profile:
    ///     ...     print(f"H-bonds: {len(profile.hydrogen_bonds)}")
    ///     ...     print(f"Salt bridges: {len(profile.salt_bridges)}")
    #[cfg(feature = "descriptors")]
    fn ligand_interactions(
        &self,
        ligand_name: &str,
    ) -> Option<crate::descriptors::PyLigandInteractionProfile> {
        self.inner
            .ligand_interactions(ligand_name)
            .map(crate::descriptors::PyLigandInteractionProfile::from)
    }

    /// Get interaction profiles for all ligands in the structure.
    ///
    /// Returns:
    ///     List of LigandInteractionProfile for each ligand
    ///
    /// Example:
    ///     >>> for profile in structure.all_ligand_interactions():
    ///     ...     print(f"{profile.ligand_name}: {profile.total_interactions()} interactions")
    #[cfg(feature = "descriptors")]
    fn all_ligand_interactions(&self) -> Vec<crate::descriptors::PyLigandInteractionProfile> {
        self.inner
            .all_ligand_interactions()
            .into_iter()
            .map(crate::descriptors::PyLigandInteractionProfile::from)
            .collect()
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
