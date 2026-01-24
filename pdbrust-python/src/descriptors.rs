//! Python bindings for structure descriptors

use pdbrust::descriptors::{
    BindingSite, ConfidenceCategory, ContactResidue, HydrophobicContact, LigandInteractionProfile,
    ProteinLigandHBond, ResidueBFactor, ResiduePlddt, SaltBridge, StructureDescriptors,
};
use pyo3::prelude::*;
use std::collections::HashMap;

#[cfg(feature = "dssp")]
use pdbrust::descriptors::{
    HBondStats, HBondType, MainchainHBond, RamachandranRegion, RamachandranStats, ResidueDihedrals,
    ResidueHBonds,
};

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

// ============================================================================
// AlphaFold/pLDDT Support
// ============================================================================

/// pLDDT confidence category
#[pyclass(name = "ConfidenceCategory")]
#[derive(Clone)]
pub struct PyConfidenceCategory {
    inner: ConfidenceCategory,
}

#[pymethods]
impl PyConfidenceCategory {
    /// Returns true if this is a reliable region (VeryHigh or Confident)
    fn is_reliable(&self) -> bool {
        self.inner.is_reliable()
    }

    /// Returns true if this region should be treated with caution
    fn needs_caution(&self) -> bool {
        self.inner.needs_caution()
    }

    fn __repr__(&self) -> String {
        match self.inner {
            ConfidenceCategory::VeryHigh => "ConfidenceCategory.VeryHigh".to_string(),
            ConfidenceCategory::Confident => "ConfidenceCategory.Confident".to_string(),
            ConfidenceCategory::Low => "ConfidenceCategory.Low".to_string(),
            ConfidenceCategory::VeryLow => "ConfidenceCategory.VeryLow".to_string(),
        }
    }

    fn __str__(&self) -> String {
        match self.inner {
            ConfidenceCategory::VeryHigh => "VeryHigh (>90)".to_string(),
            ConfidenceCategory::Confident => "Confident (70-90)".to_string(),
            ConfidenceCategory::Low => "Low (50-70)".to_string(),
            ConfidenceCategory::VeryLow => "VeryLow (<50)".to_string(),
        }
    }
}

impl From<ConfidenceCategory> for PyConfidenceCategory {
    fn from(cat: ConfidenceCategory) -> Self {
        PyConfidenceCategory { inner: cat }
    }
}

/// Per-residue pLDDT confidence score
#[pyclass(name = "ResiduePlddt")]
#[derive(Clone)]
pub struct PyResiduePlddt {
    pub(crate) inner: ResiduePlddt,
}

#[pymethods]
impl PyResiduePlddt {
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

    /// Residue name
    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    /// Mean pLDDT score
    #[getter]
    fn plddt(&self) -> f64 {
        self.inner.plddt
    }

    /// Minimum pLDDT among atoms
    #[getter]
    fn plddt_min(&self) -> f64 {
        self.inner.plddt_min
    }

    /// Maximum pLDDT among atoms
    #[getter]
    fn plddt_max(&self) -> f64 {
        self.inner.plddt_max
    }

    /// Number of atoms
    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count
    }

    /// Confidence category
    #[getter]
    fn confidence_category(&self) -> PyConfidenceCategory {
        PyConfidenceCategory::from(self.inner.confidence_category)
    }

    /// Returns true if this residue has high confidence
    fn is_confident(&self) -> bool {
        self.inner.is_confident()
    }

    /// Returns true if this residue is likely disordered
    fn is_disordered(&self) -> bool {
        self.inner.is_disordered()
    }

    fn __repr__(&self) -> String {
        format!(
            "ResiduePlddt(chain={}, seq={}, name={}, plddt={:.1})",
            self.inner.chain_id, self.inner.residue_seq, self.inner.residue_name, self.inner.plddt
        )
    }
}

impl From<ResiduePlddt> for PyResiduePlddt {
    fn from(res: ResiduePlddt) -> Self {
        PyResiduePlddt { inner: res }
    }
}

// ============================================================================
// Protein-Ligand Interaction Types
// ============================================================================

/// A residue in contact with a ligand
#[pyclass(name = "ContactResidue")]
#[derive(Clone)]
pub struct PyContactResidue {
    pub(crate) inner: ContactResidue,
}

#[pymethods]
impl PyContactResidue {
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    #[getter]
    fn min_distance(&self) -> f64 {
        self.inner.min_distance
    }

    #[getter]
    fn num_contacts(&self) -> usize {
        self.inner.num_contacts
    }

    fn __repr__(&self) -> String {
        format!(
            "ContactResidue({}{} {}, dist={:.2})",
            self.inner.chain_id,
            self.inner.residue_seq,
            self.inner.residue_name,
            self.inner.min_distance
        )
    }
}

impl From<ContactResidue> for PyContactResidue {
    fn from(res: ContactResidue) -> Self {
        PyContactResidue { inner: res }
    }
}

/// A hydrogen bond between protein and ligand
#[pyclass(name = "ProteinLigandHBond")]
#[derive(Clone)]
pub struct PyProteinLigandHBond {
    pub(crate) inner: ProteinLigandHBond,
}

#[pymethods]
impl PyProteinLigandHBond {
    #[getter]
    fn protein_chain(&self) -> &str {
        &self.inner.protein_chain
    }

    #[getter]
    fn protein_resid(&self) -> i32 {
        self.inner.protein_resid
    }

    #[getter]
    fn protein_resname(&self) -> &str {
        &self.inner.protein_resname
    }

    #[getter]
    fn protein_atom(&self) -> &str {
        &self.inner.protein_atom
    }

    #[getter]
    fn ligand_name(&self) -> &str {
        &self.inner.ligand_name
    }

    #[getter]
    fn ligand_atom(&self) -> &str {
        &self.inner.ligand_atom
    }

    #[getter]
    fn distance(&self) -> f64 {
        self.inner.distance
    }

    #[getter]
    fn is_protein_donor(&self) -> bool {
        self.inner.is_protein_donor
    }

    fn __repr__(&self) -> String {
        format!(
            "ProteinLigandHBond({}:{}-{}:{}, dist={:.2})",
            self.inner.protein_resname,
            self.inner.protein_atom,
            self.inner.ligand_name,
            self.inner.ligand_atom,
            self.inner.distance
        )
    }
}

impl From<ProteinLigandHBond> for PyProteinLigandHBond {
    fn from(hb: ProteinLigandHBond) -> Self {
        PyProteinLigandHBond { inner: hb }
    }
}

/// A salt bridge between protein and ligand
#[pyclass(name = "SaltBridge")]
#[derive(Clone)]
pub struct PySaltBridge {
    pub(crate) inner: SaltBridge,
}

#[pymethods]
impl PySaltBridge {
    #[getter]
    fn protein_chain(&self) -> &str {
        &self.inner.protein_chain
    }

    #[getter]
    fn protein_resid(&self) -> i32 {
        self.inner.protein_resid
    }

    #[getter]
    fn protein_resname(&self) -> &str {
        &self.inner.protein_resname
    }

    #[getter]
    fn protein_atom(&self) -> &str {
        &self.inner.protein_atom
    }

    #[getter]
    fn ligand_name(&self) -> &str {
        &self.inner.ligand_name
    }

    #[getter]
    fn ligand_atom(&self) -> &str {
        &self.inner.ligand_atom
    }

    #[getter]
    fn distance(&self) -> f64 {
        self.inner.distance
    }

    #[getter]
    fn protein_positive(&self) -> bool {
        self.inner.protein_positive
    }

    fn __repr__(&self) -> String {
        format!(
            "SaltBridge({}:{}-{}:{}, dist={:.2})",
            self.inner.protein_resname,
            self.inner.protein_atom,
            self.inner.ligand_name,
            self.inner.ligand_atom,
            self.inner.distance
        )
    }
}

impl From<SaltBridge> for PySaltBridge {
    fn from(sb: SaltBridge) -> Self {
        PySaltBridge { inner: sb }
    }
}

/// A hydrophobic contact between protein and ligand
#[pyclass(name = "HydrophobicContact")]
#[derive(Clone)]
pub struct PyHydrophobicContact {
    pub(crate) inner: HydrophobicContact,
}

#[pymethods]
impl PyHydrophobicContact {
    #[getter]
    fn protein_chain(&self) -> &str {
        &self.inner.protein_chain
    }

    #[getter]
    fn protein_resid(&self) -> i32 {
        self.inner.protein_resid
    }

    #[getter]
    fn protein_resname(&self) -> &str {
        &self.inner.protein_resname
    }

    #[getter]
    fn protein_atom(&self) -> &str {
        &self.inner.protein_atom
    }

    #[getter]
    fn ligand_name(&self) -> &str {
        &self.inner.ligand_name
    }

    #[getter]
    fn ligand_atom(&self) -> &str {
        &self.inner.ligand_atom
    }

    #[getter]
    fn distance(&self) -> f64 {
        self.inner.distance
    }

    fn __repr__(&self) -> String {
        format!(
            "HydrophobicContact({}:{}-{}:{}, dist={:.2})",
            self.inner.protein_resname,
            self.inner.protein_atom,
            self.inner.ligand_name,
            self.inner.ligand_atom,
            self.inner.distance
        )
    }
}

impl From<HydrophobicContact> for PyHydrophobicContact {
    fn from(hc: HydrophobicContact) -> Self {
        PyHydrophobicContact { inner: hc }
    }
}

/// Ligand interaction profile
#[pyclass(name = "LigandInteractionProfile")]
#[derive(Clone)]
pub struct PyLigandInteractionProfile {
    pub(crate) inner: LigandInteractionProfile,
}

#[pymethods]
impl PyLigandInteractionProfile {
    #[getter]
    fn ligand_name(&self) -> &str {
        &self.inner.ligand_name
    }

    #[getter]
    fn ligand_chain(&self) -> &str {
        &self.inner.ligand_chain
    }

    #[getter]
    fn ligand_resid(&self) -> i32 {
        self.inner.ligand_resid
    }

    #[getter]
    fn contact_residues(&self) -> Vec<PyContactResidue> {
        self.inner
            .contact_residues
            .iter()
            .map(|r| PyContactResidue::from(r.clone()))
            .collect()
    }

    #[getter]
    fn hydrogen_bonds(&self) -> Vec<PyProteinLigandHBond> {
        self.inner
            .hydrogen_bonds
            .iter()
            .map(|h| PyProteinLigandHBond::from(h.clone()))
            .collect()
    }

    #[getter]
    fn salt_bridges(&self) -> Vec<PySaltBridge> {
        self.inner
            .salt_bridges
            .iter()
            .map(|s| PySaltBridge::from(s.clone()))
            .collect()
    }

    #[getter]
    fn hydrophobic_contacts(&self) -> Vec<PyHydrophobicContact> {
        self.inner
            .hydrophobic_contacts
            .iter()
            .map(|h| PyHydrophobicContact::from(h.clone()))
            .collect()
    }

    /// Total number of interactions
    fn total_interactions(&self) -> usize {
        self.inner.total_interactions()
    }

    /// Returns true if any interactions were detected
    fn has_interactions(&self) -> bool {
        self.inner.has_interactions()
    }

    fn __repr__(&self) -> String {
        format!(
            "LigandInteractionProfile({}, contacts={}, hbonds={}, salt={}, hydrophobic={})",
            self.inner.ligand_name,
            self.inner.contact_residues.len(),
            self.inner.hydrogen_bonds.len(),
            self.inner.salt_bridges.len(),
            self.inner.hydrophobic_contacts.len()
        )
    }
}

impl From<LigandInteractionProfile> for PyLigandInteractionProfile {
    fn from(profile: LigandInteractionProfile) -> Self {
        PyLigandInteractionProfile { inner: profile }
    }
}

/// Binding site around a ligand
#[pyclass(name = "BindingSite")]
#[derive(Clone)]
pub struct PyBindingSite {
    pub(crate) inner: BindingSite,
}

#[pymethods]
impl PyBindingSite {
    #[getter]
    fn ligand_name(&self) -> &str {
        &self.inner.ligand_name
    }

    #[getter]
    fn ligand_chain(&self) -> &str {
        &self.inner.ligand_chain
    }

    #[getter]
    fn ligand_resid(&self) -> i32 {
        self.inner.ligand_resid
    }

    #[getter]
    fn distance_cutoff(&self) -> f64 {
        self.inner.distance_cutoff
    }

    #[getter]
    fn contact_residues(&self) -> Vec<PyContactResidue> {
        self.inner
            .contact_residues
            .iter()
            .map(|r| PyContactResidue::from(r.clone()))
            .collect()
    }

    /// Number of contact residues
    fn num_residues(&self) -> usize {
        self.inner.num_residues()
    }

    fn __repr__(&self) -> String {
        format!(
            "BindingSite({}, {} residues within {} Å)",
            self.inner.ligand_name,
            self.inner.num_residues(),
            self.inner.distance_cutoff
        )
    }
}

impl From<BindingSite> for PyBindingSite {
    fn from(site: BindingSite) -> Self {
        PyBindingSite { inner: site }
    }
}

// ============================================================================
// Dihedral/Ramachandran Types (requires dssp feature)
// ============================================================================

#[cfg(feature = "dssp")]
/// Ramachandran region classification
#[pyclass(name = "RamachandranRegion")]
#[derive(Clone)]
pub struct PyRamachandranRegion {
    inner: RamachandranRegion,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyRamachandranRegion {
    /// Returns true if this is a favorable region
    fn is_favorable(&self) -> bool {
        self.inner.is_favorable()
    }

    /// Returns true if this is an outlier
    fn is_outlier(&self) -> bool {
        self.inner.is_outlier()
    }

    fn __repr__(&self) -> String {
        match self.inner {
            RamachandranRegion::Core => "RamachandranRegion.Core".to_string(),
            RamachandranRegion::Allowed => "RamachandranRegion.Allowed".to_string(),
            RamachandranRegion::Generous => "RamachandranRegion.Generous".to_string(),
            RamachandranRegion::Outlier => "RamachandranRegion.Outlier".to_string(),
            RamachandranRegion::Glycine => "RamachandranRegion.Glycine".to_string(),
            RamachandranRegion::Proline => "RamachandranRegion.Proline".to_string(),
            RamachandranRegion::PrePro => "RamachandranRegion.PrePro".to_string(),
            RamachandranRegion::Unknown => "RamachandranRegion.Unknown".to_string(),
        }
    }
}

#[cfg(feature = "dssp")]
impl From<RamachandranRegion> for PyRamachandranRegion {
    fn from(region: RamachandranRegion) -> Self {
        PyRamachandranRegion { inner: region }
    }
}

#[cfg(feature = "dssp")]
/// Per-residue dihedral angles
#[pyclass(name = "ResidueDihedrals")]
#[derive(Clone)]
pub struct PyResidueDihedrals {
    pub(crate) inner: ResidueDihedrals,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyResidueDihedrals {
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    #[getter]
    fn residue_name(&self) -> &str {
        &self.inner.residue_name
    }

    #[getter]
    fn phi(&self) -> Option<f64> {
        self.inner.phi
    }

    #[getter]
    fn psi(&self) -> Option<f64> {
        self.inner.psi
    }

    #[getter]
    fn omega(&self) -> Option<f64> {
        self.inner.omega
    }

    #[getter]
    fn ramachandran_region(&self) -> PyRamachandranRegion {
        PyRamachandranRegion::from(self.inner.ramachandran_region)
    }

    /// Returns true if this residue has both phi and psi
    fn has_phi_psi(&self) -> bool {
        self.inner.has_phi_psi()
    }

    /// Returns true if this is a cis peptide bond
    fn is_cis_peptide(&self) -> bool {
        self.inner.is_cis_peptide()
    }

    /// Returns true if this is a trans peptide bond
    fn is_trans_peptide(&self) -> bool {
        self.inner.is_trans_peptide()
    }

    fn __repr__(&self) -> String {
        format!(
            "ResidueDihedrals({}{}, phi={:?}, psi={:?})",
            self.inner.chain_id, self.inner.residue_seq, self.inner.phi, self.inner.psi
        )
    }
}

#[cfg(feature = "dssp")]
impl From<ResidueDihedrals> for PyResidueDihedrals {
    fn from(d: ResidueDihedrals) -> Self {
        PyResidueDihedrals { inner: d }
    }
}

#[cfg(feature = "dssp")]
/// Ramachandran statistics
#[pyclass(name = "RamachandranStats")]
#[derive(Clone)]
pub struct PyRamachandranStats {
    pub(crate) inner: RamachandranStats,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyRamachandranStats {
    #[getter]
    fn total_residues(&self) -> usize {
        self.inner.total_residues
    }

    #[getter]
    fn favored_count(&self) -> usize {
        self.inner.favored_count
    }

    #[getter]
    fn allowed_count(&self) -> usize {
        self.inner.allowed_count
    }

    #[getter]
    fn outlier_count(&self) -> usize {
        self.inner.outlier_count
    }

    #[getter]
    fn favored_fraction(&self) -> f64 {
        self.inner.favored_fraction
    }

    #[getter]
    fn allowed_fraction(&self) -> f64 {
        self.inner.allowed_fraction
    }

    #[getter]
    fn outlier_fraction(&self) -> f64 {
        self.inner.outlier_fraction
    }

    #[getter]
    fn cis_peptide_count(&self) -> usize {
        self.inner.cis_peptide_count
    }

    #[getter]
    fn cis_nonpro_count(&self) -> usize {
        self.inner.cis_nonpro_count
    }

    fn __repr__(&self) -> String {
        format!(
            "RamachandranStats(favored={:.1}%, allowed={:.1}%, outliers={:.1}%)",
            self.inner.favored_fraction * 100.0,
            self.inner.allowed_fraction * 100.0,
            self.inner.outlier_fraction * 100.0
        )
    }
}

#[cfg(feature = "dssp")]
impl From<RamachandranStats> for PyRamachandranStats {
    fn from(stats: RamachandranStats) -> Self {
        PyRamachandranStats { inner: stats }
    }
}

// ============================================================================
// Hydrogen Bond Types (requires dssp feature)
// ============================================================================

#[cfg(feature = "dssp")]
/// H-bond type classification
#[pyclass(name = "HBondType")]
#[derive(Clone)]
pub struct PyHBondType {
    inner: HBondType,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyHBondType {
    /// Returns true if this is a secondary structure-forming H-bond
    fn is_secondary_structure(&self) -> bool {
        self.inner.is_secondary_structure()
    }

    fn __repr__(&self) -> String {
        match self.inner {
            HBondType::IntraHelical => "HBondType.IntraHelical".to_string(),
            HBondType::BetaSheet => "HBondType.BetaSheet".to_string(),
            HBondType::Turn => "HBondType.Turn".to_string(),
            HBondType::LongRange => "HBondType.LongRange".to_string(),
            HBondType::InterChain => "HBondType.InterChain".to_string(),
        }
    }
}

#[cfg(feature = "dssp")]
impl From<HBondType> for PyHBondType {
    fn from(t: HBondType) -> Self {
        PyHBondType { inner: t }
    }
}

#[cfg(feature = "dssp")]
/// A mainchain hydrogen bond
#[pyclass(name = "MainchainHBond")]
#[derive(Clone)]
pub struct PyMainchainHBond {
    pub(crate) inner: MainchainHBond,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyMainchainHBond {
    #[getter]
    fn donor_chain(&self) -> &str {
        &self.inner.donor_chain
    }

    #[getter]
    fn donor_resid(&self) -> i32 {
        self.inner.donor_resid
    }

    #[getter]
    fn donor_resname(&self) -> &str {
        &self.inner.donor_resname
    }

    #[getter]
    fn acceptor_chain(&self) -> &str {
        &self.inner.acceptor_chain
    }

    #[getter]
    fn acceptor_resid(&self) -> i32 {
        self.inner.acceptor_resid
    }

    #[getter]
    fn acceptor_resname(&self) -> &str {
        &self.inner.acceptor_resname
    }

    #[getter]
    fn energy(&self) -> f64 {
        self.inner.energy
    }

    #[getter]
    fn n_o_distance(&self) -> f64 {
        self.inner.n_o_distance
    }

    #[getter]
    fn sequence_separation(&self) -> i32 {
        self.inner.sequence_separation
    }

    #[getter]
    fn hbond_type(&self) -> PyHBondType {
        PyHBondType::from(self.inner.hbond_type)
    }

    /// Returns true if this is a strong H-bond
    fn is_strong(&self) -> bool {
        self.inner.is_strong()
    }

    /// Returns true if this is a helical H-bond
    fn is_helical(&self) -> bool {
        self.inner.is_helical()
    }

    /// Returns true if this is a beta-sheet H-bond
    fn is_beta_sheet(&self) -> bool {
        self.inner.is_beta_sheet()
    }

    fn __repr__(&self) -> String {
        format!(
            "MainchainHBond({}{}->{}{}, E={:.2})",
            self.inner.donor_chain,
            self.inner.donor_resid,
            self.inner.acceptor_chain,
            self.inner.acceptor_resid,
            self.inner.energy
        )
    }
}

#[cfg(feature = "dssp")]
impl From<MainchainHBond> for PyMainchainHBond {
    fn from(hb: MainchainHBond) -> Self {
        PyMainchainHBond { inner: hb }
    }
}

#[cfg(feature = "dssp")]
/// H-bonds for a specific residue
#[pyclass(name = "ResidueHBonds")]
#[derive(Clone)]
pub struct PyResidueHBonds {
    pub(crate) inner: ResidueHBonds,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyResidueHBonds {
    #[getter]
    fn donated(&self) -> Vec<PyMainchainHBond> {
        self.inner
            .donated
            .iter()
            .map(|h| PyMainchainHBond::from(h.clone()))
            .collect()
    }

    #[getter]
    fn accepted(&self) -> Vec<PyMainchainHBond> {
        self.inner
            .accepted
            .iter()
            .map(|h| PyMainchainHBond::from(h.clone()))
            .collect()
    }

    /// Total number of H-bonds
    fn total(&self) -> usize {
        self.inner.total()
    }

    /// Returns true if this residue has any H-bonds
    fn has_hbonds(&self) -> bool {
        self.inner.has_hbonds()
    }

    fn __repr__(&self) -> String {
        format!(
            "ResidueHBonds(donated={}, accepted={})",
            self.inner.donated.len(),
            self.inner.accepted.len()
        )
    }
}

#[cfg(feature = "dssp")]
impl From<ResidueHBonds> for PyResidueHBonds {
    fn from(hb: ResidueHBonds) -> Self {
        PyResidueHBonds { inner: hb }
    }
}

#[cfg(feature = "dssp")]
/// H-bond network statistics
#[pyclass(name = "HBondStats")]
#[derive(Clone)]
pub struct PyHBondStats {
    pub(crate) inner: HBondStats,
}

#[cfg(feature = "dssp")]
#[pymethods]
impl PyHBondStats {
    #[getter]
    fn total_hbonds(&self) -> usize {
        self.inner.total_hbonds
    }

    #[getter]
    fn intra_helical(&self) -> usize {
        self.inner.intra_helical
    }

    #[getter]
    fn beta_sheet(&self) -> usize {
        self.inner.beta_sheet
    }

    #[getter]
    fn turn(&self) -> usize {
        self.inner.turn
    }

    #[getter]
    fn long_range(&self) -> usize {
        self.inner.long_range
    }

    #[getter]
    fn inter_chain(&self) -> usize {
        self.inner.inter_chain
    }

    #[getter]
    fn mean_energy(&self) -> f64 {
        self.inner.mean_energy
    }

    #[getter]
    fn donor_residues(&self) -> usize {
        self.inner.donor_residues
    }

    #[getter]
    fn acceptor_residues(&self) -> usize {
        self.inner.acceptor_residues
    }

    fn __repr__(&self) -> String {
        format!(
            "HBondStats(total={}, helical={}, sheet={}, mean_E={:.2})",
            self.inner.total_hbonds,
            self.inner.intra_helical,
            self.inner.beta_sheet,
            self.inner.mean_energy
        )
    }
}

#[cfg(feature = "dssp")]
impl From<HBondStats> for PyHBondStats {
    fn from(stats: HBondStats) -> Self {
        PyHBondStats { inner: stats }
    }
}
