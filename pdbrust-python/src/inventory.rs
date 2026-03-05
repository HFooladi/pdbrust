//! Python bindings for molecular inventory

use pdbrust::inventory::{ChainInventory, ChainType, LigandInfo, MolecularInventory};
use pyo3::prelude::*;

/// Type of molecular content in a chain
#[pyclass(name = "ChainType")]
#[derive(Clone)]
pub struct PyChainType {
    inner: ChainType,
}

#[pymethods]
impl PyChainType {
    /// String representation (e.g. "Protein", "Nucleic acid")
    fn __repr__(&self) -> String {
        format!("ChainType({})", self.inner)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Whether this is a protein chain
    #[getter]
    fn is_protein(&self) -> bool {
        self.inner == ChainType::Protein
    }

    /// Whether this is a nucleic acid chain
    #[getter]
    fn is_nucleic(&self) -> bool {
        self.inner == ChainType::NucleicAcid
    }

    /// Whether this is a mixed chain
    #[getter]
    fn is_mixed(&self) -> bool {
        self.inner == ChainType::Mixed
    }

    /// Whether this is a water-only chain
    #[getter]
    fn is_water(&self) -> bool {
        self.inner == ChainType::Water
    }
}

/// Per-chain summary
#[pyclass(name = "ChainInventory")]
#[derive(Clone)]
pub struct PyChainInventory {
    inner: ChainInventory,
}

#[pymethods]
impl PyChainInventory {
    /// Chain identifier
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Dominant molecular type
    #[getter]
    fn chain_type(&self) -> PyChainType {
        PyChainType {
            inner: self.inner.chain_type,
        }
    }

    /// Total atoms in this chain
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    /// Number of standard amino-acid residues
    #[getter]
    fn protein_residues(&self) -> usize {
        self.inner.protein_residues
    }

    /// Number of standard nucleotide residues
    #[getter]
    fn nucleic_residues(&self) -> usize {
        self.inner.nucleic_residues
    }

    /// Number of water molecules
    #[getter]
    fn water_molecules(&self) -> usize {
        self.inner.water_molecules
    }

    /// Number of non-water HETATM residues
    #[getter]
    fn het_residues(&self) -> usize {
        self.inner.het_residues
    }

    fn __repr__(&self) -> String {
        format!(
            "ChainInventory(chain={}, type={}, atoms={})",
            self.inner.chain_id, self.inner.chain_type, self.inner.num_atoms
        )
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

/// Description of a single ligand instance
#[pyclass(name = "LigandInfo")]
#[derive(Clone)]
pub struct PyLigandInfo {
    inner: LigandInfo,
}

#[pymethods]
impl PyLigandInfo {
    /// Residue name (e.g. "ATP")
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Chain where the ligand resides
    #[getter]
    fn chain_id(&self) -> &str {
        &self.inner.chain_id
    }

    /// Residue sequence number
    #[getter]
    fn residue_seq(&self) -> i32 {
        self.inner.residue_seq
    }

    /// Number of atoms
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms
    }

    fn __repr__(&self) -> String {
        format!(
            "LigandInfo(name={}, chain={}, resid={}, atoms={})",
            self.inner.name, self.inner.chain_id, self.inner.residue_seq, self.inner.num_atoms
        )
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

/// Complete molecular inventory of a PDB structure
///
/// Provides a one-call breakdown of all molecular entities:
/// protein chains, nucleic acid chains, ligands, water, and ions.
///
/// Example:
///     >>> structure = pdbrust.parse_pdb_file("complex.pdb")
///     >>> inv = structure.molecular_inventory()
///     >>> print(inv)
///     >>> print(f"Protein chains: {inv.num_protein_chains}")
///     >>> for lig in inv.ligands:
///     ...     print(f"  {lig.name} in chain {lig.chain_id}")
#[pyclass(name = "MolecularInventory")]
#[derive(Clone)]
pub struct PyMolecularInventory {
    inner: MolecularInventory,
}

#[pymethods]
impl PyMolecularInventory {
    /// Per-chain breakdown
    #[getter]
    fn chains(&self) -> Vec<PyChainInventory> {
        self.inner
            .chains
            .iter()
            .map(|c| PyChainInventory { inner: c.clone() })
            .collect()
    }

    /// All ligand instances
    #[getter]
    fn ligands(&self) -> Vec<PyLigandInfo> {
        self.inner
            .ligands
            .iter()
            .map(|l| PyLigandInfo { inner: l.clone() })
            .collect()
    }

    /// Total number of atoms
    #[getter]
    fn total_atoms(&self) -> usize {
        self.inner.total_atoms
    }

    /// Total protein atoms
    #[getter]
    fn total_protein_atoms(&self) -> usize {
        self.inner.total_protein_atoms
    }

    /// Total nucleic acid atoms
    #[getter]
    fn total_nucleic_atoms(&self) -> usize {
        self.inner.total_nucleic_atoms
    }

    /// Total water atoms
    #[getter]
    fn total_water_atoms(&self) -> usize {
        self.inner.total_water_atoms
    }

    /// Total non-water HETATM atoms
    #[getter]
    fn total_het_atoms(&self) -> usize {
        self.inner.total_het_atoms
    }

    /// Number of protein chains
    #[getter]
    fn num_protein_chains(&self) -> usize {
        self.inner.num_protein_chains
    }

    /// Number of nucleic acid chains
    #[getter]
    fn num_nucleic_chains(&self) -> usize {
        self.inner.num_nucleic_chains
    }

    /// Total water molecules
    #[getter]
    fn num_water_molecules(&self) -> usize {
        self.inner.num_water_molecules
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularInventory(atoms={}, protein_chains={}, nucleic_chains={}, ligands={}, waters={})",
            self.inner.total_atoms,
            self.inner.num_protein_chains,
            self.inner.num_nucleic_chains,
            self.inner.ligands.len(),
            self.inner.num_water_molecules
        )
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

impl From<MolecularInventory> for PyMolecularInventory {
    fn from(inv: MolecularInventory) -> Self {
        PyMolecularInventory { inner: inv }
    }
}
