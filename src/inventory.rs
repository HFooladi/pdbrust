//! Molecular inventory: a one-call breakdown of structure contents.
//!
//! Provides [`MolecularInventory`] which summarises what is inside a PDB
//! structure — protein chains, nucleic acid chains, ligands, water molecules,
//! and ions — without requiring any optional feature flags.
//!
//! # Example
//!
//! ```ignore
//! use pdbrust::parse_pdb_file;
//!
//! let structure = parse_pdb_file("complex.pdb")?;
//! let inv = structure.molecular_inventory();
//!
//! println!("{}", inv);
//! // Structure inventory:
//! //   Chain A: Protein (99 residues, 757 atoms)
//! //   Chain B: Protein (99 residues, 757 atoms)
//! //   Ligands: MK1 (chain B, res 902, 45 atoms)
//! //   Water: 127 molecules
//! //   Total: 1686 atoms
//! ```

use crate::classify::{COMMON_IONS, is_standard_amino_acid, is_standard_nucleotide, is_water};
use crate::core::PdbStructure;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// The dominant molecular type in a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChainType {
    /// Mostly standard amino acids.
    Protein,
    /// Mostly standard nucleotides.
    NucleicAcid,
    /// Contains both protein and nucleic residues.
    Mixed,
    /// Only water molecules.
    Water,
    /// Only ions or other HETATM.
    Other,
}

impl fmt::Display for ChainType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChainType::Protein => write!(f, "Protein"),
            ChainType::NucleicAcid => write!(f, "Nucleic acid"),
            ChainType::Mixed => write!(f, "Mixed"),
            ChainType::Water => write!(f, "Water"),
            ChainType::Other => write!(f, "Other"),
        }
    }
}

/// Per-chain summary.
#[derive(Debug, Clone)]
pub struct ChainInventory {
    /// Chain identifier.
    pub chain_id: String,
    /// Dominant molecular type.
    pub chain_type: ChainType,
    /// Total atoms in this chain.
    pub num_atoms: usize,
    /// Number of standard amino-acid residues.
    pub protein_residues: usize,
    /// Number of standard nucleotide residues.
    pub nucleic_residues: usize,
    /// Number of water molecules.
    pub water_molecules: usize,
    /// Number of non-water HETATM residues (ligands, ions, etc.).
    pub het_residues: usize,
}

impl fmt::Display for ChainInventory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.chain_type {
            ChainType::Protein => write!(
                f,
                "Chain {}: Protein ({} residues, {} atoms)",
                self.chain_id, self.protein_residues, self.num_atoms
            ),
            ChainType::NucleicAcid => write!(
                f,
                "Chain {}: Nucleic acid ({} residues, {} atoms)",
                self.chain_id, self.nucleic_residues, self.num_atoms
            ),
            ChainType::Mixed => write!(
                f,
                "Chain {}: Mixed ({} protein + {} nucleic residues, {} atoms)",
                self.chain_id, self.protein_residues, self.nucleic_residues, self.num_atoms
            ),
            ChainType::Water => write!(
                f,
                "Chain {}: Water ({} molecules)",
                self.chain_id, self.water_molecules
            ),
            ChainType::Other => write!(
                f,
                "Chain {}: Other ({} residues, {} atoms)",
                self.chain_id, self.het_residues, self.num_atoms
            ),
        }
    }
}

/// Description of a single ligand instance.
#[derive(Debug, Clone)]
pub struct LigandInfo {
    /// Residue name (3-letter code, e.g. "ATP").
    pub name: String,
    /// Chain where the ligand resides.
    pub chain_id: String,
    /// Residue sequence number.
    pub residue_seq: i32,
    /// Number of atoms in this ligand instance.
    pub num_atoms: usize,
}

impl fmt::Display for LigandInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} (chain {}, res {}, {} atoms)",
            self.name, self.chain_id, self.residue_seq, self.num_atoms
        )
    }
}

/// Complete molecular inventory of a PDB structure.
#[derive(Debug, Clone)]
pub struct MolecularInventory {
    /// Per-chain breakdown (ordered by first appearance).
    pub chains: Vec<ChainInventory>,
    /// All ligand instances found.
    pub ligands: Vec<LigandInfo>,
    /// Total number of atoms.
    pub total_atoms: usize,
    /// Total protein atoms.
    pub total_protein_atoms: usize,
    /// Total nucleic acid atoms.
    pub total_nucleic_atoms: usize,
    /// Total water atoms.
    pub total_water_atoms: usize,
    /// Total ligand atoms (non-water HETATM).
    pub total_het_atoms: usize,
    /// Number of chains that are predominantly protein.
    pub num_protein_chains: usize,
    /// Number of chains that are predominantly nucleic acid.
    pub num_nucleic_chains: usize,
    /// Total water molecules (unique residues).
    pub num_water_molecules: usize,
}

impl fmt::Display for MolecularInventory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Structure inventory:")?;
        for chain in &self.chains {
            writeln!(f, "  {}", chain)?;
        }
        if !self.ligands.is_empty() {
            let lig_strs: Vec<String> = self.ligands.iter().map(|l| l.to_string()).collect();
            writeln!(f, "  Ligands: {}", lig_strs.join(", "))?;
        }
        if self.num_water_molecules > 0 {
            writeln!(f, "  Water: {} molecules", self.num_water_molecules)?;
        }
        write!(f, "  Total: {} atoms", self.total_atoms)
    }
}

/// Residue key for deduplication: (chain_id, residue_seq, ins_code).
type ResKey = (String, i32, Option<char>);

impl PdbStructure {
    /// Compute a one-call breakdown of all molecular entities in the structure.
    ///
    /// Returns a [`MolecularInventory`] describing every chain, ligand, water
    /// group and ion present, with atom and residue counts.
    ///
    /// This method requires no feature flags and works on any `PdbStructure`.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let structure = pdbrust::parse_pdb_file("1HSG.pdb")?;
    /// let inv = structure.molecular_inventory();
    ///
    /// println!("Protein chains: {}", inv.num_protein_chains);
    /// for lig in &inv.ligands {
    ///     println!("  Ligand: {}", lig);
    /// }
    /// println!("{}", inv);
    /// ```
    pub fn molecular_inventory(&self) -> MolecularInventory {
        // Ordered chain ids (by first appearance)
        let mut chain_order: Vec<String> = Vec::new();

        // Per-chain atom counts
        let mut chain_atoms: BTreeMap<String, usize> = BTreeMap::new();
        // Per-chain unique residue sets by category
        let mut chain_protein_res: BTreeMap<String, BTreeSet<ResKey>> = BTreeMap::new();
        let mut chain_nucleic_res: BTreeMap<String, BTreeSet<ResKey>> = BTreeMap::new();
        let mut chain_water_res: BTreeMap<String, BTreeSet<ResKey>> = BTreeMap::new();
        let mut chain_het_res: BTreeMap<String, BTreeSet<ResKey>> = BTreeMap::new();

        // Global counters
        let mut total_protein_atoms = 0usize;
        let mut total_nucleic_atoms = 0usize;
        let mut total_water_atoms = 0usize;
        let mut total_het_atoms = 0usize;

        // Ligand tracking: (name, chain, resid) -> atom count
        let mut ligand_counts: BTreeMap<(String, String, i32), usize> = BTreeMap::new();
        // Preserve ligand insertion order
        let mut ligand_order: Vec<(String, String, i32)> = Vec::new();
        let mut ligand_seen: BTreeSet<(String, String, i32)> = BTreeSet::new();

        for atom in &self.atoms {
            let cid = &atom.chain_id;

            // Track chain order
            if !chain_order.contains(cid) {
                chain_order.push(cid.clone());
            }
            *chain_atoms.entry(cid.clone()).or_insert(0) += 1;

            let res_key: ResKey = (cid.clone(), atom.residue_seq, atom.ins_code);
            let resname = atom.residue_name.trim();

            if is_water(resname) {
                chain_water_res
                    .entry(cid.clone())
                    .or_default()
                    .insert(res_key);
                total_water_atoms += 1;
            } else if !atom.is_hetatm && is_standard_amino_acid(resname) {
                chain_protein_res
                    .entry(cid.clone())
                    .or_default()
                    .insert(res_key);
                total_protein_atoms += 1;
            } else if !atom.is_hetatm && is_standard_nucleotide(resname) {
                chain_nucleic_res
                    .entry(cid.clone())
                    .or_default()
                    .insert(res_key);
                total_nucleic_atoms += 1;
            } else {
                // HETATM non-water: ligand, ion, modified residue, etc.
                chain_het_res
                    .entry(cid.clone())
                    .or_default()
                    .insert(res_key.clone());
                total_het_atoms += 1;

                // Track individual ligand instances (skip ions)
                if atom.is_hetatm && !COMMON_IONS.contains(&resname) {
                    let lig_key = (resname.to_string(), cid.clone(), atom.residue_seq);
                    *ligand_counts.entry(lig_key.clone()).or_insert(0) += 1;
                    if !ligand_seen.contains(&lig_key) {
                        ligand_seen.insert(lig_key.clone());
                        ligand_order.push(lig_key);
                    }
                }
            }
        }

        // Build per-chain inventory
        let mut chains = Vec::new();
        for cid in &chain_order {
            let num_atoms = chain_atoms.get(cid).copied().unwrap_or(0);
            let protein_residues = chain_protein_res.get(cid).map(|s| s.len()).unwrap_or(0);
            let nucleic_residues = chain_nucleic_res.get(cid).map(|s| s.len()).unwrap_or(0);
            let water_molecules = chain_water_res.get(cid).map(|s| s.len()).unwrap_or(0);
            let het_residues = chain_het_res.get(cid).map(|s| s.len()).unwrap_or(0);

            let chain_type = if protein_residues > 0 && nucleic_residues > 0 {
                ChainType::Mixed
            } else if protein_residues > 0 {
                ChainType::Protein
            } else if nucleic_residues > 0 {
                ChainType::NucleicAcid
            } else if water_molecules > 0 && het_residues == 0 {
                ChainType::Water
            } else {
                ChainType::Other
            };

            chains.push(ChainInventory {
                chain_id: cid.clone(),
                chain_type,
                num_atoms,
                protein_residues,
                nucleic_residues,
                water_molecules,
                het_residues,
            });
        }

        // Build ligand list
        let ligands: Vec<LigandInfo> = ligand_order
            .iter()
            .map(|(name, chain, resid)| LigandInfo {
                name: name.clone(),
                chain_id: chain.clone(),
                residue_seq: *resid,
                num_atoms: ligand_counts
                    .get(&(name.clone(), chain.clone(), *resid))
                    .copied()
                    .unwrap_or(0),
            })
            .collect();

        let num_protein_chains = chains
            .iter()
            .filter(|c| c.chain_type == ChainType::Protein || c.chain_type == ChainType::Mixed)
            .count();
        let num_nucleic_chains = chains
            .iter()
            .filter(|c| c.chain_type == ChainType::NucleicAcid || c.chain_type == ChainType::Mixed)
            .count();
        let num_water_molecules: usize = chains.iter().map(|c| c.water_molecules).sum();

        MolecularInventory {
            chains,
            ligands,
            total_atoms: self.atoms.len(),
            total_protein_atoms,
            total_nucleic_atoms,
            total_water_atoms,
            total_het_atoms,
            num_protein_chains,
            num_nucleic_chains,
            num_water_molecules,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn atom(
        chain: &str,
        resid: i32,
        resname: &str,
        name: &str,
        element: &str,
        is_hetatm: bool,
    ) -> Atom {
        Atom {
            serial: resid,
            name: name.to_string(),
            alt_loc: None,
            residue_name: resname.to_string(),
            chain_id: chain.to_string(),
            residue_seq: resid,
            ins_code: None,
            is_hetatm,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    #[test]
    fn test_empty_structure() {
        let s = PdbStructure::new();
        let inv = s.molecular_inventory();
        assert_eq!(inv.total_atoms, 0);
        assert!(inv.chains.is_empty());
        assert!(inv.ligands.is_empty());
    }

    #[test]
    fn test_protein_only() {
        let mut s = PdbStructure::new();
        s.atoms.push(atom("A", 1, "ALA", "N", "N", false));
        s.atoms.push(atom("A", 1, "ALA", "CA", "C", false));
        s.atoms.push(atom("A", 2, "GLY", "N", "N", false));
        s.atoms.push(atom("A", 2, "GLY", "CA", "C", false));

        let inv = s.molecular_inventory();
        assert_eq!(inv.total_atoms, 4);
        assert_eq!(inv.total_protein_atoms, 4);
        assert_eq!(inv.num_protein_chains, 1);
        assert_eq!(inv.chains.len(), 1);
        assert_eq!(inv.chains[0].chain_type, ChainType::Protein);
        assert_eq!(inv.chains[0].protein_residues, 2);
    }

    #[test]
    fn test_protein_with_ligand_and_water() {
        let mut s = PdbStructure::new();
        // Protein
        s.atoms.push(atom("A", 1, "ALA", "CA", "C", false));
        s.atoms.push(atom("A", 2, "GLY", "CA", "C", false));
        // Ligand
        s.atoms.push(atom("A", 100, "ATP", "C1", "C", true));
        s.atoms.push(atom("A", 100, "ATP", "N1", "N", true));
        // Water
        s.atoms.push(atom("A", 200, "HOH", "O", "O", true));
        s.atoms.push(atom("A", 201, "HOH", "O", "O", true));

        let inv = s.molecular_inventory();
        assert_eq!(inv.total_atoms, 6);
        assert_eq!(inv.total_protein_atoms, 2);
        assert_eq!(inv.total_water_atoms, 2);
        assert_eq!(inv.total_het_atoms, 2);
        assert_eq!(inv.num_water_molecules, 2);
        assert_eq!(inv.ligands.len(), 1);
        assert_eq!(inv.ligands[0].name, "ATP");
        assert_eq!(inv.ligands[0].num_atoms, 2);
        assert_eq!(inv.num_protein_chains, 1);
    }

    #[test]
    fn test_multi_chain() {
        let mut s = PdbStructure::new();
        s.atoms.push(atom("A", 1, "ALA", "CA", "C", false));
        s.atoms.push(atom("B", 1, "GLY", "CA", "C", false));
        s.atoms.push(atom("C", 1, "DA", "P", "P", false));

        let inv = s.molecular_inventory();
        assert_eq!(inv.chains.len(), 3);
        assert_eq!(inv.num_protein_chains, 2);
        assert_eq!(inv.num_nucleic_chains, 1);
        assert_eq!(inv.chains[0].chain_type, ChainType::Protein);
        assert_eq!(inv.chains[1].chain_type, ChainType::Protein);
        assert_eq!(inv.chains[2].chain_type, ChainType::NucleicAcid);
    }

    #[test]
    fn test_ions_not_listed_as_ligands() {
        let mut s = PdbStructure::new();
        s.atoms.push(atom("A", 1, "ALA", "CA", "C", false));
        s.atoms.push(atom("A", 500, "ZN", "ZN", "ZN", true));

        let inv = s.molecular_inventory();
        assert!(inv.ligands.is_empty());
        assert_eq!(inv.total_het_atoms, 1);
    }

    #[test]
    fn test_display() {
        let mut s = PdbStructure::new();
        s.atoms.push(atom("A", 1, "ALA", "CA", "C", false));
        s.atoms.push(atom("A", 100, "ATP", "C1", "C", true));
        s.atoms.push(atom("A", 200, "HOH", "O", "O", true));

        let inv = s.molecular_inventory();
        let text = inv.to_string();
        assert!(text.contains("Protein"));
        assert!(text.contains("ATP"));
        assert!(text.contains("Water"));
    }
}
