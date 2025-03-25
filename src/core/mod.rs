use std::path::Path;
use crate::error::PdbError;
use crate::records::{Atom, Model, SeqRes, Conect, SSBond, Remark};

/// Main structure representing a PDB file.
#[derive(Debug, Clone)]
pub struct PdbStructure {
    /// List of all atoms in the structure.
    pub atoms: Vec<Atom>,
    /// Header information from the PDB file.
    pub header: Option<String>,
    /// Title of the structure.
    pub title: Option<String>,
    /// Sequence information from SEQRES records.
    pub seqres: Vec<SeqRes>,
    /// Connectivity information from CONECT records.
    pub connects: Vec<Conect>,
    /// Disulfide bond information from SSBOND records.
    pub ssbonds: Vec<SSBond>,
    /// Additional information from REMARK records.
    pub remarks: Vec<Remark>,
    /// List of models in case of multi-model structures.
    pub models: Vec<Model>,
    /// Currently active model number (if any).
    pub current_model: Option<i32>,
}

impl PdbStructure {
    /// Creates a new empty PDB structure.
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            header: None,
            title: None,
            seqres: Vec::new(),
            connects: Vec::new(),
            ssbonds: Vec::new(),
            remarks: Vec::new(),
            models: Vec::new(),
            current_model: None,
        }
    }

    /// Creates a PDB structure from a file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, PdbError> {
        crate::parser::parse_pdb_file(path)
    }

    /// Gets all unique chain IDs in the structure.
    pub fn get_chain_ids(&self) -> Vec<String> {
        let mut chain_ids = std::collections::HashSet::new();
        for atom in &self.atoms {
            chain_ids.insert(atom.chain_id.clone());
        }
        let mut chain_ids: Vec<String> = chain_ids.into_iter().collect();
        chain_ids.sort();
        chain_ids
    }

    /// Gets all residues for a specific chain.
    pub fn get_residues_for_chain(&self, chain_id: &str) -> Vec<(i32, String)> {
        let mut residues = std::collections::HashSet::new();
        for atom in &self.atoms {
            if atom.chain_id == chain_id {
                residues.insert((atom.residue_seq, atom.residue_name.clone()));
            }
        }
        let mut residues: Vec<(i32, String)> = residues.into_iter().collect();
        residues.sort_by_key(|&(num, _)| num);
        residues
    }

    /// Gets the sequence for a specific chain.
    pub fn get_sequence(&self, chain_id: &str) -> Vec<String> {
        let mut sequence = Vec::new();
        for seqres in &self.seqres {
            if seqres.chain_id == chain_id {
                sequence.extend(seqres.residues.clone());
            }
        }
        sequence
    }

    /// Gets all remarks with a specific number.
    pub fn get_remarks_by_number(&self, number: i32) -> Vec<&Remark> {
        self.remarks.iter().filter(|r| r.number == number).collect()
    }

    /// Gets all atoms connected to a specific atom.
    pub fn get_connected_atoms(&self, atom_serial: i32) -> Vec<&Atom> {
        let mut connected = Vec::new();
        for conect in &self.connects {
            if conect.atom1 == atom_serial {
                if let Some(atom) = self.atoms.iter().find(|a| a.serial == conect.atom2) {
                    connected.push(atom);
                }
            }
            if conect.atom2 == atom_serial {
                if let Some(atom) = self.atoms.iter().find(|a| a.serial == conect.atom1) {
                    connected.push(atom);
                }
            }
        }
        connected
    }

    /// Translates all atoms in the structure by the given coordinates.
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.atoms {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }

    /// Writes the structure to a file.
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), PdbError> {
        crate::writer::write_pdb_file(self, path)
    }
}