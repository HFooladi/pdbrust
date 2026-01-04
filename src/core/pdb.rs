//! Core PDB structure representation and manipulation.
//!
//! This module provides the main structure for representing and manipulating PDB (Protein Data Bank) files.
//! The `PdbStructure` struct is the central data structure that holds all information from a PDB file,
//! including atoms, models, connectivity information, and metadata.
//!
//! # Examples
//! ```ignore
//! use pdbrust::core::PdbStructure;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Create a new empty structure
//!     let structure = PdbStructure::new();
//!
//!     // Load a structure from a file
//!     let structure = PdbStructure::from_file("path/to/structure.pdb")?;
//!
//!     // Get all chain IDs
//!     let chain_ids = structure.get_chain_ids();
//!
//!     // Get residues for a specific chain
//!     let residues = structure.get_residues_for_chain("A");
//!     Ok(())
//! }
//! ```

use crate::error::PdbError;
use crate::records::{Atom, Conect, Model, Remark, SSBond, SeqRes};
use std::path::Path;

/// Main structure representing a PDB file.
///
/// This structure holds all the information from a PDB file, including:
/// - Atomic coordinates and properties
/// - Structure metadata (header, title)
/// - Sequence information
/// - Connectivity data
/// - Disulfide bonds
/// - Remarks and annotations
/// - Multiple models (for NMR structures)
///
/// # Fields
/// * `atoms` - List of all atoms in the structure
/// * `header` - Header information from the PDB file
/// * `title` - Title of the structure
/// * `seqres` - Sequence information from SEQRES records
/// * `connects` - Connectivity information from CONECT records
/// * `ssbonds` - Disulfide bond information from SSBOND records
/// * `remarks` - Additional information from REMARK records
/// * `models` - List of models in case of multi-model structures
/// * `current_model` - Currently active model number (if any)
///
/// # Examples
/// ```ignore
/// use pdbrust::core::PdbStructure;
///
/// // Create a new structure
/// let structure = PdbStructure::new();
///
/// // Load from file
/// let structure = PdbStructure::from_file("example.pdb")?;
/// ```
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
    ///
    /// Returns a `PdbStructure` with empty vectors for all collections and `None` for optional fields.
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::new();
    /// assert!(structure.atoms.is_empty());
    /// assert!(structure.header.is_none());
    /// ```
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
    ///
    /// # Arguments
    /// * `path` - Path to the PDB file to load
    ///
    /// # Returns
    /// * `Result<Self, PdbError>` - The parsed PDB structure or an error if parsing fails
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// ```
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, PdbError> {
        crate::parser::parse_pdb_file(path)
    }

    /// Get number of atoms in the structure.
    ///
    /// # Returns
    /// * `usize` - Number of atoms in the structure
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let num_atoms = structure.get_num_atoms();
    /// println!("Number of atoms: {}", num_atoms);
    /// ```
    pub fn get_num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get number of chains in the structure.
    ///
    /// # Returns
    /// * `usize` - Number of chains in the structure
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let num_chains = structure.get_num_chains();
    /// println!("Number of chains: {}", num_chains);
    /// ```
    pub fn get_num_chains(&self) -> usize {
        let mut chain_ids = std::collections::HashSet::new();
        for atom in &self.atoms {
            chain_ids.insert(atom.chain_id.clone());
        }
        chain_ids.len()
    }

    /// Gets all unique chain IDs in the structure.
    ///
    /// Returns a sorted vector of all unique chain IDs found in the structure's atoms.
    ///
    /// # Returns
    /// * `Vec<String>` - Sorted list of unique chain IDs
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let chain_ids = structure.get_chain_ids();
    /// println!("Found chains: {:?}", chain_ids);
    /// ```
    pub fn get_chain_ids(&self) -> Vec<String> {
        let mut chain_ids = std::collections::HashSet::new();
        for atom in &self.atoms {
            chain_ids.insert(atom.chain_id.clone());
        }
        let mut chain_ids: Vec<String> = chain_ids.into_iter().collect();
        chain_ids.sort();
        chain_ids
    }

    /// Get number of residues in the structure.           
    ///
    /// # Returns
    /// * `usize` - Number of residues in the structure
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;        
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let num_residues = structure.get_num_residues();
    /// println!("Number of residues: {}", num_residues);
    /// ```
    pub fn get_num_residues(&self) -> usize {
        let mut residues = std::collections::HashSet::new();
        for atom in &self.atoms {
            residues.insert((atom.residue_seq, atom.residue_name.clone()));
        }
        residues.len()
    }

    /// Get all residues in the structure.
    ///
    /// Returns a sorted vector of tuples containing residue sequence numbers and residue names.
    ///
    /// # Returns
    /// * `Vec<(i32, String)>` - Sorted list of (sequence number, residue name) tuples
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let residues = structure.get_residues();    
    /// for (seq_num, res_name) in residues {
    ///     println!("Residue {}: {}", seq_num, res_name);
    /// }
    /// ```
    pub fn get_residues(&self) -> Vec<(i32, String)> {
        let mut residues = std::collections::HashSet::new();
        for atom in &self.atoms {
            residues.insert((atom.residue_seq, atom.residue_name.clone()));
        }
        let mut residues: Vec<(i32, String)> = residues.into_iter().collect();
        residues.sort_by_key(|&(num, _)| num);
        residues
    }

    /// Gets all residues for a specific chain.
    ///
    /// Returns a sorted vector of tuples containing residue sequence numbers and residue names
    /// for the specified chain.
    ///
    /// # Arguments
    /// * `chain_id` - The chain ID to get residues for
    ///
    /// # Returns
    /// * `Vec<(i32, String)>` - Sorted list of (sequence number, residue name) tuples
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let residues = structure.get_residues_for_chain("A");
    /// for (seq_num, res_name) in residues {
    ///     println!("Residue {}: {}", seq_num, res_name);
    /// }
    /// ```
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
    ///
    /// Returns a vector of residue names representing the sequence of the specified chain,
    /// based on SEQRES records.
    ///
    /// # Arguments
    /// * `chain_id` - The chain ID to get the sequence for
    ///
    /// # Returns
    /// * `Vec<String>` - List of residue names in sequence order
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let sequence = structure.get_sequence("A");
    /// println!("Chain A sequence: {:?}", sequence);
    /// ```
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
    ///
    /// Returns a vector of references to remarks that match the specified remark number.
    ///
    /// # Arguments
    /// * `number` - The remark number to filter by
    ///
    /// # Returns
    /// * `Vec<&Remark>` - List of remarks with the specified number
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let resolution_remarks = structure.get_remarks_by_number(2);
    /// for remark in resolution_remarks {
    ///     println!("Resolution remark: {}", remark.text);
    /// }
    /// ```
    pub fn get_remarks_by_number(&self, number: i32) -> Vec<&Remark> {
        self.remarks.iter().filter(|r| r.number == number).collect()
    }

    /// Gets all atoms connected to a specific atom.
    ///
    /// Returns a vector of references to atoms that are connected to the specified atom
    /// through CONECT records.
    ///
    /// # Arguments
    /// * `atom_serial` - The serial number of the atom to find connections for
    ///
    /// # Returns
    /// * `Vec<&Atom>` - List of connected atoms
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let connected = structure.get_connected_atoms(1);
    /// println!("Atoms connected to atom 1: {:?}", connected);
    /// ```
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
    ///
    /// Moves all atoms in the structure by the specified x, y, and z offsets.
    ///
    /// # Arguments
    /// * `dx` - X-axis translation
    /// * `dy` - Y-axis translation
    /// * `dz` - Z-axis translation
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let mut structure = PdbStructure::from_file("example.pdb")?;
    /// // Move structure 1.0 Å in each direction
    /// structure.translate(1.0, 1.0, 1.0);
    /// ```
    pub fn translate(&mut self, dx: f64, dy: f64, dz: f64) {
        for atom in &mut self.atoms {
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
        }
    }

    /// Writes the structure to a file.
    ///
    /// Saves the current structure to a PDB file at the specified path.
    ///
    /// # Arguments
    /// * `path` - Path where the PDB file should be written
    ///
    /// # Returns
    /// * `Result<(), PdbError>` - Success or error if writing fails
    ///
    /// # Examples
    /// ```ignore
    /// use pdbrust::core::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("input.pdb")?;
    /// structure.to_file("output.pdb")?;
    /// ```
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), PdbError> {
        crate::writer::write_pdb_file(self, path)
    }
}

/// Implements the Default trait for PdbStructure.
///
/// Creates a new empty PDB structure when using `PdbStructure::default()`.
impl Default for PdbStructure {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::{Atom, Conect, Remark, SSBond, SeqRes};

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();

        // Add some test atoms
        structure.atoms = vec![
            Atom {
                serial: 1,
                name: "CA".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 0.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            Atom {
                serial: 2,
                name: "CA".to_string(),
                alt_loc: None,
                residue_name: "GLY".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 2,
                ins_code: None,
                x: 1.0,
                y: 1.0,
                z: 1.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
        ];

        // Add test connectivity
        structure.connects = vec![Conect {
            atom1: 1,
            atom2: 2,
            atom3: None,
            atom4: None,
        }];

        // Add test sequence
        structure.seqres = vec![SeqRes {
            serial: 1,
            chain_id: "A".to_string(),
            num_residues: 2,
            residues: vec!["ALA".to_string(), "GLY".to_string()],
        }];

        // Add test remarks
        structure.remarks = vec![Remark {
            number: 2,
            content: "RESOLUTION. 2.0 ANGSTROMS.".to_string(),
        }];

        // Add test SSBOND
        structure.ssbonds = vec![SSBond {
            serial: 1,
            residue1_name: "CYS".to_string(),
            chain1_id: "A".to_string(),
            residue1_seq: 1,
            icode1: None,
            residue2_name: "CYS".to_string(),
            chain2_id: "A".to_string(),
            residue2_seq: 2,
            icode2: None,
            sym1: 1,
            sym2: 1,
            length: 2.0,
        }];

        structure
    }

    #[test]
    fn test_new_structure() {
        let structure = PdbStructure::new();
        assert!(structure.atoms.is_empty());
        assert!(structure.header.is_none());
        assert!(structure.title.is_none());
        assert!(structure.seqres.is_empty());
        assert!(structure.connects.is_empty());
        assert!(structure.ssbonds.is_empty());
        assert!(structure.remarks.is_empty());
        assert!(structure.models.is_empty());
        assert!(structure.current_model.is_none());
    }

    #[test]
    fn test_get_chain_ids() {
        let structure = create_test_structure();
        let chain_ids = structure.get_chain_ids();
        assert_eq!(chain_ids, vec!["A"]);
    }

    #[test]
    fn test_get_residues_for_chain() {
        let structure = create_test_structure();
        let residues = structure.get_residues_for_chain("A");
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0], (1, "ALA".to_string()));
        assert_eq!(residues[1], (2, "GLY".to_string()));
    }

    #[test]
    fn test_get_sequence() {
        let structure = create_test_structure();
        let sequence = structure.get_sequence("A");
        assert_eq!(sequence, vec!["ALA", "GLY"]);
    }

    #[test]
    fn test_get_remarks_by_number() {
        let structure = create_test_structure();
        let remarks = structure.get_remarks_by_number(2);
        assert_eq!(remarks.len(), 1);
        assert!(remarks[0].content.contains("RESOLUTION"));
    }

    #[test]
    fn test_get_connected_atoms() {
        let structure = create_test_structure();
        let connected = structure.get_connected_atoms(1);
        assert_eq!(connected.len(), 1);
        assert_eq!(connected[0].serial, 2);
    }

    #[test]
    fn test_translate() {
        let mut structure = create_test_structure();
        let original_x = structure.atoms[0].x;
        let original_y = structure.atoms[0].y;
        let original_z = structure.atoms[0].z;

        structure.translate(1.0, 2.0, 3.0);

        assert_eq!(structure.atoms[0].x, original_x + 1.0);
        assert_eq!(structure.atoms[0].y, original_y + 2.0);
        assert_eq!(structure.atoms[0].z, original_z + 3.0);
    }

    #[test]
    fn test_default() {
        let structure = PdbStructure::default();
        assert!(structure.atoms.is_empty());
        assert!(structure.header.is_none());
        assert!(structure.title.is_none());
        assert!(structure.seqres.is_empty());
        assert!(structure.connects.is_empty());
        assert!(structure.ssbonds.is_empty());
        assert!(structure.remarks.is_empty());
        assert!(structure.models.is_empty());
        assert!(structure.current_model.is_none());
    }

    #[test]
    fn test_clone() {
        let structure = create_test_structure();
        let cloned = structure.clone();

        assert_eq!(structure.atoms.len(), cloned.atoms.len());
        assert_eq!(structure.connects.len(), cloned.connects.len());
        assert_eq!(structure.seqres.len(), cloned.seqres.len());
        assert_eq!(structure.remarks.len(), cloned.remarks.len());
        assert_eq!(structure.ssbonds.len(), cloned.ssbonds.len());

        // Test deep clone by modifying the original
        let mut structure = structure;
        structure.atoms[0].x = 100.0;
        assert_ne!(structure.atoms[0].x, cloned.atoms[0].x);
    }

    #[test]
    fn test_debug() {
        let structure = create_test_structure();
        let debug_string = format!("{:?}", structure);
        assert!(debug_string.contains("atoms"));
        assert!(debug_string.contains("connects"));
        assert!(debug_string.contains("seqres"));
        assert!(debug_string.contains("remarks"));
        assert!(debug_string.contains("ssbonds"));
    }
}
