//! Structure cleaning and normalization utilities.
//!
//! Functions for modifying PDB structures in-place to normalize
//! chain identifiers, residue numbering, and coordinate systems.

use crate::core::PdbStructure;
use std::collections::HashMap;

impl PdbStructure {
    /// Normalize chain identifiers to consecutive letters (A, B, C, ...).
    ///
    /// Remaps all chain IDs to start from 'A' and continue alphabetically.
    /// Original chain ordering is preserved. This is useful for structures
    /// with non-standard or discontinuous chain identifiers.
    ///
    /// # Returns
    ///
    /// A mutable reference to self for method chaining.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let mut structure = PdbStructure::from_file("example.pdb")?;
    ///
    /// // Chains might be "X", "Y", "Z" -> become "A", "B", "C"
    /// structure.normalize_chain_ids();
    ///
    /// let chains = structure.get_chain_ids();
    /// assert_eq!(chains[0], "A");
    /// ```
    ///
    /// # Panics
    ///
    /// May produce unexpected results if there are more than 26 chains.
    pub fn normalize_chain_ids(&mut self) -> &mut Self {
        // Collect unique chain IDs in order of first appearance
        let mut chain_order: Vec<String> = Vec::new();
        for atom in &self.atoms {
            if !chain_order.contains(&atom.chain_id) {
                chain_order.push(atom.chain_id.clone());
            }
        }

        // Create mapping from old chain ID to new chain ID
        let mut chain_map: HashMap<String, String> = HashMap::new();
        for (i, old_chain) in chain_order.iter().enumerate() {
            let new_chain = if i < 26 {
                ((b'A' + i as u8) as char).to_string()
            } else {
                // For > 26 chains, use AA, AB, etc.
                let first = (b'A' + (i / 26 - 1) as u8) as char;
                let second = (b'A' + (i % 26) as u8) as char;
                format!("{}{}", first, second)
            };
            chain_map.insert(old_chain.clone(), new_chain);
        }

        // Apply mapping to atoms
        for atom in &mut self.atoms {
            if let Some(new_chain) = chain_map.get(&atom.chain_id) {
                atom.chain_id = new_chain.clone();
            }
        }

        // Apply mapping to SEQRES records
        for seqres in &mut self.seqres {
            if let Some(new_chain) = chain_map.get(&seqres.chain_id) {
                seqres.chain_id = new_chain.clone();
            }
        }

        // Apply mapping to SSBOND records
        for ssbond in &mut self.ssbonds {
            if let Some(new_chain) = chain_map.get(&ssbond.chain1_id) {
                ssbond.chain1_id = new_chain.clone();
            }
            if let Some(new_chain) = chain_map.get(&ssbond.chain2_id) {
                ssbond.chain2_id = new_chain.clone();
            }
        }

        // Apply mapping to models
        for model in &mut self.models {
            for atom in &mut model.atoms {
                if let Some(new_chain) = chain_map.get(&atom.chain_id) {
                    atom.chain_id = new_chain.clone();
                }
            }
        }

        self
    }

    /// Reindex residue numbers to be continuous starting from 1.
    ///
    /// Assigns new residue sequence numbers starting from 1 for each chain,
    /// eliminating gaps and non-standard numbering. Insertion codes are
    /// also cleared.
    ///
    /// # Returns
    ///
    /// A mutable reference to self for method chaining.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let mut structure = PdbStructure::from_file("example.pdb")?;
    ///
    /// // Residues might be 5, 10, 15 -> become 1, 2, 3
    /// structure.reindex_residues();
    /// ```
    pub fn reindex_residues(&mut self) -> &mut Self {
        // Group atoms by chain
        let chains = self.get_chain_ids();

        for chain_id in chains {
            // Collect unique residue identifiers for this chain in order
            let mut residue_order: Vec<(i32, Option<char>)> = Vec::new();
            for atom in &self.atoms {
                if atom.chain_id == chain_id {
                    let key = (atom.residue_seq, atom.ins_code);
                    if !residue_order.contains(&key) {
                        residue_order.push(key);
                    }
                }
            }

            // Create mapping from (old_seq, ins_code) to new_seq
            let mut residue_map: HashMap<(i32, Option<char>), i32> = HashMap::new();
            for (new_seq, old_key) in residue_order.iter().enumerate() {
                residue_map.insert(*old_key, (new_seq + 1) as i32);
            }

            // Apply mapping to atoms
            for atom in &mut self.atoms {
                if atom.chain_id == chain_id {
                    let key = (atom.residue_seq, atom.ins_code);
                    if let Some(&new_seq) = residue_map.get(&key) {
                        atom.residue_seq = new_seq;
                        atom.ins_code = None; // Clear insertion code
                    }
                }
            }
        }

        self
    }

    /// Center the structure at the origin (0, 0, 0).
    ///
    /// Calculates the centroid of all Cα atoms and translates all atoms
    /// so that the centroid is at the origin. If no Cα atoms are present,
    /// uses the centroid of all atoms instead.
    ///
    /// # Returns
    ///
    /// A mutable reference to self for method chaining.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let mut structure = PdbStructure::from_file("example.pdb")?;
    /// structure.center_structure();
    ///
    /// // The centroid should now be at approximately (0, 0, 0)
    /// ```
    pub fn center_structure(&mut self) -> &mut Self {
        // Calculate centroid from CA atoms
        let ca_coords = self.get_ca_coords(None);

        let (cx, cy, cz) = if !ca_coords.is_empty() {
            // Use CA centroid
            let n = ca_coords.len() as f64;
            (
                ca_coords.iter().map(|c| c.0).sum::<f64>() / n,
                ca_coords.iter().map(|c| c.1).sum::<f64>() / n,
                ca_coords.iter().map(|c| c.2).sum::<f64>() / n,
            )
        } else if !self.atoms.is_empty() {
            // Fallback to all atoms
            let n = self.atoms.len() as f64;
            (
                self.atoms.iter().map(|a| a.x).sum::<f64>() / n,
                self.atoms.iter().map(|a| a.y).sum::<f64>() / n,
                self.atoms.iter().map(|a| a.z).sum::<f64>() / n,
            )
        } else {
            // Empty structure, nothing to do
            return self;
        };

        // Translate to origin
        self.translate(-cx, -cy, -cz);

        self
    }

    /// Get the centroid (center of mass) of the structure.
    ///
    /// Calculates the geometric center of all atoms in the structure.
    ///
    /// # Returns
    ///
    /// A tuple (x, y, z) representing the centroid coordinates,
    /// or (0, 0, 0) for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    /// let (cx, cy, cz) = structure.get_centroid();
    /// println!("Centroid: ({}, {}, {})", cx, cy, cz);
    /// ```
    pub fn get_centroid(&self) -> (f64, f64, f64) {
        if self.atoms.is_empty() {
            return (0.0, 0.0, 0.0);
        }

        let n = self.atoms.len() as f64;
        (
            self.atoms.iter().map(|a| a.x).sum::<f64>() / n,
            self.atoms.iter().map(|a| a.y).sum::<f64>() / n,
            self.atoms.iter().map(|a| a.z).sum::<f64>() / n,
        )
    }

    /// Get the centroid of Cα atoms only.
    ///
    /// # Returns
    ///
    /// A tuple (x, y, z) representing the Cα centroid,
    /// or (0, 0, 0) if no Cα atoms present.
    pub fn get_ca_centroid(&self) -> (f64, f64, f64) {
        let ca_coords = self.get_ca_coords(None);

        if ca_coords.is_empty() {
            return (0.0, 0.0, 0.0);
        }

        let n = ca_coords.len() as f64;
        (
            ca_coords.iter().map(|c| c.0).sum::<f64>() / n,
            ca_coords.iter().map(|c| c.1).sum::<f64>() / n,
            ca_coords.iter().map(|c| c.2).sum::<f64>() / n,
        )
    }

    /// Renumber all atom serial numbers to be continuous starting from 1.
    ///
    /// # Returns
    ///
    /// A mutable reference to self for method chaining.
    pub fn renumber_atoms(&mut self) -> &mut Self {
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            atom.serial = (i + 1) as i32;
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            // Chain X - residue 10
            Atom {
                serial: 100,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "X".to_string(),
                residue_seq: 10,
                ins_code: None,
                x: 0.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            Atom {
                serial: 101,
                name: " N  ".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "X".to_string(),
                residue_seq: 10,
                ins_code: None,
                x: 1.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "N".to_string(),
            },
            // Chain X - residue 20 (gap in numbering)
            Atom {
                serial: 102,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "GLY".to_string(),
                chain_id: "X".to_string(),
                residue_seq: 20,
                ins_code: None,
                x: 4.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            // Chain Y - residue 5
            Atom {
                serial: 200,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "VAL".to_string(),
                chain_id: "Y".to_string(),
                residue_seq: 5,
                ins_code: None,
                x: 10.0,
                y: 10.0,
                z: 10.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
        ];

        structure
    }

    #[test]
    fn test_normalize_chain_ids() {
        let mut structure = create_test_structure();

        // Original chains are X and Y
        let original_chains = structure.get_chain_ids();
        assert!(original_chains.contains(&"X".to_string()));
        assert!(original_chains.contains(&"Y".to_string()));

        structure.normalize_chain_ids();

        // After normalization, should be A and B
        let new_chains = structure.get_chain_ids();
        assert!(new_chains.contains(&"A".to_string()));
        assert!(new_chains.contains(&"B".to_string()));
        assert!(!new_chains.contains(&"X".to_string()));
        assert!(!new_chains.contains(&"Y".to_string()));
    }

    #[test]
    fn test_reindex_residues() {
        let mut structure = create_test_structure();

        // Original residue numbers are 10, 20 (chain X) and 5 (chain Y)
        structure.reindex_residues();

        // After reindexing, chain X should have residues 1, 2 and chain Y should have residue 1
        let chain_x_atoms: Vec<_> = structure
            .atoms
            .iter()
            .filter(|a| a.chain_id == "X")
            .collect();

        assert_eq!(chain_x_atoms[0].residue_seq, 1);
        assert_eq!(chain_x_atoms[1].residue_seq, 1); // Same residue as first
        assert_eq!(chain_x_atoms[2].residue_seq, 2);

        let chain_y_atoms: Vec<_> = structure
            .atoms
            .iter()
            .filter(|a| a.chain_id == "Y")
            .collect();
        assert_eq!(chain_y_atoms[0].residue_seq, 1);
    }

    #[test]
    fn test_center_structure() {
        let mut structure = create_test_structure();

        structure.center_structure();

        // After centering, CA centroid should be at approximately (0, 0, 0)
        let (cx, cy, cz) = structure.get_ca_centroid();

        assert!(cx.abs() < 1e-10, "cx = {} should be ~0", cx);
        assert!(cy.abs() < 1e-10, "cy = {} should be ~0", cy);
        assert!(cz.abs() < 1e-10, "cz = {} should be ~0", cz);
    }

    #[test]
    fn test_get_centroid() {
        let structure = create_test_structure();
        let (cx, cy, cz) = structure.get_centroid();

        // Manual calculation: (0+1+4+10)/4=3.75, (0+0+0+10)/4=2.5, (0+0+0+10)/4=2.5
        assert!((cx - 3.75).abs() < 1e-10);
        assert!((cy - 2.5).abs() < 1e-10);
        assert!((cz - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_get_ca_centroid() {
        let structure = create_test_structure();
        let (cx, cy, cz) = structure.get_ca_centroid();

        // Only CA atoms: (0, 0, 0), (4, 0, 0), (10, 10, 10)
        // Average: (14/3, 10/3, 10/3)
        let expected_x = 14.0 / 3.0;
        let expected_y = 10.0 / 3.0;
        let expected_z = 10.0 / 3.0;

        assert!((cx - expected_x).abs() < 1e-10);
        assert!((cy - expected_y).abs() < 1e-10);
        assert!((cz - expected_z).abs() < 1e-10);
    }

    #[test]
    fn test_renumber_atoms() {
        let mut structure = create_test_structure();

        // Original serial numbers are 100, 101, 102, 200
        structure.renumber_atoms();

        // After renumbering, should be 1, 2, 3, 4
        assert_eq!(structure.atoms[0].serial, 1);
        assert_eq!(structure.atoms[1].serial, 2);
        assert_eq!(structure.atoms[2].serial, 3);
        assert_eq!(structure.atoms[3].serial, 4);
    }

    #[test]
    fn test_chained_cleaning_operations() {
        let mut structure = create_test_structure();

        structure
            .normalize_chain_ids()
            .reindex_residues()
            .renumber_atoms()
            .center_structure();

        // Verify all operations were applied
        let chains = structure.get_chain_ids();
        assert!(chains.contains(&"A".to_string()));
        assert!(chains.contains(&"B".to_string()));

        assert_eq!(structure.atoms[0].serial, 1);

        let (cx, cy, cz) = structure.get_ca_centroid();
        assert!(cx.abs() < 1e-10);
        assert!(cy.abs() < 1e-10);
        assert!(cz.abs() < 1e-10);
    }

    #[test]
    fn test_center_empty_structure() {
        let mut structure = PdbStructure::new();
        structure.center_structure();
        // Should not panic on empty structure
        assert!(structure.atoms.is_empty());
    }
}
