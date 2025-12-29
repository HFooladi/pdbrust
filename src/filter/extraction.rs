//! Extraction utilities for PDB structures.
//!
//! Functions for extracting specific subsets of atoms or creating
//! filtered versions of structures.

use crate::core::PdbStructure;
use super::is_standard_residue;

impl PdbStructure {
    /// Extract Cα (alpha-carbon) coordinates from the structure.
    ///
    /// Returns a vector of (x, y, z) coordinate tuples for all Cα atoms.
    /// Optionally filter by chain ID.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - Optional chain identifier to filter by. If `None`, returns
    ///   Cα coordinates from all chains.
    ///
    /// # Returns
    ///
    /// A vector of (x, y, z) tuples representing Cα atom positions in Angstroms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("example.pdb")?;
    ///
    /// // Get all Cα coordinates
    /// let all_ca = structure.get_ca_coords(None);
    ///
    /// // Get Cα coordinates for chain A only
    /// let chain_a_ca = structure.get_ca_coords(Some("A"));
    /// ```
    pub fn get_ca_coords(&self, chain_id: Option<&str>) -> Vec<(f64, f64, f64)> {
        self.atoms
            .iter()
            .filter(|atom| {
                atom.name.trim() == "CA"
                    && chain_id.is_none_or(|c| atom.chain_id == c)
            })
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect()
    }

    /// Get all Cα atoms from the structure.
    ///
    /// Similar to `get_ca_coords` but returns full Atom references instead of just coordinates.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - Optional chain identifier to filter by.
    ///
    /// # Returns
    ///
    /// A vector of references to Cα atoms.
    pub fn get_ca_atoms(&self, chain_id: Option<&str>) -> Vec<&crate::records::Atom> {
        self.atoms
            .iter()
            .filter(|atom| {
                atom.name.trim() == "CA"
                    && chain_id.is_none_or(|c| atom.chain_id == c)
            })
            .collect()
    }

    /// Remove ligands and heteroatoms, keeping only standard biological residues.
    ///
    /// Creates a new `PdbStructure` containing only atoms from standard amino acids
    /// and nucleotides. Water molecules (HOH), ions, and other ligands are removed.
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` containing only protein/nucleic acid atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("complex.pdb")?;
    /// let protein_only = structure.remove_ligands();
    ///
    /// // Save the cleaned structure
    /// protein_only.to_file("protein_only.pdb")?;
    /// ```
    pub fn remove_ligands(&self) -> Self {
        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| is_standard_residue(&atom.residue_name))
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: self.seqres.clone(),
            connects: Vec::new(), // Connectivity may be invalid after filtering
            ssbonds: self.ssbonds.clone(),
            remarks: self.remarks.clone(),
            models: Vec::new(), // Models need special handling
            current_model: self.current_model,
        }
    }

    /// Extract a single chain from the structure.
    ///
    /// Creates a new `PdbStructure` containing only atoms from the specified chain.
    ///
    /// # Arguments
    ///
    /// * `chain_id` - The chain identifier to extract (e.g., "A", "B").
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` containing only atoms from the specified chain.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("complex.pdb")?;
    /// let chain_a = structure.keep_only_chain("A");
    ///
    /// println!("Chain A has {} atoms", chain_a.get_num_atoms());
    /// ```
    pub fn keep_only_chain(&self, chain_id: &str) -> Self {
        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| atom.chain_id == chain_id)
            .cloned()
            .collect();

        // Filter SEQRES to only include the specified chain
        let filtered_seqres: Vec<_> = self.seqres
            .iter()
            .filter(|s| s.chain_id == chain_id)
            .cloned()
            .collect();

        // Filter SSBONDs that involve this chain
        let filtered_ssbonds: Vec<_> = self.ssbonds
            .iter()
            .filter(|s| s.chain1_id == chain_id || s.chain2_id == chain_id)
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: filtered_seqres,
            connects: Vec::new(), // Rebuild if needed
            ssbonds: filtered_ssbonds,
            remarks: self.remarks.clone(),
            models: Vec::new(),
            current_model: self.current_model,
        }
    }

    /// Keep only Cα atoms in the structure.
    ///
    /// Creates a new `PdbStructure` containing only alpha-carbon atoms.
    /// Useful for coarse-grained analysis and reducing structural complexity.
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` with only Cα atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let ca_only = structure.keep_only_ca();
    ///
    /// // Much smaller structure
    /// assert!(ca_only.get_num_atoms() < structure.get_num_atoms());
    /// ```
    pub fn keep_only_ca(&self) -> Self {
        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: self.seqres.clone(),
            connects: Vec::new(),
            ssbonds: self.ssbonds.clone(),
            remarks: self.remarks.clone(),
            models: Vec::new(),
            current_model: self.current_model,
        }
    }

    /// Keep only backbone atoms (N, CA, C, O).
    ///
    /// Creates a new `PdbStructure` containing only backbone atoms.
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` with only backbone atoms.
    pub fn keep_only_backbone(&self) -> Self {
        let backbone_names = ["N", "CA", "C", "O"];

        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| backbone_names.contains(&atom.name.trim()))
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: self.seqres.clone(),
            connects: Vec::new(),
            ssbonds: self.ssbonds.clone(),
            remarks: self.remarks.clone(),
            models: Vec::new(),
            current_model: self.current_model,
        }
    }

    /// Remove hydrogen atoms from the structure.
    ///
    /// Creates a new `PdbStructure` without hydrogen atoms.
    /// Useful for analysis that focuses on heavy atoms only.
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` without hydrogen atoms.
    pub fn remove_hydrogens(&self) -> Self {
        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| !atom.is_hydrogen())
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: self.seqres.clone(),
            connects: Vec::new(),
            ssbonds: self.ssbonds.clone(),
            remarks: self.remarks.clone(),
            models: Vec::new(),
            current_model: self.current_model,
        }
    }

    /// Filter atoms by a custom predicate.
    ///
    /// Creates a new `PdbStructure` containing only atoms that satisfy
    /// the given predicate function.
    ///
    /// # Arguments
    ///
    /// * `predicate` - A function that takes an `&Atom` and returns `true`
    ///   to keep the atom, `false` to remove it.
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` with filtered atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    ///
    /// // Keep only atoms with temperature factor < 50
    /// let low_bfactor = structure.filter_atoms(|atom| atom.temp_factor < 50.0);
    /// ```
    pub fn filter_atoms<F>(&self, predicate: F) -> Self
    where
        F: Fn(&crate::records::Atom) -> bool,
    {
        let filtered_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| predicate(atom))
            .cloned()
            .collect();

        Self {
            atoms: filtered_atoms,
            header: self.header.clone(),
            title: self.title.clone(),
            seqres: self.seqres.clone(),
            connects: Vec::new(),
            ssbonds: self.ssbonds.clone(),
            remarks: self.remarks.clone(),
            models: Vec::new(),
            current_model: self.current_model,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            // Chain A - ALA residue (standard amino acid)
            Atom {
                serial: 1,
                name: " N  ".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 0.0, y: 0.0, z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "N".to_string(),
            },
            Atom {
                serial: 2,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 1.0, y: 1.0, z: 1.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            Atom {
                serial: 3,
                name: " C  ".to_string(),
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 2.0, y: 2.0, z: 2.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "C".to_string(),
            },
            // Chain A - GLY residue
            Atom {
                serial: 4,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "GLY".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 2,
                ins_code: None,
                x: 3.0, y: 3.0, z: 3.0,
                occupancy: 1.0,
                temp_factor: 25.0,
                element: "C".to_string(),
            },
            // Chain B - VAL residue
            Atom {
                serial: 5,
                name: " CA ".to_string(),
                alt_loc: None,
                residue_name: "VAL".to_string(),
                chain_id: "B".to_string(),
                residue_seq: 1,
                ins_code: None,
                x: 10.0, y: 10.0, z: 10.0,
                occupancy: 1.0,
                temp_factor: 30.0,
                element: "C".to_string(),
            },
            // Water molecule (HETATM)
            Atom {
                serial: 6,
                name: " O  ".to_string(),
                alt_loc: None,
                residue_name: "HOH".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 100,
                ins_code: None,
                x: 5.0, y: 5.0, z: 5.0,
                occupancy: 1.0,
                temp_factor: 50.0,
                element: "O".to_string(),
            },
            // Ligand
            Atom {
                serial: 7,
                name: " C1 ".to_string(),
                alt_loc: None,
                residue_name: "ATP".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 101,
                ins_code: None,
                x: 6.0, y: 6.0, z: 6.0,
                occupancy: 1.0,
                temp_factor: 40.0,
                element: "C".to_string(),
            },
        ];

        structure
    }

    #[test]
    fn test_get_ca_coords_all() {
        let structure = create_test_structure();
        let coords = structure.get_ca_coords(None);

        assert_eq!(coords.len(), 3); // 3 CA atoms total
        assert_eq!(coords[0], (1.0, 1.0, 1.0));
        assert_eq!(coords[1], (3.0, 3.0, 3.0));
        assert_eq!(coords[2], (10.0, 10.0, 10.0));
    }

    #[test]
    fn test_get_ca_coords_chain_a() {
        let structure = create_test_structure();
        let coords = structure.get_ca_coords(Some("A"));

        assert_eq!(coords.len(), 2); // 2 CA atoms in chain A
        assert_eq!(coords[0], (1.0, 1.0, 1.0));
        assert_eq!(coords[1], (3.0, 3.0, 3.0));
    }

    #[test]
    fn test_get_ca_coords_chain_b() {
        let structure = create_test_structure();
        let coords = structure.get_ca_coords(Some("B"));

        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0], (10.0, 10.0, 10.0));
    }

    #[test]
    fn test_get_ca_coords_nonexistent_chain() {
        let structure = create_test_structure();
        let coords = structure.get_ca_coords(Some("Z"));

        assert!(coords.is_empty());
    }

    #[test]
    fn test_remove_ligands() {
        let structure = create_test_structure();
        let cleaned = structure.remove_ligands();

        // Should remove HOH and ATP
        assert_eq!(cleaned.get_num_atoms(), 5);

        // Check that all remaining atoms are standard amino acids
        for atom in &cleaned.atoms {
            assert!(
                super::super::is_standard_residue(&atom.residue_name),
                "Found non-standard residue: {}",
                atom.residue_name
            );
        }
    }

    #[test]
    fn test_keep_only_chain() {
        let structure = create_test_structure();
        let chain_a = structure.keep_only_chain("A");

        assert_eq!(chain_a.get_num_chains(), 1);

        for atom in &chain_a.atoms {
            assert_eq!(atom.chain_id, "A");
        }
    }

    #[test]
    fn test_keep_only_ca() {
        let structure = create_test_structure();
        let ca_only = structure.keep_only_ca();

        assert_eq!(ca_only.get_num_atoms(), 3);

        for atom in &ca_only.atoms {
            assert_eq!(atom.name.trim(), "CA");
        }
    }

    #[test]
    fn test_keep_only_backbone() {
        let structure = create_test_structure();
        let backbone = structure.keep_only_backbone();

        // We have N, CA, C from ALA (3) + CA from GLY (1) + CA from VAL (1) + O from HOH (1)
        // = 6 backbone atom names (note: water O matches backbone "O" name)
        // To get protein-only backbone, use remove_ligands() first
        assert_eq!(backbone.get_num_atoms(), 6);

        for atom in &backbone.atoms {
            let name = atom.name.trim();
            assert!(
                ["N", "CA", "C", "O"].contains(&name),
                "Found non-backbone atom: {}",
                name
            );
        }
    }

    #[test]
    fn test_keep_only_backbone_protein_only() {
        let structure = create_test_structure();
        // Chain remove_ligands + keep_only_backbone for protein-only backbone
        let backbone = structure.remove_ligands().keep_only_backbone();

        // After removing ligands: N, CA, C from ALA (3) + CA from GLY (1) + CA from VAL (1)
        // = 5 backbone atoms
        assert_eq!(backbone.get_num_atoms(), 5);
    }

    #[test]
    fn test_filter_atoms_custom() {
        let structure = create_test_structure();

        // Filter atoms with temp_factor < 30
        let low_bfactor = structure.filter_atoms(|atom| atom.temp_factor < 30.0);

        for atom in &low_bfactor.atoms {
            assert!(atom.temp_factor < 30.0);
        }
    }

    #[test]
    fn test_chained_operations() {
        let structure = create_test_structure();

        // Chain operations: remove ligands, then keep only CA
        let result = structure.remove_ligands().keep_only_ca();

        assert_eq!(result.get_num_atoms(), 3);

        for atom in &result.atoms {
            assert_eq!(atom.name.trim(), "CA");
            assert!(super::super::is_standard_residue(&atom.residue_name));
        }
    }
}
