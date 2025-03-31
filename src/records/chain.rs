use super::atom::Atom;
use super::residue::Residue;
/// Represents a single chain within a protein structure.
///
/// A Chain contains multiple residues identified by a unique chain identifier.
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Chain {
    /// Chain identifier (usually a single character).
    pub id: String,
    /// Residues in this chain, mapped by residue ID.
    residues: HashMap<String, Residue>,
    /// List of residue IDs in sequential order.
    residue_order: Vec<String>,
}

impl Chain {
    /// Creates a new empty chain with the given identifier.
    pub fn new(id: String) -> Self {
        Self {
            id,
            residues: HashMap::new(),
            residue_order: Vec::new(),
        }
    }

    /// Adds an atom to the chain, creating a new residue if needed.
    ///
    /// This is the primary method for building a chain from individual atom records.
    pub fn add_atom(&mut self, atom: Atom, is_hetero: bool) {
        let residue_id = match atom.ins_code {
            Some(code) => format!("{}{}", atom.residue_seq, code),
            None => atom.residue_seq.to_string(),
        };

        if !self.residues.contains_key(&residue_id) {
            let residue = Residue::new(
                atom.residue_name.clone(),
                atom.residue_seq,
                atom.ins_code,
                is_hetero,
            );
            self.residues.insert(residue_id.clone(), residue);
            self.residue_order.push(residue_id.clone());
        }

        if let Some(residue) = self.residues.get_mut(&residue_id) {
            residue.atoms.push(atom);
        }
    }

    /// Gets the number of residues in the chain.
    pub fn get_residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Gets the number of atoms in the chain.
    pub fn get_atom_count(&self) -> usize {
        self.residues.values().map(|r| r.atoms.len()).sum()
    }

    /// Gets a residue by its ID, if present.
    pub fn get_residue(&self, residue_id: &str) -> Option<&Residue> {
        self.residues.get(residue_id)
    }

    /// Gets a residue by its number and optional insertion code, if present.
    pub fn get_residue_by_number(&self, number: i32, ins_code: Option<char>) -> Option<&Residue> {
        let id = match ins_code {
            Some(code) => format!("{}{}", number, code),
            None => number.to_string(),
        };
        self.residues.get(&id)
    }

    /// Gets all residues in the chain in sequential order.
    pub fn get_residues(&self) -> Vec<&Residue> {
        self.residue_order
            .iter()
            .filter_map(|id| self.residues.get(id))
            .collect()
    }

    /// Gets all standard (non-hetero) residues in the chain.
    pub fn get_standard_residues(&self) -> Vec<&Residue> {
        self.get_residues()
            .into_iter()
            .filter(|r| !r.is_hetero)
            .collect()
    }

    /// Gets all hetero residues in the chain.
    pub fn get_hetero_residues(&self) -> Vec<&Residue> {
        self.get_residues()
            .into_iter()
            .filter(|r| r.is_hetero)
            .collect()
    }

    /// Gets the sequence of the chain as a vector of residue names.
    pub fn get_sequence(&self) -> Vec<String> {
        self.get_standard_residues()
            .iter()
            .map(|r| r.name.clone())
            .collect()
    }

    /// Gets the center of mass of the chain.
    pub fn get_center_of_mass(&self) -> (f64, f64, f64) {
        let mut x_sum = 0.0;
        let mut y_sum = 0.0;
        let mut z_sum = 0.0;
        let mut total_atoms = 0;

        for residue in self.residues.values() {
            for atom in &residue.atoms {
                x_sum += atom.x;
                y_sum += atom.y;
                z_sum += atom.z;
                total_atoms += 1;
            }
        }

        if total_atoms == 0 {
            return (0.0, 0.0, 0.0);
        }

        let count = total_atoms as f64;
        (x_sum / count, y_sum / count, z_sum / count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::atom::Atom;

    fn create_test_atom(residue_seq: i32, residue_name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: "A".to_string(),
            residue_seq,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        }
    }

    #[test]
    fn test_chain_creation() {
        let chain = Chain::new("A".to_string());
        assert_eq!(chain.id, "A");
        assert_eq!(chain.get_residue_count(), 0);
        assert_eq!(chain.get_atom_count(), 0);
    }

    #[test]
    fn test_add_atom() {
        let mut chain = Chain::new("A".to_string());
        let atom = create_test_atom(1, "ALA", 1.0, 2.0, 3.0);

        chain.add_atom(atom.clone(), false);
        assert_eq!(chain.get_residue_count(), 1);
        assert_eq!(chain.get_atom_count(), 1);

        let residue = chain.get_residue("1").unwrap();
        assert_eq!(residue.name, "ALA");
        assert_eq!(residue.atoms.len(), 1);
    }

    #[test]
    fn test_add_multiple_atoms() {
        let mut chain = Chain::new("A".to_string());

        // Add atoms for two residues
        chain.add_atom(create_test_atom(1, "ALA", 1.0, 2.0, 3.0), false);
        chain.add_atom(create_test_atom(2, "GLY", 4.0, 5.0, 6.0), false);

        assert_eq!(chain.get_residue_count(), 2);
        assert_eq!(chain.get_atom_count(), 2);

        let residues = chain.get_residues();
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].name, "ALA");
        assert_eq!(residues[1].name, "GLY");
    }

    #[test]
    fn test_hetero_atoms() {
        let mut chain = Chain::new("A".to_string());

        // Add standard and hetero atoms
        chain.add_atom(create_test_atom(1, "ALA", 1.0, 2.0, 3.0), false);
        chain.add_atom(create_test_atom(2, "HOH", 4.0, 5.0, 6.0), true);

        assert_eq!(chain.get_standard_residues().len(), 1);
        assert_eq!(chain.get_hetero_residues().len(), 1);
    }

    #[test]
    fn test_residue_lookup() {
        let mut chain = Chain::new("A".to_string());
        chain.add_atom(create_test_atom(1, "ALA", 1.0, 2.0, 3.0), false);

        // Test by ID
        assert!(chain.get_residue("1").is_some());
        assert!(chain.get_residue("2").is_none());

        // Test by number
        assert!(chain.get_residue_by_number(1, None).is_some());
        assert!(chain.get_residue_by_number(2, None).is_none());
    }

    #[test]
    fn test_sequence() {
        let mut chain = Chain::new("A".to_string());

        // Add some residues
        chain.add_atom(create_test_atom(1, "ALA", 1.0, 2.0, 3.0), false);
        chain.add_atom(create_test_atom(2, "GLY", 4.0, 5.0, 6.0), false);
        chain.add_atom(create_test_atom(3, "HOH", 7.0, 8.0, 9.0), true);

        let sequence = chain.get_sequence();
        assert_eq!(sequence, vec!["ALA", "GLY"]);
    }

    #[test]
    fn test_center_of_mass() {
        let mut chain = Chain::new("A".to_string());

        // Add atoms at known positions
        chain.add_atom(create_test_atom(1, "ALA", 0.0, 0.0, 0.0), false);
        chain.add_atom(create_test_atom(2, "GLY", 2.0, 2.0, 2.0), false);

        let (x, y, z) = chain.get_center_of_mass();
        assert!((x - 1.0).abs() < 1e-10);
        assert!((y - 1.0).abs() < 1e-10);
        assert!((z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_empty_chain_center_of_mass() {
        let chain = Chain::new("A".to_string());
        let (x, y, z) = chain.get_center_of_mass();
        assert_eq!((x, y, z), (0.0, 0.0, 0.0));
    }
}
