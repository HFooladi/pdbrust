//! Residue record structure and implementations
//!
//! This module defines the `Residue` struct, which represents a residue within a chain.
//! It includes fields for the residue name, sequence number, insertion code, and atoms.
//! The residue can be either a standard amino acid or a hetero atom.

use std::collections::HashMap;
use super::Atom;

/// Represents a residue within a chain.
#[derive(Debug, Clone)]
pub struct Residue {
    /// Residue name (e.g., "ALA", "LYS").
    pub name: String,
    /// Residue sequence number.
    pub number: i32,
    /// Insertion code (if any).
    pub ins_code: Option<char>,
    /// Atoms in this residue.
    pub atoms: Vec<Atom>,
    /// Whether this is a hetero residue (from HETATM records).
    pub is_hetero: bool,
}


impl Residue {
    /// Creates a new empty residue.
    pub fn new(name: String, number: i32, ins_code: Option<char>, is_hetero: bool) -> Self {
        Self {
            name,
            number,
            ins_code,
            atoms: Vec::new(),
            is_hetero,
        }
    }

    /// Returns a unique identifier for this residue.
    ///
    /// This combines the residue number and insertion code to create a unique key
    /// that can be used for lookups.
    pub fn id(&self) -> String {
        match self.ins_code {
            Some(code) => format!("{}{}", self.number, code),
            None => self.number.to_string(),
        }
    }

    /// Gets an atom by name, if present.
    pub fn get_atom_by_name(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|atom| atom.name == name)
    }

    /// Gets the CA atom, if present.
    pub fn get_ca_atom(&self) -> Option<&Atom> {
        self.get_atom_by_name("CA")
    }

    /// Gets the center of mass of the residue.
    pub fn center_of_mass(&self) -> (f64, f64, f64) {
        if self.atoms.is_empty() {
            return (0.0, 0.0, 0.0);
        }

        let mut x_sum = 0.0;
        let mut y_sum = 0.0;
        let mut z_sum = 0.0;

        for atom in &self.atoms {
            x_sum += atom.x;
            y_sum += atom.y;
            z_sum += atom.z;
        }

        let count = self.atoms.len() as f64;
        (x_sum / count, y_sum / count, z_sum / count)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    #[test]
    fn test_residue_creation() {
        let residue = Residue::new("ALA".to_string(), 1, None, false);
        
        assert_eq!(residue.name, "ALA");
        assert_eq!(residue.number, 1);
        assert_eq!(residue.ins_code, None);
        assert!(residue.atoms.is_empty());
        assert!(!residue.is_hetero);
    }

    #[test]
    fn test_residue_with_insertion_code() {
        let residue = Residue::new("ALA".to_string(), 1, Some('A'), false);
        
        assert_eq!(residue.name, "ALA");
        assert_eq!(residue.number, 1);
        assert_eq!(residue.ins_code, Some('A'));
        assert!(residue.atoms.is_empty());
        assert!(!residue.is_hetero);
    }

    #[test]
    fn test_residue_id() {
        let residue1 = Residue::new("ALA".to_string(), 1, None, false);
        let residue2 = Residue::new("ALA".to_string(), 1, Some('A'), false);
        
        assert_eq!(residue1.id(), "1");
        assert_eq!(residue2.id(), "1A");
    }

    #[test]
    fn test_residue_atom_management() {
        let mut residue = Residue::new("ALA".to_string(), 1, None, false);
        
        // Create test atoms
        let ca_atom = Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            1.0,
            2.0,
            3.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        
        let cb_atom = Atom::new(
            2,
            "CB".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            2.0,
            3.0,
            4.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        
        residue.atoms.push(ca_atom.clone());
        residue.atoms.push(cb_atom.clone());
        
        // Test atom retrieval
        assert_eq!(residue.get_atom_by_name("CA"), Some(&ca_atom));
        assert_eq!(residue.get_atom_by_name("CB"), Some(&cb_atom));
        assert_eq!(residue.get_atom_by_name("N"), None);
        
        // Test CA atom retrieval
        assert_eq!(residue.get_ca_atom(), Some(&ca_atom));
    }

    #[test]
    fn test_residue_center_of_mass() {
        let mut residue = Residue::new("ALA".to_string(), 1, None, false);
        
        // Add atoms at different positions
        residue.atoms.push(Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        ));
        
        residue.atoms.push(Atom::new(
            2,
            "CB".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            1.0,
            1.0,
            1.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        ));
        
        let (x, y, z) = residue.center_of_mass();
        assert!((x - 0.5).abs() < 1e-6);
        assert!((y - 0.5).abs() < 1e-6);
        assert!((z - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_empty_residue_center_of_mass() {
        let residue = Residue::new("ALA".to_string(), 1, None, false);
        let (x, y, z) = residue.center_of_mass();
        assert_eq!((x, y, z), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_hetero_residue() {
        let residue = Residue::new("HOH".to_string(), 1, None, true);
        
        assert_eq!(residue.name, "HOH");
        assert!(residue.is_hetero);
    }

    #[test]
    fn test_residue_with_multiple_atoms() {
        let mut residue = Residue::new("ALA".to_string(), 1, None, false);
        
        // Add multiple atoms
        for i in 0..3 {
            residue.atoms.push(Atom::new(
                i + 1,
                format!("ATOM{}", i + 1),
                None,
                "ALA".to_string(),
                "A".to_string(),
                1,
                i as f64,
                i as f64,
                i as f64,
                1.0,
                20.0,
                "C".to_string(),
                None,
            ));
        }
        
        assert_eq!(residue.atoms.len(), 3);
        
        // Test center of mass with multiple atoms
        let (x, y, z) = residue.center_of_mass();
        assert!((x - 1.0).abs() < 1e-6);
        assert!((y - 1.0).abs() < 1e-6);
        assert!((z - 1.0).abs() < 1e-6);
    }
}
