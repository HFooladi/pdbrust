//! Interface contact detection and fnat/fnonnat/F1 calculation.
//!
//! Detects residue-residue contacts at protein-protein interfaces
//! and computes the fraction of native contacts preserved in a model.

use std::collections::HashSet;

use crate::core::PdbStructure;
use crate::records::Atom;

/// A residue key for contact identification: (chain_id, residue_seq).
pub type ResidueKey = (String, i32);

/// A contact pair between two residues from different chains.
pub type ContactPair = (ResidueKey, ResidueKey);

/// Result of contact comparison between model and native.
#[derive(Debug, Clone)]
pub struct ContactResult {
    /// Fraction of native contacts preserved in the model.
    pub fnat: f64,
    /// Fraction of model contacts that are non-native.
    pub fnonnat: f64,
    /// F1 score combining precision and recall of contacts.
    pub f1: f64,
    /// Number of contacts in the native structure.
    pub num_native_contacts: usize,
    /// Number of contacts in the model structure.
    pub num_model_contacts: usize,
}

/// Find residue-residue contacts at an interface between two chains.
///
/// A contact exists when any pair of heavy atoms from residues in different
/// chains are within the distance threshold.
pub fn find_interface_contacts(
    structure: &PdbStructure,
    chain_a: &str,
    chain_b: &str,
    threshold: f64,
) -> HashSet<ContactPair> {
    let threshold_sq = threshold * threshold;

    let atoms_a: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_a && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let atoms_b: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_b && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let mut contacts = HashSet::new();

    for atom_a in &atoms_a {
        for atom_b in &atoms_b {
            let dx = atom_a.x - atom_b.x;
            let dy = atom_a.y - atom_b.y;
            let dz = atom_a.z - atom_b.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq < threshold_sq {
                let key_a = (atom_a.chain_id.clone(), atom_a.residue_seq);
                let key_b = (atom_b.chain_id.clone(), atom_b.residue_seq);
                // Normalize order: smaller chain first
                let pair = if key_a <= key_b {
                    (key_a, key_b)
                } else {
                    (key_b, key_a)
                };
                contacts.insert(pair);
            }
        }
    }

    contacts
}

/// Count atomic clashes (atom pairs closer than 2.0 A) at an interface.
pub fn count_clashes(structure: &PdbStructure, chain_a: &str, chain_b: &str) -> usize {
    let clash_threshold_sq = 2.0 * 2.0;

    let atoms_a: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_a && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let atoms_b: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_b && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let mut count = 0;
    for atom_a in &atoms_a {
        for atom_b in &atoms_b {
            let dx = atom_a.x - atom_b.x;
            let dy = atom_a.y - atom_b.y;
            let dz = atom_a.z - atom_b.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < clash_threshold_sq {
                count += 1;
            }
        }
    }
    count
}

/// Compare interface contacts between model and native structures.
///
/// The `residue_mapping` maps (model_chain, model_resid) -> (native_chain, native_resid)
/// so that contacts can be compared across different numbering schemes.
pub fn compare_contacts(
    native_contacts: &HashSet<ContactPair>,
    model_contacts: &HashSet<ContactPair>,
    residue_mapping: &std::collections::HashMap<ResidueKey, ResidueKey>,
) -> ContactResult {
    if native_contacts.is_empty() {
        return ContactResult {
            fnat: 0.0,
            fnonnat: 0.0,
            f1: 0.0,
            num_native_contacts: 0,
            num_model_contacts: model_contacts.len(),
        };
    }

    // Map model contacts to native residue keys
    let mapped_model_contacts: HashSet<ContactPair> = model_contacts
        .iter()
        .filter_map(|(ka, kb)| {
            let mapped_a = residue_mapping.get(ka)?;
            let mapped_b = residue_mapping.get(kb)?;
            let pair = if mapped_a <= mapped_b {
                (mapped_a.clone(), mapped_b.clone())
            } else {
                (mapped_b.clone(), mapped_a.clone())
            };
            Some(pair)
        })
        .collect();

    let true_positives = native_contacts.intersection(&mapped_model_contacts).count();
    let false_negatives = native_contacts.len() - true_positives;
    let false_positives = mapped_model_contacts
        .iter()
        .filter(|c| !native_contacts.contains(c))
        .count();

    let fnat = true_positives as f64 / native_contacts.len() as f64;

    let fnonnat = if mapped_model_contacts.is_empty() {
        0.0
    } else {
        false_positives as f64 / mapped_model_contacts.len() as f64
    };

    let f1 = if true_positives == 0 {
        0.0
    } else {
        2.0 * true_positives as f64
            / (2.0 * true_positives as f64 + false_positives as f64 + false_negatives as f64)
    };

    ContactResult {
        fnat,
        fnonnat,
        f1,
        num_native_contacts: native_contacts.len(),
        num_model_contacts: model_contacts.len(),
    }
}

/// Find interface residues: residues from one chain within a distance of another chain.
pub fn find_interface_residues(
    structure: &PdbStructure,
    chain_a: &str,
    chain_b: &str,
    threshold: f64,
) -> HashSet<ResidueKey> {
    let threshold_sq = threshold * threshold;

    let atoms_a: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_a && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let atoms_b: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|a| a.chain_id == chain_b && !a.is_hetatm && !is_hydrogen(a))
        .collect();

    let mut interface_residues = HashSet::new();

    for atom_a in &atoms_a {
        for atom_b in &atoms_b {
            let dx = atom_a.x - atom_b.x;
            let dy = atom_a.y - atom_b.y;
            let dz = atom_a.z - atom_b.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq < threshold_sq {
                interface_residues.insert((atom_a.chain_id.clone(), atom_a.residue_seq));
                interface_residues.insert((atom_b.chain_id.clone(), atom_b.residue_seq));
            }
        }
    }

    interface_residues
}

fn is_hydrogen(atom: &Atom) -> bool {
    atom.element.trim() == "H" || (atom.element.is_empty() && atom.name.trim().starts_with('H'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn make_atom(
        serial: i32,
        name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_find_contacts_close() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_atom(1, "CA", "A", 1, 0.0, 0.0, 0.0),
            make_atom(2, "CA", "B", 1, 3.0, 0.0, 0.0), // 3.0 A apart
        ];

        let contacts = find_interface_contacts(&structure, "A", "B", 5.0);
        assert_eq!(contacts.len(), 1);
    }

    #[test]
    fn test_find_contacts_far() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_atom(1, "CA", "A", 1, 0.0, 0.0, 0.0),
            make_atom(2, "CA", "B", 1, 10.0, 0.0, 0.0), // 10.0 A apart
        ];

        let contacts = find_interface_contacts(&structure, "A", "B", 5.0);
        assert_eq!(contacts.len(), 0);
    }

    #[test]
    fn test_compare_identical_contacts() {
        let mut contacts = HashSet::new();
        contacts.insert((("A".to_string(), 1), ("B".to_string(), 1)));
        contacts.insert((("A".to_string(), 2), ("B".to_string(), 2)));

        // Identity mapping
        let mut mapping = std::collections::HashMap::new();
        mapping.insert(("A".to_string(), 1), ("A".to_string(), 1));
        mapping.insert(("A".to_string(), 2), ("A".to_string(), 2));
        mapping.insert(("B".to_string(), 1), ("B".to_string(), 1));
        mapping.insert(("B".to_string(), 2), ("B".to_string(), 2));

        let result = compare_contacts(&contacts, &contacts, &mapping);
        assert!((result.fnat - 1.0).abs() < 1e-10);
        assert!((result.fnonnat - 0.0).abs() < 1e-10);
        assert!((result.f1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_interface_residues() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_atom(1, "CA", "A", 1, 0.0, 0.0, 0.0),
            make_atom(2, "CA", "A", 2, 3.8, 0.0, 0.0),
            make_atom(3, "CA", "B", 1, 4.0, 0.0, 0.0),
            make_atom(4, "CA", "B", 2, 20.0, 0.0, 0.0),
        ];

        let iface = find_interface_residues(&structure, "A", "B", 10.0);
        // A:1 is within 4.0 of B:1, A:2 is within 0.2 of B:1
        assert!(iface.contains(&("A".to_string(), 1)));
        assert!(iface.contains(&("A".to_string(), 2)));
        assert!(iface.contains(&("B".to_string(), 1)));
        // B:2 is 20 A away from everything in A
        assert!(!iface.contains(&("B".to_string(), 2)));
    }
}
