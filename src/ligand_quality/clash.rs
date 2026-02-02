//! Steric clash detection for protein-ligand complexes.
//!
//! This module implements van der Waals-based clash detection following
//! the PoseBusters methodology.

use super::CLASH_VDW_MULTIPLIER;
use super::radii::{is_metal, min_contact_distance, vdw_radius};
use crate::records::Atom;

/// Represents a steric clash between two atoms.
///
/// A clash occurs when two non-bonded atoms are closer than the expected
/// minimum distance based on their van der Waals radii.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomClash {
    /// Serial number of the protein/receptor atom.
    pub protein_atom_serial: i32,
    /// Serial number of the ligand atom.
    pub ligand_atom_serial: i32,
    /// Chain ID of the protein atom.
    pub protein_chain_id: String,
    /// Residue name of the protein atom.
    pub protein_residue_name: String,
    /// Residue sequence number of the protein atom.
    pub protein_residue_seq: i32,
    /// Atom name of the protein atom.
    pub protein_atom_name: String,
    /// Element of the protein atom.
    pub protein_element: String,
    /// Atom name of the ligand atom.
    pub ligand_atom_name: String,
    /// Element of the ligand atom.
    pub ligand_element: String,
    /// Actual distance between atoms in Angstroms.
    pub distance: f64,
    /// Expected minimum distance (0.75 × sum of vdW radii).
    pub expected_min_distance: f64,
    /// Severity ratio: expected_min_distance / distance.
    /// Higher values indicate more severe clashes.
    pub severity: f64,
}

impl AtomClash {
    /// Create a new AtomClash from two atoms and their distance.
    pub fn new(
        protein_atom: &Atom,
        ligand_atom: &Atom,
        distance: f64,
        expected_min_distance: f64,
    ) -> Self {
        Self {
            protein_atom_serial: protein_atom.serial,
            ligand_atom_serial: ligand_atom.serial,
            protein_chain_id: protein_atom.chain_id.clone(),
            protein_residue_name: protein_atom.residue_name.clone(),
            protein_residue_seq: protein_atom.residue_seq,
            protein_atom_name: protein_atom.name.clone(),
            protein_element: protein_atom.element.clone(),
            ligand_atom_name: ligand_atom.name.clone(),
            ligand_element: ligand_atom.element.clone(),
            distance,
            expected_min_distance,
            severity: expected_min_distance / distance,
        }
    }
}

/// Detect steric clashes between ligand atoms and protein atoms.
///
/// Uses van der Waals radii to determine the expected minimum distance
/// between atom pairs. A clash is detected when the actual distance is
/// less than `CLASH_VDW_MULTIPLIER` (0.75) times the sum of vdW radii.
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand
/// * `protein_atoms` - Atoms belonging to the protein
/// * `connected_pairs` - Set of (serial1, serial2) pairs that are bonded
///   (used to exclude covalent ligands)
///
/// # Returns
///
/// A vector of `AtomClash` structs describing each detected clash.
pub fn detect_clashes(
    ligand_atoms: &[&Atom],
    protein_atoms: &[&Atom],
    connected_pairs: &std::collections::HashSet<(i32, i32)>,
) -> Vec<AtomClash> {
    let mut clashes = Vec::new();

    for lig_atom in ligand_atoms {
        for prot_atom in protein_atoms {
            // Skip if atoms are bonded (covalent ligand case)
            if is_bonded(lig_atom.serial, prot_atom.serial, connected_pairs) {
                continue;
            }

            let distance = calculate_distance(lig_atom, prot_atom);

            // Calculate expected minimum distance
            let expected_min =
                min_contact_distance(&lig_atom.element, &prot_atom.element, CLASH_VDW_MULTIPLIER);

            // Check for clash
            if distance < expected_min {
                clashes.push(AtomClash::new(prot_atom, lig_atom, distance, expected_min));
            }
        }
    }

    // Sort clashes by severity (most severe first)
    clashes.sort_by(|a, b| {
        b.severity
            .partial_cmp(&a.severity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    clashes
}

/// Detect clashes between ligand atoms and cofactor/other HETATM atoms.
///
/// Uses the same methodology as protein clash detection but distinguishes
/// between organic cofactors (vdW-based) and metal ions (covalent-based
/// for coordination).
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand being evaluated
/// * `cofactor_atoms` - Other HETATM atoms (excluding water and the ligand itself)
/// * `connected_pairs` - Set of (serial1, serial2) pairs that are bonded
///
/// # Returns
///
/// A vector of `AtomClash` structs describing each detected clash.
pub fn detect_cofactor_clashes(
    ligand_atoms: &[&Atom],
    cofactor_atoms: &[&Atom],
    connected_pairs: &std::collections::HashSet<(i32, i32)>,
) -> Vec<AtomClash> {
    let mut clashes = Vec::new();

    for lig_atom in ligand_atoms {
        for cof_atom in cofactor_atoms {
            // Skip if atoms are bonded
            if is_bonded(lig_atom.serial, cof_atom.serial, connected_pairs) {
                continue;
            }

            let distance = calculate_distance(lig_atom, cof_atom);

            // For metal cofactors, use a more permissive threshold
            // as coordination distances can be shorter
            let scale = if is_metal(&cof_atom.element) || is_metal(&lig_atom.element) {
                0.5 // More permissive for metals
            } else {
                CLASH_VDW_MULTIPLIER
            };

            let expected_min = min_contact_distance(&lig_atom.element, &cof_atom.element, scale);

            if distance < expected_min {
                clashes.push(AtomClash::new(cof_atom, lig_atom, distance, expected_min));
            }
        }
    }

    // Sort by severity
    clashes.sort_by(|a, b| {
        b.severity
            .partial_cmp(&a.severity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    clashes
}

/// Find the minimum distance between any ligand atom and any protein atom.
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand
/// * `protein_atoms` - Atoms belonging to the protein
///
/// # Returns
///
/// The minimum distance in Angstroms, or `f64::INFINITY` if either set is empty.
pub fn find_min_distance(ligand_atoms: &[&Atom], protein_atoms: &[&Atom]) -> f64 {
    let mut min_dist = f64::INFINITY;

    for lig_atom in ligand_atoms {
        for prot_atom in protein_atoms {
            let distance = calculate_distance(lig_atom, prot_atom);
            if distance < min_dist {
                min_dist = distance;
            }
        }
    }

    min_dist
}

/// Calculate the Euclidean distance between two atoms.
#[inline]
fn calculate_distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let dx = atom1.x - atom2.x;
    let dy = atom1.y - atom2.y;
    let dz = atom1.z - atom2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Check if two atoms are bonded based on CONECT records.
#[inline]
fn is_bonded(
    serial1: i32,
    serial2: i32,
    connected_pairs: &std::collections::HashSet<(i32, i32)>,
) -> bool {
    connected_pairs.contains(&(serial1, serial2)) || connected_pairs.contains(&(serial2, serial1))
}

/// Get the radius for volume calculations (scaled vdW radius).
pub fn get_volume_radius(element: &str, scale: f64) -> f64 {
    vdw_radius(element) * scale
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(
        serial: i32,
        x: f64,
        y: f64,
        z: f64,
        element: &str,
        residue_name: &str,
    ) -> Atom {
        Atom {
            serial,
            name: "C".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    #[test]
    fn test_atom_clash_creation() {
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let lig_atom = create_test_atom(2, 2.0, 0.0, 0.0, "C", "LIG");

        let clash = AtomClash::new(&prot_atom, &lig_atom, 2.0, 2.55);

        assert_eq!(clash.protein_atom_serial, 1);
        assert_eq!(clash.ligand_atom_serial, 2);
        assert_eq!(clash.distance, 2.0);
        assert_eq!(clash.expected_min_distance, 2.55);
        assert!((clash.severity - 1.275).abs() < 1e-10); // 2.55 / 2.0
    }

    #[test]
    fn test_detect_clashes_no_clash() {
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let lig_atom = create_test_atom(2, 5.0, 0.0, 0.0, "C", "LIG"); // 5 Å apart

        let prot_atoms: Vec<&Atom> = vec![&prot_atom];
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let connected = std::collections::HashSet::new();

        let clashes = detect_clashes(&lig_atoms, &prot_atoms, &connected);

        assert!(clashes.is_empty());
    }

    #[test]
    fn test_detect_clashes_with_clash() {
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let lig_atom = create_test_atom(2, 2.0, 0.0, 0.0, "C", "LIG"); // 2 Å apart

        // C-C expected min: (1.70 + 1.70) * 0.75 = 2.55 Å
        // 2.0 < 2.55, so this should be a clash

        let prot_atoms: Vec<&Atom> = vec![&prot_atom];
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let connected = std::collections::HashSet::new();

        let clashes = detect_clashes(&lig_atoms, &prot_atoms, &connected);

        assert_eq!(clashes.len(), 1);
        assert!((clashes[0].distance - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_detect_clashes_skip_bonded() {
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let lig_atom = create_test_atom(2, 1.5, 0.0, 0.0, "C", "LIG"); // Very close but bonded

        let prot_atoms: Vec<&Atom> = vec![&prot_atom];
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];

        let mut connected = std::collections::HashSet::new();
        connected.insert((1, 2)); // Mark as bonded

        let clashes = detect_clashes(&lig_atoms, &prot_atoms, &connected);

        assert!(clashes.is_empty());
    }

    #[test]
    fn test_find_min_distance() {
        let prot_atom1 = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let prot_atom2 = create_test_atom(2, 10.0, 0.0, 0.0, "C", "ALA");
        let lig_atom = create_test_atom(3, 3.0, 0.0, 0.0, "C", "LIG");

        let prot_atoms: Vec<&Atom> = vec![&prot_atom1, &prot_atom2];
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];

        let min_dist = find_min_distance(&lig_atoms, &prot_atoms);

        assert!((min_dist - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_min_distance_empty() {
        let prot_atoms: Vec<&Atom> = vec![];
        let lig_atoms: Vec<&Atom> = vec![];

        let min_dist = find_min_distance(&lig_atoms, &prot_atoms);

        assert!(min_dist.is_infinite());
    }

    #[test]
    fn test_clash_severity_ordering() {
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C", "ALA");
        let lig_atom1 = create_test_atom(2, 2.0, 0.0, 0.0, "C", "LIG"); // Moderate clash
        let lig_atom2 = create_test_atom(3, 1.5, 0.0, 0.0, "C", "LIG"); // Severe clash

        let prot_atoms: Vec<&Atom> = vec![&prot_atom];
        let lig_atoms: Vec<&Atom> = vec![&lig_atom1, &lig_atom2];
        let connected = std::collections::HashSet::new();

        let clashes = detect_clashes(&lig_atoms, &prot_atoms, &connected);

        // Should be sorted by severity (most severe first)
        assert_eq!(clashes.len(), 2);
        assert!(clashes[0].severity > clashes[1].severity);
    }
}
