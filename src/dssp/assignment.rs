//! Main DSSP secondary structure assignment algorithm.
//!
//! This module implements the complete DSSP 4 algorithm:
//! 1. Extract backbone atoms
//! 2. Compute virtual hydrogens where needed
//! 3. Detect hydrogen bonds using Kabsch-Sander energy
//! 4. Identify helix patterns (H, G, I)
//! 5. Identify beta-sheet patterns (E, B)
//! 6. Detect turns (T) and bends (S)
//! 7. Detect PPII/κ-helix (P) in coil regions
//! 8. Assign coil (C) to remaining residues

use super::dihedral::{BackboneDihedrals, calculate_all_dihedrals};
use super::hbond::{
    BackboneAtoms, compute_all_virtual_hydrogens, create_hbond_matrix, detect_hydrogen_bonds,
    extract_backbone_atoms,
};
use super::patterns::{
    HelixType, detect_bends, detect_beta_bridges, detect_helix_patterns, detect_n_turns,
    detect_ppii_regions, group_bridges_into_ladders,
};
use super::types::{ResidueSSAssignment, SecondaryStructure, SecondaryStructureAssignment};
use crate::records::Atom;

/// Assigns secondary structure to a protein structure using the DSSP algorithm.
///
/// This function performs the complete DSSP analysis on a list of atoms
/// and returns per-residue secondary structure assignments.
///
/// # Arguments
/// * `atoms` - Slice of atoms from a PDB structure
///
/// # Returns
/// A `SecondaryStructureAssignment` containing per-residue assignments and statistics.
pub fn assign_secondary_structure(atoms: &[Atom]) -> SecondaryStructureAssignment {
    let mut result = SecondaryStructureAssignment::new();

    if atoms.is_empty() {
        result.add_warning("No atoms provided".to_string());
        return result;
    }

    // Step 1: Extract backbone atoms
    let mut residues = extract_backbone_atoms(atoms);

    if residues.is_empty() {
        result.add_warning("No residues found with backbone atoms".to_string());
        return result;
    }

    // Check for CA-only structures
    let complete_count = residues.iter().filter(|r| r.is_complete()).count();
    if complete_count == 0 {
        result.add_warning("Structure appears to be CA-only; returning all coil".to_string());
        for res in &residues {
            result.residue_assignments.push(ResidueSSAssignment::new(
                res.chain_id.clone(),
                res.residue_seq,
                res.residue_name.clone(),
                res.ins_code,
                SecondaryStructure::Coil,
            ));
        }
        result.compute_statistics();
        return result;
    }

    let num_residues = residues.len();

    // Step 2: Compute virtual hydrogens
    compute_all_virtual_hydrogens(&mut residues);

    // Step 3: Detect hydrogen bonds
    let hbonds = detect_hydrogen_bonds(&residues);
    let hbond_matrix = create_hbond_matrix(&hbonds);

    // Step 4: Calculate backbone dihedrals for PPII detection
    let dihedrals = calculate_all_dihedrals(&residues);

    // Initialize all residues to Coil
    let mut assignments = vec![SecondaryStructure::Coil; num_residues];

    // Step 5: Detect and assign helices (priority: α > π > 3₁₀)
    // Higher priority helices can override lower priority ones
    assign_helices(
        &mut assignments,
        &hbond_matrix,
        num_residues,
        &residues,
        HelixType::Alpha,
    );
    assign_helices(
        &mut assignments,
        &hbond_matrix,
        num_residues,
        &residues,
        HelixType::Pi,
    );
    assign_helices(
        &mut assignments,
        &hbond_matrix,
        num_residues,
        &residues,
        HelixType::Helix310,
    );

    // Step 6: Detect and assign beta-sheets
    assign_beta_structures(&mut assignments, &hbond_matrix, num_residues, &residues);

    // Step 7: Detect and assign turns
    assign_turns(&mut assignments, &hbond_matrix, num_residues, &residues);

    // Step 8: Detect and assign bends
    assign_bends(&mut assignments, &residues);

    // Step 9: Detect and assign PPII/κ-helix (only to coil residues)
    assign_ppii(&mut assignments, &dihedrals);

    // Build the result
    for (i, res) in residues.iter().enumerate() {
        result.residue_assignments.push(ResidueSSAssignment::new(
            res.chain_id.clone(),
            res.residue_seq,
            res.residue_name.clone(),
            res.ins_code,
            assignments[i],
        ));
    }

    result.compute_statistics();
    result
}

/// Assigns helix secondary structure based on detected patterns.
#[allow(clippy::needless_range_loop)]
fn assign_helices(
    assignments: &mut [SecondaryStructure],
    hbond_matrix: &std::collections::HashMap<(usize, usize), f64>,
    num_residues: usize,
    residues: &[BackboneAtoms],
    helix_type: HelixType,
) {
    let segments = detect_helix_patterns(helix_type, hbond_matrix, num_residues, residues);
    let ss = helix_type.to_secondary_structure();

    for segment in segments {
        for i in segment.start..=segment.end.min(num_residues - 1) {
            // Only assign if current assignment is Coil or a lower-priority helix
            match assignments[i] {
                SecondaryStructure::Coil | SecondaryStructure::Turn | SecondaryStructure::Bend => {
                    assignments[i] = ss;
                }
                SecondaryStructure::Helix310 => {
                    // α and π helices override 3₁₀
                    if helix_type != HelixType::Helix310 {
                        assignments[i] = ss;
                    }
                }
                _ => {
                    // Don't override existing assignments
                }
            }
        }
    }
}

/// Assigns beta-sheet secondary structure.
fn assign_beta_structures(
    assignments: &mut [SecondaryStructure],
    hbond_matrix: &std::collections::HashMap<(usize, usize), f64>,
    num_residues: usize,
    residues: &[BackboneAtoms],
) {
    let bridges = detect_beta_bridges(hbond_matrix, num_residues, residues);
    let ladders = group_bridges_into_ladders(&bridges);

    // Assign E (extended strand) for residues in ladders with >1 bridge
    // Assign B (isolated bridge) for single bridges
    for ladder in &ladders {
        let is_ladder = ladder.bridges.len() > 1;

        for bridge in &ladder.bridges {
            let ss = if is_ladder {
                SecondaryStructure::ExtendedStrand
            } else {
                SecondaryStructure::BetaBridge
            };

            // Only assign if not already assigned to helix
            if !assignments[bridge.residue1].is_helix() {
                assignments[bridge.residue1] = ss;
            }
            if !assignments[bridge.residue2].is_helix() {
                assignments[bridge.residue2] = ss;
            }
        }
    }
}

/// Assigns turn secondary structure.
fn assign_turns(
    assignments: &mut [SecondaryStructure],
    hbond_matrix: &std::collections::HashMap<(usize, usize), f64>,
    num_residues: usize,
    residues: &[BackboneAtoms],
) {
    // Detect turns of length 3, 4, and 5
    for n in [3, 4, 5] {
        let turns = detect_n_turns(n, hbond_matrix, num_residues);

        for i in 0..num_residues {
            if turns[i] {
                // Assign T to residues in the turn (i to i+n-1)
                for j in i..(i + n).min(num_residues) {
                    // Check chain continuity
                    if j < residues.len() && residues[j].chain_id != residues[i].chain_id {
                        break;
                    }

                    // Only assign if currently Coil
                    if assignments[j] == SecondaryStructure::Coil {
                        assignments[j] = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }
}

/// Assigns bend secondary structure.
fn assign_bends(assignments: &mut [SecondaryStructure], residues: &[BackboneAtoms]) {
    let bends = detect_bends(residues);

    for (i, &is_bend) in bends.iter().enumerate() {
        if is_bend {
            // Only assign if currently Coil
            if assignments[i] == SecondaryStructure::Coil {
                assignments[i] = SecondaryStructure::Bend;
            }
        }
    }
}

/// Assigns PPII/κ-helix secondary structure.
fn assign_ppii(assignments: &mut [SecondaryStructure], dihedrals: &[BackboneDihedrals]) {
    let ppii = detect_ppii_regions(dihedrals, assignments);

    for (i, &is_ppii) in ppii.iter().enumerate() {
        if is_ppii {
            // Only assign if currently Coil
            if assignments[i] == SecondaryStructure::Coil {
                assignments[i] = SecondaryStructure::KappaHelix;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    #[allow(clippy::too_many_arguments)]
    fn create_test_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
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
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: name.chars().next().unwrap().to_string(),
        }
    }

    #[test]
    fn test_empty_atoms() {
        let atoms: Vec<Atom> = vec![];
        let result = assign_secondary_structure(&atoms);
        assert!(result.is_empty());
        assert!(result.warnings.iter().any(|w| w.contains("No atoms")));
    }

    #[test]
    fn test_ca_only_structure() {
        let atoms = vec![
            create_test_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_test_atom(2, "CA", "GLY", "A", 2, 3.8, 0.0, 0.0),
            create_test_atom(3, "CA", "VAL", "A", 3, 7.6, 0.0, 0.0),
        ];

        let result = assign_secondary_structure(&atoms);
        assert_eq!(result.len(), 3);
        assert!(result.warnings.iter().any(|w| w.contains("CA-only")));

        // All should be coil for CA-only
        for assignment in &result.residue_assignments {
            assert_eq!(assignment.ss, SecondaryStructure::Coil);
        }
    }

    #[test]
    fn test_single_residue() {
        let atoms = vec![
            create_test_atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_test_atom(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0),
            create_test_atom(3, "C", "ALA", "A", 1, 2.5, 1.0, 0.0),
            create_test_atom(4, "O", "ALA", "A", 1, 2.5, 2.2, 0.0),
        ];

        let result = assign_secondary_structure(&atoms);
        assert_eq!(result.len(), 1);
        assert_eq!(result.residue_assignments[0].residue_name, "ALA");
    }

    #[test]
    fn test_statistics_calculation() {
        let mut result = SecondaryStructureAssignment::new();

        // Add 3 helix, 2 sheet, 5 coil
        for i in 0..3 {
            result.residue_assignments.push(ResidueSSAssignment::new(
                "A".to_string(),
                i,
                "ALA".to_string(),
                None,
                SecondaryStructure::AlphaHelix,
            ));
        }
        for i in 3..5 {
            result.residue_assignments.push(ResidueSSAssignment::new(
                "A".to_string(),
                i,
                "VAL".to_string(),
                None,
                SecondaryStructure::ExtendedStrand,
            ));
        }
        for i in 5..10 {
            result.residue_assignments.push(ResidueSSAssignment::new(
                "A".to_string(),
                i,
                "GLY".to_string(),
                None,
                SecondaryStructure::Coil,
            ));
        }

        result.compute_statistics();

        assert_eq!(result.helix_count, 3);
        assert_eq!(result.sheet_count, 2);
        assert_eq!(result.coil_count, 5);
        assert!((result.helix_fraction - 0.3).abs() < 0.01);
        assert!((result.sheet_fraction - 0.2).abs() < 0.01);
        assert!((result.coil_fraction - 0.5).abs() < 0.01);
    }
}
