//! Main DSSP secondary structure assignment algorithm.
//!
//! This module implements the DSSP 4 algorithm, following mkdssp's rules
//! and precedence:
//! 1. Extract backbone atoms
//! 2. Compute virtual hydrogens where needed (H = N + unit(C_prev − O_prev))
//! 3. Detect hydrogen bonds using Kabsch-Sander energy
//! 4. Identify beta bridges and ladders, joined across beta-bulges (E, B)
//! 5. Identify helices from consecutive n-turns: α (H), then 3₁₀ (G) and
//!    π (I) where free; π may replace α, as in DSSP 4
//! 6. Detect turns (T) and bends (S) in the remaining coil
//! 7. Detect PPII/κ-helix (P) in coil regions
//! 8. Assign coil (C) to remaining residues

use super::dihedral::{BackboneDihedrals, calculate_all_dihedrals};
use super::hbond::{
    BackboneAtoms, compute_all_virtual_hydrogens, create_hbond_matrix, detect_hydrogen_bonds,
    extract_backbone_atoms,
};
use super::patterns::{
    HelixType, detect_bends, detect_beta_bridges, detect_helix_starts, detect_ppii_regions,
    group_bridges_into_ladders,
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

    // Steps 5-7 follow mkdssp's order and precedence:
    // E/B first, then H (which may replace them), then G and I only where the
    // residues are still free (I may also replace H, as in DSSP 4), then T.
    assign_beta_structures(&mut assignments, &hbond_matrix, num_residues, &residues);

    let starts = |n| detect_helix_starts(n, &hbond_matrix, &residues);
    let (starts3, starts4, starts5) = (starts(3), starts(4), starts(5));

    assign_helices(&mut assignments, &starts4, HelixType::Alpha, |_| true);
    assign_helices(&mut assignments, &starts3, HelixType::Helix310, |ss| {
        matches!(ss, SecondaryStructure::Coil | SecondaryStructure::Helix310)
    });
    assign_helices(&mut assignments, &starts5, HelixType::Pi, |ss| {
        matches!(
            ss,
            SecondaryStructure::Coil | SecondaryStructure::PiHelix | SecondaryStructure::AlphaHelix
        )
    });

    assign_turns(&mut assignments, [&starts3, &starts4, &starts5]);

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

/// Assigns a helix type from n-turn start flags (DSSP rule).
///
/// Two consecutive n-turns starting at i-1 and i make residues i..i+n-1 a
/// helix, provided `can_replace` accepts every residue currently there.
fn assign_helices(
    assignments: &mut [SecondaryStructure],
    starts: &[bool],
    helix_type: HelixType,
    can_replace: impl Fn(SecondaryStructure) -> bool,
) {
    let n = helix_type.offset();
    let ss = helix_type.to_secondary_structure();
    for i in 1..starts.len() {
        if !(starts[i - 1] && starts[i]) || i + n > assignments.len() {
            continue;
        }
        let span = i..i + n;
        if assignments[span.clone()].iter().all(|&s| can_replace(s)) {
            assignments[span].fill(ss);
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
    let ladders = group_bridges_into_ladders(&bridges, residues);

    // A ladder with more than one bridge (including bulge-joined ladders)
    // makes every residue across both strands E; a lone bridge is B.
    for ladder in &ladders {
        let ss = if ladder.bridge_count() > 1 {
            SecondaryStructure::ExtendedStrand
        } else {
            SecondaryStructure::BetaBridge
        };
        let (i_first, i_last) = ladder.i_span();
        let (j_first, j_last) = ladder.j_span();
        for r in (i_first..=i_last).chain(j_first..=j_last) {
            // A residue in a ladder stays E even if it also has an isolated bridge.
            if assignments[r] != SecondaryStructure::ExtendedStrand {
                assignments[r] = ss;
            }
        }
    }
}

/// Assigns turns (T): a residue still in coil lies inside an n-turn, i.e. an
/// n-turn starts at one of the n-1 residues before it.
fn assign_turns(assignments: &mut [SecondaryStructure], starts: [&Vec<bool>; 3]) {
    for i in 0..assignments.len() {
        if assignments[i] != SecondaryStructure::Coil {
            continue;
        }
        let in_turn = starts
            .iter()
            .zip(3..=5usize)
            .any(|(flags, n)| (1..n).any(|k| i >= k && flags[i - k]));
        if in_turn {
            assignments[i] = SecondaryStructure::Turn;
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

    fn flags(len: usize, at: &[usize]) -> Vec<bool> {
        (0..len).map(|i| at.contains(&i)).collect()
    }

    #[test]
    fn test_helix_needs_two_consecutive_turns() {
        // 4-turns starting at 2 and 3 make residues 3..=6 an alpha helix.
        let mut ss = vec![SecondaryStructure::Coil; 12];
        assign_helices(&mut ss, &flags(12, &[2, 3]), HelixType::Alpha, |_| true);
        let helix: Vec<usize> = (0..12).filter(|&i| ss[i].is_helix()).collect();
        assert_eq!(helix, vec![3, 4, 5, 6]);

        // A single turn is not a helix.
        let mut ss = vec![SecondaryStructure::Coil; 12];
        assign_helices(&mut ss, &flags(12, &[2]), HelixType::Alpha, |_| true);
        assert!(ss.iter().all(|s| *s == SecondaryStructure::Coil));
    }

    #[test]
    fn test_turn_covers_residues_inside_the_turn() {
        // A single 4-turn starting at 2 marks residues 3, 4 and 5 as T.
        let mut ss = vec![SecondaryStructure::Coil; 12];
        let none = flags(12, &[]);
        assign_turns(&mut ss, [&none, &flags(12, &[2]), &none]);
        let turn: Vec<usize> = (0..12)
            .filter(|&i| ss[i] == SecondaryStructure::Turn)
            .collect();
        assert_eq!(turn, vec![3, 4, 5]);
    }

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
