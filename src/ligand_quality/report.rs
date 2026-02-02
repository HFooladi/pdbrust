//! Ligand pose quality report generation.
//!
//! This module combines clash detection and volume overlap calculation
//! into comprehensive quality reports.

use super::clash::{AtomClash, detect_clashes, detect_cofactor_clashes, find_min_distance};
use super::overlap::{calculate_volume_overlap, filter_nearby_protein_atoms};
use super::{MAX_VOLUME_OVERLAP_PCT, WATER_RESIDUES};
use crate::core::PdbStructure;
use crate::records::Atom;
use std::collections::HashSet;

/// Comprehensive report on ligand pose quality.
///
/// Contains results from PoseBusters-style geometry checks including
/// steric clash detection and volume overlap calculation.
#[derive(Debug, Clone)]
pub struct LigandPoseReport {
    /// Residue name of the ligand (e.g., "LIG", "ATP").
    pub ligand_name: String,

    /// Chain ID where the ligand is located.
    pub ligand_chain_id: String,

    /// Residue sequence number of the ligand.
    pub ligand_residue_seq: i32,

    /// Number of atoms in the ligand.
    pub ligand_atom_count: usize,

    // Distance-based clash detection
    /// Minimum distance between any ligand atom and any protein atom (Å).
    pub min_protein_ligand_distance: f64,

    /// List of detected steric clashes with protein atoms.
    pub clashes: Vec<AtomClash>,

    /// Whether any protein clashes were detected.
    pub has_protein_clash: bool,

    /// Total number of protein-ligand clashes.
    pub num_clashes: usize,

    /// Severity of the worst clash (expected_min / actual distance).
    /// Higher values indicate more severe clashes.
    pub worst_clash_severity: f64,

    // Volume overlap
    /// Percentage of ligand volume overlapping with protein (0-100%).
    pub protein_volume_overlap_pct: f64,

    // Cofactor checks
    /// List of detected clashes with other HETATM atoms (cofactors).
    pub cofactor_clashes: Vec<AtomClash>,

    /// Whether any cofactor clashes were detected.
    pub has_cofactor_clash: bool,

    /// Number of clashes with cofactors.
    pub num_cofactor_clashes: usize,

    // Overall verdicts
    /// Whether the pose passes the minimum distance check.
    /// True if min_protein_ligand_distance > 0.75 × sum(vdW radii).
    pub passes_distance_check: bool,

    /// Whether the pose passes the volume overlap check.
    /// True if protein_volume_overlap_pct < 7.5%.
    pub passes_overlap_check: bool,

    /// Whether the pose passes all geometry checks.
    /// True only if both distance and overlap checks pass.
    pub is_geometry_valid: bool,
}

impl LigandPoseReport {
    /// Create a new report with default (empty) values.
    pub fn new(ligand_name: &str, chain_id: &str, residue_seq: i32) -> Self {
        Self {
            ligand_name: ligand_name.to_string(),
            ligand_chain_id: chain_id.to_string(),
            ligand_residue_seq: residue_seq,
            ligand_atom_count: 0,
            min_protein_ligand_distance: f64::INFINITY,
            clashes: Vec::new(),
            has_protein_clash: false,
            num_clashes: 0,
            worst_clash_severity: 0.0,
            protein_volume_overlap_pct: 0.0,
            cofactor_clashes: Vec::new(),
            has_cofactor_clash: false,
            num_cofactor_clashes: 0,
            passes_distance_check: true,
            passes_overlap_check: true,
            is_geometry_valid: true,
        }
    }

    /// Get a summary string of the report.
    pub fn summary(&self) -> String {
        let status = if self.is_geometry_valid {
            "PASS"
        } else {
            "FAIL"
        };

        format!(
            "{} ({}{}): {} - {} clashes, {:.1}% overlap, min_dist={:.2}Å",
            self.ligand_name,
            self.ligand_chain_id,
            self.ligand_residue_seq,
            status,
            self.num_clashes,
            self.protein_volume_overlap_pct,
            self.min_protein_ligand_distance
        )
    }
}

/// Compute a ligand pose quality report for a specific ligand.
///
/// # Arguments
///
/// * `structure` - The PDB structure containing the protein-ligand complex
/// * `ligand_name` - The residue name of the ligand to analyze
///
/// # Returns
///
/// `Some(LigandPoseReport)` if the ligand exists, `None` otherwise.
pub fn compute_ligand_pose_report(
    structure: &PdbStructure,
    ligand_name: &str,
) -> Option<LigandPoseReport> {
    // Find all atoms belonging to this ligand
    let ligand_atoms: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|atom| atom.is_hetatm && atom.residue_name == ligand_name)
        .collect();

    if ligand_atoms.is_empty() {
        return None;
    }

    // Get ligand info from first atom
    let first_lig = ligand_atoms[0];
    let mut report = LigandPoseReport::new(ligand_name, &first_lig.chain_id, first_lig.residue_seq);
    report.ligand_atom_count = ligand_atoms.len();

    // Get protein atoms (non-HETATM)
    let protein_atoms: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|atom| !atom.is_hetatm)
        .collect();

    // Get cofactor atoms (HETATM but not the ligand or water)
    let cofactor_atoms: Vec<&Atom> = structure
        .atoms
        .iter()
        .filter(|atom| {
            atom.is_hetatm
                && atom.residue_name != ligand_name
                && !WATER_RESIDUES.contains(&atom.residue_name.as_str())
        })
        .collect();

    // Build connectivity set from CONECT records
    let connected_pairs = build_connectivity_set(structure);

    // Calculate minimum distance
    report.min_protein_ligand_distance = find_min_distance(&ligand_atoms, &protein_atoms);

    // Detect protein-ligand clashes
    report.clashes = detect_clashes(&ligand_atoms, &protein_atoms, &connected_pairs);
    report.num_clashes = report.clashes.len();
    report.has_protein_clash = !report.clashes.is_empty();
    report.worst_clash_severity = report
        .clashes
        .iter()
        .map(|c| c.severity)
        .fold(0.0f64, f64::max);

    // Detect cofactor clashes
    report.cofactor_clashes =
        detect_cofactor_clashes(&ligand_atoms, &cofactor_atoms, &connected_pairs);
    report.num_cofactor_clashes = report.cofactor_clashes.len();
    report.has_cofactor_clash = !report.cofactor_clashes.is_empty();

    // Calculate volume overlap (filter to nearby atoms for efficiency)
    let nearby_protein = filter_nearby_protein_atoms(&ligand_atoms, &protein_atoms, 8.0);
    let nearby_refs: Vec<&Atom> = nearby_protein.to_vec();
    report.protein_volume_overlap_pct = calculate_volume_overlap(&ligand_atoms, &nearby_refs);

    // Determine verdicts
    report.passes_distance_check = !report.has_protein_clash;
    report.passes_overlap_check = report.protein_volume_overlap_pct < MAX_VOLUME_OVERLAP_PCT;
    report.is_geometry_valid = report.passes_distance_check && report.passes_overlap_check;

    Some(report)
}

/// Compute quality reports for all ligands in a structure.
///
/// # Arguments
///
/// * `structure` - The PDB structure containing protein-ligand complex(es)
///
/// # Returns
///
/// A vector of `LigandPoseReport` for each unique ligand found.
pub fn compute_all_ligand_reports(structure: &PdbStructure) -> Vec<LigandPoseReport> {
    // Get unique ligand identifiers (name + chain + residue_seq)
    let mut seen: HashSet<(String, String, i32)> = HashSet::new();
    let mut ligand_ids: Vec<(String, String, i32)> = Vec::new();

    for atom in &structure.atoms {
        if atom.is_hetatm && !WATER_RESIDUES.contains(&atom.residue_name.as_str()) {
            let id = (
                atom.residue_name.clone(),
                atom.chain_id.clone(),
                atom.residue_seq,
            );
            if !seen.contains(&id) {
                seen.insert(id.clone());
                ligand_ids.push(id);
            }
        }
    }

    // Sort by chain, then residue sequence
    ligand_ids.sort_by(|a, b| match a.1.cmp(&b.1) {
        std::cmp::Ordering::Equal => a.2.cmp(&b.2),
        other => other,
    });

    // Compute reports for each unique ligand instance
    ligand_ids
        .iter()
        .filter_map(|(name, _chain, _seq)| compute_ligand_pose_report(structure, name))
        .collect()
}

/// Build a set of connected atom pairs from CONECT records.
fn build_connectivity_set(structure: &PdbStructure) -> HashSet<(i32, i32)> {
    let mut connected: HashSet<(i32, i32)> = HashSet::new();

    for conect in &structure.connects {
        // Add the primary bond
        connected.insert((conect.atom1, conect.atom2));
        connected.insert((conect.atom2, conect.atom1));

        // Add optional bonds
        if let Some(atom3) = conect.atom3 {
            connected.insert((conect.atom1, atom3));
            connected.insert((atom3, conect.atom1));
        }
        if let Some(atom4) = conect.atom4 {
            connected.insert((conect.atom1, atom4));
            connected.insert((atom4, conect.atom1));
        }
    }

    connected
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
        element: &str,
        is_hetatm: bool,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    fn create_valid_complex() -> PdbStructure {
        let mut structure = PdbStructure::new();

        // Protein backbone
        structure.atoms.push(create_test_atom(
            1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", false,
        ));
        structure.atoms.push(create_test_atom(
            2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, "C", false,
        ));
        structure.atoms.push(create_test_atom(
            3, "C", "ALA", "A", 1, 3.0, 0.0, 0.0, "C", false,
        ));
        structure.atoms.push(create_test_atom(
            4, "O", "ALA", "A", 1, 3.0, 1.2, 0.0, "O", false,
        ));

        // Ligand well away from protein
        structure.atoms.push(create_test_atom(
            10, "C1", "LIG", "A", 100, 10.0, 10.0, 10.0, "C", true,
        ));
        structure.atoms.push(create_test_atom(
            11, "O1", "LIG", "A", 100, 11.0, 10.0, 10.0, "O", true,
        ));
        structure.atoms.push(create_test_atom(
            12, "N1", "LIG", "A", 100, 10.0, 11.0, 10.0, "N", true,
        ));

        structure
    }

    fn create_clashing_complex() -> PdbStructure {
        let mut structure = PdbStructure::new();

        // Protein
        structure.atoms.push(create_test_atom(
            1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", false,
        ));

        // Ligand with atom clashing with protein
        structure.atoms.push(create_test_atom(
            10, "C1", "BAD", "A", 100, 1.5, 0.0, 0.0, "C", true,
        ));

        structure
    }

    #[test]
    fn test_report_new() {
        let report = LigandPoseReport::new("LIG", "A", 100);

        assert_eq!(report.ligand_name, "LIG");
        assert_eq!(report.ligand_chain_id, "A");
        assert_eq!(report.ligand_residue_seq, 100);
        assert!(report.is_geometry_valid);
    }

    #[test]
    fn test_compute_valid_pose_report() {
        let structure = create_valid_complex();
        let report = compute_ligand_pose_report(&structure, "LIG");

        assert!(report.is_some());
        let report = report.unwrap();

        assert_eq!(report.ligand_name, "LIG");
        assert_eq!(report.ligand_atom_count, 3);
        assert!(!report.has_protein_clash);
        assert_eq!(report.num_clashes, 0);
        assert!(report.passes_distance_check);
        assert!(report.is_geometry_valid);
    }

    #[test]
    fn test_compute_clashing_pose_report() {
        let structure = create_clashing_complex();
        let report = compute_ligand_pose_report(&structure, "BAD");

        assert!(report.is_some());
        let report = report.unwrap();

        assert!(report.has_protein_clash);
        assert!(report.num_clashes > 0);
        assert!(!report.passes_distance_check);
        assert!(!report.is_geometry_valid);
    }

    #[test]
    fn test_compute_nonexistent_ligand() {
        let structure = create_valid_complex();
        let report = compute_ligand_pose_report(&structure, "XYZ");

        assert!(report.is_none());
    }

    #[test]
    fn test_compute_all_reports() {
        let structure = create_valid_complex();
        let reports = compute_all_ligand_reports(&structure);

        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].ligand_name, "LIG");
    }

    #[test]
    fn test_report_summary() {
        let structure = create_valid_complex();
        let report = compute_ligand_pose_report(&structure, "LIG").unwrap();

        let summary = report.summary();
        assert!(summary.contains("LIG"));
        assert!(summary.contains("PASS") || summary.contains("FAIL"));
    }

    #[test]
    fn test_connectivity_excludes_bonded() {
        let mut structure = create_clashing_complex();

        // Add CONECT record to mark atoms as bonded
        structure.connects.push(crate::records::Conect {
            atom1: 1,
            atom2: 10,
            atom3: None,
            atom4: None,
        });

        let report = compute_ligand_pose_report(&structure, "BAD");
        assert!(report.is_some());
        let report = report.unwrap();

        // Should have no clashes because atoms are marked as bonded
        assert!(!report.has_protein_clash);
        assert_eq!(report.num_clashes, 0);
    }

    #[test]
    fn test_water_excluded() {
        let mut structure = create_valid_complex();

        // Add water
        structure.atoms.push(create_test_atom(
            20, "O", "HOH", "A", 200, 5.0, 5.0, 5.0, "O", true,
        ));

        let reports = compute_all_ligand_reports(&structure);

        // Water should not appear as a ligand
        assert_eq!(reports.len(), 1);
        assert_eq!(reports[0].ligand_name, "LIG");
    }
}
