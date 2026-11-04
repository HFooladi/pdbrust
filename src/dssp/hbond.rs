//! Hydrogen bond detection for DSSP secondary structure assignment.
//!
//! This module implements the Kabsch-Sander hydrogen bond energy calculation
//! based on the electrostatic model from the original DSSP paper.
//!
//! The energy is calculated as:
//! E = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) kcal/mol
//!
//! A hydrogen bond exists if E < -0.5 kcal/mol.

use crate::records::Atom;
use std::collections::HashMap;

/// Energy threshold for hydrogen bond detection (kcal/mol).
/// Hydrogen bonds are assigned when E < HBOND_THRESHOLD.
pub const HBOND_THRESHOLD: f64 = -0.5;

/// Electrostatic constant for H-bond energy calculation.
/// This gives energy in kcal/mol when distances are in Angstroms.
const ELECTROSTATIC_CONSTANT: f64 = 27.888;

/// N-H distance used for the virtual amide hydrogen (Angstroms).
/// Kabsch & Sander (and mkdssp) use 1.0 Å.
const NH_BOND_LENGTH: f64 = 1.0;

/// Maximum distance between C and N atoms to consider H-bond (Angstroms).
/// This is an initial filter before computing the full energy.
const MAX_CN_DISTANCE: f64 = 6.0;

/// Backbone atom coordinates for a single residue.
#[derive(Debug, Clone)]
pub struct BackboneAtoms {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name
    pub residue_name: String,
    /// Insertion code
    pub ins_code: Option<char>,
    /// Nitrogen atom coordinates (donor)
    pub n: Option<[f64; 3]>,
    /// C-alpha atom coordinates
    pub ca: Option<[f64; 3]>,
    /// Carbonyl carbon atom coordinates
    pub c: Option<[f64; 3]>,
    /// Carbonyl oxygen atom coordinates (acceptor)
    pub o: Option<[f64; 3]>,
    /// Amide hydrogen coordinates (virtual if not present)
    pub h: Option<[f64; 3]>,
}

impl BackboneAtoms {
    /// Creates a new empty BackboneAtoms struct.
    pub fn new(
        chain_id: String,
        residue_seq: i32,
        residue_name: String,
        ins_code: Option<char>,
    ) -> Self {
        Self {
            chain_id,
            residue_seq,
            residue_name,
            ins_code,
            n: None,
            ca: None,
            c: None,
            o: None,
            h: None,
        }
    }

    /// Returns true if all required backbone atoms are present for H-bond donor.
    pub fn is_complete_donor(&self) -> bool {
        self.n.is_some() && self.h.is_some()
    }

    /// Returns true if all required backbone atoms are present for H-bond acceptor.
    pub fn is_complete_acceptor(&self) -> bool {
        self.c.is_some() && self.o.is_some()
    }

    /// Returns true if all backbone atoms are present.
    pub fn is_complete(&self) -> bool {
        self.n.is_some() && self.ca.is_some() && self.c.is_some() && self.o.is_some()
    }

    /// Returns true if this is a proline residue (cannot donate H-bond).
    pub fn is_proline(&self) -> bool {
        self.residue_name == "PRO"
    }
}

/// Represents a hydrogen bond between two residues.
#[derive(Debug, Clone)]
pub struct HydrogenBond {
    /// Donor residue index
    pub donor_index: usize,
    /// Acceptor residue index
    pub acceptor_index: usize,
    /// H-bond energy (kcal/mol)
    pub energy: f64,
}

impl HydrogenBond {
    /// Creates a new hydrogen bond.
    pub fn new(donor_index: usize, acceptor_index: usize, energy: f64) -> Self {
        Self {
            donor_index,
            acceptor_index,
            energy,
        }
    }
}

/// Calculates the Euclidean distance between two 3D points.
#[inline]
fn distance(p1: &[f64; 3], p2: &[f64; 3]) -> f64 {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dz = p2[2] - p1[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Calculates the squared distance between two 3D points.
#[inline]
fn distance_squared(p1: &[f64; 3], p2: &[f64; 3]) -> f64 {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let dz = p2[2] - p1[2];
    dx * dx + dy * dy + dz * dz
}

/// Normalizes a 3D vector to unit length.
#[inline]
fn normalize(v: &[f64; 3]) -> Option<[f64; 3]> {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-10 {
        return None;
    }
    Some([v[0] / len, v[1] / len, v[2] / len])
}

/// Computes a virtual hydrogen position for a residue.
///
/// Following Kabsch & Sander (and mkdssp), the amide H is placed
/// 1.0 Å from N along the previous residue's O→C direction, so
/// that N–H is antiparallel to the preceding C=O:
/// `H = N + 1.0 Å · (C_prev − O_prev) / |C_prev − O_prev|`.
///
/// Returns None for proline, for the first residue in a chain, or when the
/// previous residue lacks C or O.
pub fn compute_virtual_hydrogen(
    backbone: &BackboneAtoms,
    prev_backbone: Option<&BackboneAtoms>,
) -> Option<[f64; 3]> {
    // Proline cannot donate H-bonds
    if backbone.is_proline() {
        return None;
    }

    let n = backbone.n?;
    let prev = prev_backbone?;
    let (prev_c, prev_o) = (prev.c?, prev.o?);

    let oc_dir = normalize(&[
        prev_c[0] - prev_o[0],
        prev_c[1] - prev_o[1],
        prev_c[2] - prev_o[2],
    ])?;

    Some([
        n[0] + NH_BOND_LENGTH * oc_dir[0],
        n[1] + NH_BOND_LENGTH * oc_dir[1],
        n[2] + NH_BOND_LENGTH * oc_dir[2],
    ])
}

/// Calculates the Kabsch-Sander hydrogen bond energy.
///
/// Energy (kcal/mol) = 27.888 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
///
/// Where:
/// - O is the carbonyl oxygen of the acceptor
/// - C is the carbonyl carbon of the acceptor
/// - N is the nitrogen of the donor
/// - H is the amide hydrogen of the donor
///
/// Returns None if the energy calculation is not possible (missing atoms
/// or too far apart).
pub fn calculate_hbond_energy(donor: &BackboneAtoms, acceptor: &BackboneAtoms) -> Option<f64> {
    // Get required atom positions
    let n = donor.n?;
    let h = donor.h?;
    let c = acceptor.c?;
    let o = acceptor.o?;

    // Quick distance check
    let cn_dist_sq = distance_squared(&c, &n);
    if cn_dist_sq > MAX_CN_DISTANCE * MAX_CN_DISTANCE {
        return None;
    }

    // Calculate all required distances
    let r_on = distance(&o, &n);
    let r_ch = distance(&c, &h);
    let r_oh = distance(&o, &h);
    let r_cn = cn_dist_sq.sqrt();

    // Avoid division by zero or very small distances
    if r_on < 0.5 || r_ch < 0.5 || r_oh < 0.5 || r_cn < 0.5 {
        return None;
    }

    // Calculate energy
    let energy = ELECTROSTATIC_CONSTANT * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);

    Some(energy)
}

/// Extracts backbone atoms from a PDB structure's atom list.
///
/// Groups atoms by residue and extracts N, CA, C, O positions.
pub fn extract_backbone_atoms(atoms: &[Atom]) -> Vec<BackboneAtoms> {
    // Group atoms by (chain_id, residue_seq, ins_code)
    let mut residue_map: HashMap<(String, i32, Option<char>), BackboneAtoms> = HashMap::new();

    for atom in atoms {
        // Skip HETATM records for standard residues
        let name = atom.name.trim();
        let key = (atom.chain_id.clone(), atom.residue_seq, atom.ins_code);

        let backbone = residue_map.entry(key).or_insert_with(|| {
            BackboneAtoms::new(
                atom.chain_id.clone(),
                atom.residue_seq,
                atom.residue_name.clone(),
                atom.ins_code,
            )
        });

        let coords = [atom.x, atom.y, atom.z];
        match name {
            "N" => backbone.n = Some(coords),
            "CA" => backbone.ca = Some(coords),
            "C" => backbone.c = Some(coords),
            // Use O preferentially, only use OXT if O is not set
            "O" | "OXT" if backbone.o.is_none() || name == "O" => {
                backbone.o = Some(coords);
            }
            // Use explicit H if present
            "H" | "HN" | "H1" if backbone.h.is_none() => {
                backbone.h = Some(coords);
            }
            _ => {}
        }
    }

    // Sort by chain and residue number
    let mut residues: Vec<BackboneAtoms> = residue_map.into_values().collect();
    residues.sort_by(|a, b| {
        a.chain_id
            .cmp(&b.chain_id)
            .then_with(|| a.residue_seq.cmp(&b.residue_seq))
            .then_with(|| a.ins_code.cmp(&b.ins_code))
    });

    residues
}

/// Maximum C(i-1)–N(i) distance for two residues to count as peptide-bonded
/// (Angstroms), as in mkdssp.
const MAX_PEPTIDE_BOND_LENGTH: f64 = 2.5;

/// Whether `residues[i]` and `residues[i + 1]` are consecutive residues of one
/// chain, i.e. there is no chain break between them.
///
/// Like mkdssp, a break is detected physically: the residues must share a
/// chain and C(i)–N(i+1) must be at most 2.5 Å. If either atom is missing the
/// residue numbering decides instead (consecutive or same number).
pub fn residues_linked(residues: &[BackboneAtoms], i: usize) -> bool {
    let (Some(a), Some(b)) = (residues.get(i), residues.get(i + 1)) else {
        return false;
    };
    if a.chain_id != b.chain_id {
        return false;
    }
    match (a.c, b.n) {
        (Some(c), Some(n)) => distance(&c, &n) <= MAX_PEPTIDE_BOND_LENGTH,
        _ => (0..=1).contains(&(b.residue_seq - a.residue_seq)),
    }
}

/// Whether residues `from..=to` form one unbroken stretch of chain.
pub fn no_chain_break(residues: &[BackboneAtoms], from: usize, to: usize) -> bool {
    (from..to).all(|i| residues_linked(residues, i))
}

/// Computes virtual hydrogens for all residues where needed.
///
/// Modifies the backbone atoms in place to add virtual H positions
/// where explicit H is not present.
pub fn compute_all_virtual_hydrogens(residues: &mut [BackboneAtoms]) {
    for i in 1..residues.len() {
        // Only compute virtual H if no explicit H is present, and only when
        // the previous residue is peptide-bonded to this one.
        if residues[i].h.is_some() || !residues_linked(residues, i - 1) {
            continue;
        }

        let prev = &residues[i - 1];
        if let Some(virtual_h) = compute_virtual_hydrogen(&residues[i], Some(prev)) {
            residues[i].h = Some(virtual_h);
        }
    }
}

/// Detects all hydrogen bonds between residues.
///
/// Returns a vector of HydrogenBond structs for all valid H-bonds.
#[allow(clippy::needless_range_loop)]
pub fn detect_hydrogen_bonds(residues: &[BackboneAtoms]) -> Vec<HydrogenBond> {
    let mut hbonds = Vec::new();

    for i in 0..residues.len() {
        let donor = &residues[i];

        // Skip if donor is not complete or is proline
        if !donor.is_complete_donor() || donor.is_proline() {
            continue;
        }

        for j in 0..residues.len() {
            // A residue cannot H-bond to itself or its immediate neighbors
            // (need at least i+2 for any meaningful H-bond)
            if (i as isize - j as isize).abs() < 2 {
                continue;
            }

            let acceptor = &residues[j];

            // Skip if acceptor is not complete
            if !acceptor.is_complete_acceptor() {
                continue;
            }

            // Skip if different chains (inter-chain H-bonds are handled separately for sheets)
            // For now, we process all pairs

            // Calculate H-bond energy
            if let Some(energy) = calculate_hbond_energy(donor, acceptor) {
                if energy < HBOND_THRESHOLD {
                    hbonds.push(HydrogenBond::new(i, j, energy));
                }
            }
        }
    }

    hbonds
}

/// Creates an H-bond matrix for fast lookup.
///
/// Returns a map where key is (donor_index, acceptor_index) and value is energy.
pub fn create_hbond_matrix(hbonds: &[HydrogenBond]) -> HashMap<(usize, usize), f64> {
    hbonds
        .iter()
        .map(|hb| ((hb.donor_index, hb.acceptor_index), hb.energy))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance() {
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [3.0, 4.0, 0.0];
        assert!((distance(&p1, &p2) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalize() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize(&v).unwrap();
        assert!((n[0] - 0.6).abs() < 1e-10);
        assert!((n[1] - 0.8).abs() < 1e-10);
        assert!(n[2].abs() < 1e-10);

        // Zero vector should return None
        let zero = [0.0, 0.0, 0.0];
        assert!(normalize(&zero).is_none());
    }

    #[test]
    fn test_backbone_atoms() {
        let backbone = BackboneAtoms::new("A".to_string(), 1, "ALA".to_string(), None);
        assert!(!backbone.is_complete());
        assert!(!backbone.is_proline());

        let proline = BackboneAtoms::new("A".to_string(), 2, "PRO".to_string(), None);
        assert!(proline.is_proline());
    }

    #[test]
    fn test_virtual_hydrogen_is_antiparallel_to_previous_carbonyl() {
        // Kabsch & Sander: H = N + unit(C_prev - O_prev), i.e. N-H points
        // opposite to the previous residue's C=O, 1.0 A from N.
        let mut prev = BackboneAtoms::new("A".to_string(), 1, "ALA".to_string(), None);
        prev.c = Some([0.0, 0.0, 0.0]);
        prev.o = Some([0.0, 1.23, 0.0]);
        let mut cur = BackboneAtoms::new("A".to_string(), 2, "GLY".to_string(), None);
        cur.n = Some([1.33, 0.0, 0.0]);

        let h = compute_virtual_hydrogen(&cur, Some(&prev)).unwrap();
        let expected = [1.33, -1.0, 0.0];
        for k in 0..3 {
            assert!((h[k] - expected[k]).abs() < 1e-9, "H = {h:?}");
        }
    }

    #[test]
    fn test_virtual_hydrogen_needs_previous_carbonyl_oxygen() {
        let mut prev = BackboneAtoms::new("A".to_string(), 1, "ALA".to_string(), None);
        prev.c = Some([0.0, 0.0, 0.0]);
        let mut cur = BackboneAtoms::new("A".to_string(), 2, "GLY".to_string(), None);
        cur.n = Some([1.33, 0.0, 0.0]);
        assert!(compute_virtual_hydrogen(&cur, Some(&prev)).is_none());
    }

    #[test]
    fn test_hbond_energy_threshold() {
        // Just verify the constant is correct
        assert!((HBOND_THRESHOLD - (-0.5)).abs() < 1e-10);
    }
}
