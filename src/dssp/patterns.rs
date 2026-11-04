//! Secondary structure pattern recognition.
//!
//! This module implements pattern recognition for:
//! - Helices (α-helix, 3₁₀-helix, π-helix)
//! - Beta-sheets (parallel and antiparallel)
//! - Turns (H-bonded turns)
//! - Bends (high backbone curvature)
//! - PPII/κ-helix (dihedral-based)

use super::dihedral::BackboneDihedrals;
use super::hbond::{BackboneAtoms, no_chain_break, residues_linked};
use super::types::SecondaryStructure;
use std::collections::HashMap;

/// Threshold angle for bend detection (degrees).
pub const BEND_ANGLE_THRESHOLD: f64 = 70.0;

/// Helix pattern types based on H-bond offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HelixType {
    /// 3₁₀-helix: i → i+3 H-bond
    Helix310,
    /// α-helix: i → i+4 H-bond
    Alpha,
    /// π-helix: i → i+5 H-bond
    Pi,
}

impl HelixType {
    /// Returns the H-bond offset for this helix type.
    pub fn offset(&self) -> usize {
        match self {
            HelixType::Helix310 => 3,
            HelixType::Alpha => 4,
            HelixType::Pi => 5,
        }
    }

    /// Returns the corresponding SecondaryStructure enum.
    pub fn to_secondary_structure(self) -> SecondaryStructure {
        match self {
            HelixType::Helix310 => SecondaryStructure::Helix310,
            HelixType::Alpha => SecondaryStructure::AlphaHelix,
            HelixType::Pi => SecondaryStructure::PiHelix,
        }
    }
}

/// Beta-sheet pattern types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BetaType {
    /// Parallel beta-sheet
    Parallel,
    /// Antiparallel beta-sheet
    Antiparallel,
}

/// Represents a beta-bridge between two residues.
#[derive(Debug, Clone)]
pub struct BetaBridge {
    /// First residue index
    pub residue1: usize,
    /// Second residue index
    pub residue2: usize,
    /// Bridge type (parallel or antiparallel)
    pub bridge_type: BetaType,
}

impl BetaBridge {
    /// Creates a new beta-bridge.
    pub fn new(residue1: usize, residue2: usize, bridge_type: BetaType) -> Self {
        Self {
            residue1,
            residue2,
            bridge_type,
        }
    }
}

/// A beta ladder: consecutive bridges of one type, possibly joined across
/// beta-bulges.
#[derive(Debug, Clone)]
pub struct BetaLadder {
    /// Ladder type (parallel or antiparallel)
    pub ladder_type: BetaType,
    /// Residue indices on the first strand, ascending
    pub i: Vec<usize>,
    /// Residue indices on the partner strand, ascending
    pub j: Vec<usize>,
}

impl BetaLadder {
    /// Creates a ladder from a single bridge.
    pub fn new(bridge: &BetaBridge) -> Self {
        Self {
            ladder_type: bridge.bridge_type,
            i: vec![bridge.residue1],
            j: vec![bridge.residue2],
        }
    }

    /// First and last residue of the first strand.
    pub fn i_span(&self) -> (usize, usize) {
        (self.i[0], self.i[self.i.len() - 1])
    }

    /// First and last residue of the partner strand.
    pub fn j_span(&self) -> (usize, usize) {
        (self.j[0], self.j[self.j.len() - 1])
    }

    /// Number of bridges in the ladder.
    pub fn bridge_count(&self) -> usize {
        self.i.len()
    }
}

/// Marks the residues that start an n-turn, as mkdssp's helix "start" flags.
///
/// `starts[i]` is true when the N-H of residue i+n bonds to the C=O of residue
/// i and residues i..=i+n contain no chain break.
pub fn detect_helix_starts(
    n: usize,
    hbond_matrix: &HashMap<(usize, usize), f64>,
    residues: &[BackboneAtoms],
) -> Vec<bool> {
    (0..residues.len())
        .map(|i| {
            i + n < residues.len()
                && hbond_matrix.contains_key(&(i + n, i))
                && no_chain_break(residues, i, i + n)
        })
        .collect()
}

/// Detects beta-bridges between residue pairs.
///
/// Uses the Kabsch & Sander definitions. Writing Hbond(x, y) for an H-bond from
/// the C=O of residue x to the N-H of residue y (matrix key `(y, x)`), a bridge
/// between residues i and j (j ≥ i + 3) is
///
/// - antiparallel if [Hbond(i, j) and Hbond(j, i)] or
///   [Hbond(i-1, j+1) and Hbond(j-1, i+1)];
/// - parallel if [Hbond(i-1, j) and Hbond(j, i+1)] or
///   [Hbond(j-1, i) and Hbond(i, j+1)].
///
/// Both residues need intact neighbours: i-1, i+1, j-1 and j+1 must exist and
/// be peptide-bonded to i and j (see [`residues_linked`]).
pub fn detect_beta_bridges(
    hbond_matrix: &HashMap<(usize, usize), f64>,
    num_residues: usize,
    residues: &[BackboneAtoms],
) -> Vec<BetaBridge> {
    let mut bridges = Vec::new();
    let num_residues = num_residues.min(residues.len());

    // C=O of `acceptor` bonds to N-H of `donor`.
    let hbond = |acceptor: usize, donor: usize| hbond_matrix.contains_key(&(donor, acceptor));

    // Residue i has neighbours i-1 and i+1 with no chain break between them.
    let has_neighbours = |i: usize| {
        i > 0
            && i + 1 < num_residues
            && residues_linked(residues, i - 1)
            && residues_linked(residues, i)
    };

    for i in 0..num_residues {
        if !has_neighbours(i) {
            continue;
        }
        for j in (i + 3)..num_residues {
            if !has_neighbours(j) {
                continue;
            }

            let antiparallel =
                (hbond(i, j) && hbond(j, i)) || (hbond(i - 1, j + 1) && hbond(j - 1, i + 1));
            if antiparallel {
                bridges.push(BetaBridge::new(i, j, BetaType::Antiparallel));
                continue;
            }

            let parallel =
                (hbond(i - 1, j) && hbond(j, i + 1)) || (hbond(j - 1, i) && hbond(i, j + 1));
            if parallel {
                bridges.push(BetaBridge::new(i, j, BetaType::Parallel));
            }
        }
    }

    bridges
}

/// Groups beta-bridges into ladders and joins ladders across beta-bulges,
/// following mkdssp.
///
/// A bridge (i, j) extends a ladder whose last bridge is (i-1, j-1) for
/// parallel or (i-1, j+1) for antiparallel ladders. Two ladders of the same
/// type are then joined when the gap between them is at most 1 residue on one
/// strand and at most 4 on the other, with no chain break across either gap.
/// `bridges` must be ordered by first residue, as [`detect_beta_bridges`]
/// returns them.
pub fn group_bridges_into_ladders(
    bridges: &[BetaBridge],
    residues: &[BackboneAtoms],
) -> Vec<BetaLadder> {
    let mut ladders: Vec<BetaLadder> = Vec::new();

    for bridge in bridges {
        let (i, j) = (bridge.residue1, bridge.residue2);
        let extends = |l: &BetaLadder| {
            let (_, last_i) = l.i_span();
            let (first_j, last_j) = l.j_span();
            l.ladder_type == bridge.bridge_type
                && last_i + 1 == i
                && match bridge.bridge_type {
                    BetaType::Parallel => last_j + 1 == j,
                    BetaType::Antiparallel => first_j == j + 1,
                }
        };
        match ladders.iter_mut().find(|l| extends(l)) {
            Some(ladder) => {
                ladder.i.push(i);
                match bridge.bridge_type {
                    BetaType::Parallel => ladder.j.push(j),
                    BetaType::Antiparallel => ladder.j.insert(0, j),
                }
            }
            None => ladders.push(BetaLadder::new(bridge)),
        }
    }

    // Join ladders separated by a beta-bulge.
    let mut a = 0;
    while a < ladders.len() {
        let mut b = a + 1;
        while b < ladders.len() {
            if is_bulge(&ladders[a], &ladders[b], residues) {
                let other = ladders.remove(b);
                let ladder = &mut ladders[a];
                ladder.i.extend(other.i);
                match ladder.ladder_type {
                    BetaType::Parallel => ladder.j.extend(other.j),
                    BetaType::Antiparallel => {
                        let mut j = other.j;
                        j.append(&mut ladder.j);
                        ladder.j = j;
                    }
                }
            } else {
                b += 1;
            }
        }
        a += 1;
    }

    ladders
}

/// Whether ladder `second` follows ladder `first` across a beta-bulge
/// (mkdssp's criterion).
fn is_bulge(first: &BetaLadder, second: &BetaLadder, residues: &[BackboneAtoms]) -> bool {
    if first.ladder_type != second.ladder_type {
        return false;
    }
    let (ibi, iei) = first.i_span();
    let (jbi, jei) = first.j_span();
    let (ibj, iej) = second.i_span();
    let (jbj, jej) = second.j_span();

    if !no_chain_break(residues, ibi.min(ibj), iei.max(iej))
        || !no_chain_break(residues, jbi.min(jbj), jei.max(jej))
    {
        return false;
    }
    // The second ladder must start after the first ends on the first strand.
    if ibj <= iei || ibj - iei >= 6 {
        return false;
    }
    let i_gap = ibj - iei;
    // Gap on the partner strand; it must also run in the ladder's direction.
    let j_gap = match first.ladder_type {
        BetaType::Parallel if jbj > jei => jbj - jei,
        BetaType::Antiparallel if jbi > jej => jbi - jej,
        _ => return false,
    };
    (j_gap < 6 && i_gap < 3) || j_gap < 3
}

/// Detects bends based on backbone curvature.
///
/// A bend exists at residue i if the angle between
/// CA(i-2)-CA(i) and CA(i)-CA(i+2) is greater than 70°.
pub fn detect_bends(residues: &[BackboneAtoms]) -> Vec<bool> {
    let n = residues.len();
    let mut bends = vec![false; n];

    if n < 5 {
        return bends;
    }

    for i in 2..(n - 2) {
        // Check chain continuity
        if residues[i - 2].chain_id != residues[i].chain_id
            || residues[i].chain_id != residues[i + 2].chain_id
        {
            continue;
        }

        // Get CA coordinates
        let ca_im2 = match residues[i - 2].ca {
            Some(ca) => ca,
            None => continue,
        };
        let ca_i = match residues[i].ca {
            Some(ca) => ca,
            None => continue,
        };
        let ca_ip2 = match residues[i + 2].ca {
            Some(ca) => ca,
            None => continue,
        };

        // Calculate vectors
        let v1 = [
            ca_i[0] - ca_im2[0],
            ca_i[1] - ca_im2[1],
            ca_i[2] - ca_im2[2],
        ];
        let v2 = [
            ca_ip2[0] - ca_i[0],
            ca_ip2[1] - ca_i[1],
            ca_ip2[2] - ca_i[2],
        ];

        // Calculate angle
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();

        if len1 < 1e-10 || len2 < 1e-10 {
            continue;
        }

        let cos_angle = dot / (len1 * len2);
        let angle_deg = cos_angle.clamp(-1.0, 1.0).acos() * 180.0 / std::f64::consts::PI;

        if angle_deg > BEND_ANGLE_THRESHOLD {
            bends[i] = true;
        }
    }

    bends
}

/// Detects PPII/κ-helix regions based on backbone dihedrals.
///
/// PPII helix is assigned to coil residues with:
/// - φ ≈ -75° ± 29°
/// - ψ ≈ +145° ± 29°
/// - At least 2 consecutive residues meeting these criteria
///
/// Returns indices of residues that should be assigned to PPII.
#[allow(clippy::needless_range_loop)]
pub fn detect_ppii_regions(
    dihedrals: &[BackboneDihedrals],
    current_assignments: &[SecondaryStructure],
) -> Vec<bool> {
    let n = dihedrals.len();
    let mut ppii = vec![false; n];

    if n < super::dihedral::PPII_MIN_CONSECUTIVE {
        return ppii;
    }

    // Find consecutive PPII-like residues that are currently coil
    let ppii_like: Vec<bool> = dihedrals
        .iter()
        .zip(current_assignments.iter())
        .map(|(d, &ss)| d.is_ppii_like() && ss == SecondaryStructure::Coil)
        .collect();

    // Find stretches of consecutive PPII-like residues
    let mut start = 0;
    let mut in_stretch = false;

    for i in 0..n {
        if ppii_like[i] {
            if !in_stretch {
                start = i;
                in_stretch = true;
            }
        } else {
            if in_stretch {
                let stretch_len = i - start;
                if stretch_len >= super::dihedral::PPII_MIN_CONSECUTIVE {
                    for j in start..i {
                        ppii[j] = true;
                    }
                }
            }
            in_stretch = false;
        }
    }

    // Handle stretch at end
    if in_stretch {
        let stretch_len = n - start;
        if stretch_len >= super::dihedral::PPII_MIN_CONSECUTIVE {
            for j in start..n {
                ppii[j] = true;
            }
        }
    }

    ppii
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_helix_type_offset() {
        assert_eq!(HelixType::Helix310.offset(), 3);
        assert_eq!(HelixType::Alpha.offset(), 4);
        assert_eq!(HelixType::Pi.offset(), 5);
    }

    #[test]
    fn test_helix_starts_detection() {
        let mut hbond_matrix = HashMap::new();
        // 4-turns at 0 and 1: N-H of 4 -> C=O of 0, N-H of 5 -> C=O of 1
        hbond_matrix.insert((4, 0), -1.0);
        hbond_matrix.insert((5, 1), -1.0);

        let starts = detect_helix_starts(4, &hbond_matrix, &chain(10));
        assert_eq!(&starts[..3], &[true, true, false]);

        // A chain break inside the turn removes it.
        let mut broken = chain(10);
        for r in broken.iter_mut().skip(3) {
            r.residue_seq += 10;
        }
        let starts = detect_helix_starts(4, &hbond_matrix, &broken);
        assert!(!starts[0] && !starts[1]);
    }

    #[test]
    fn test_beta_bridge_creation() {
        let bridge = BetaBridge::new(5, 15, BetaType::Antiparallel);
        assert_eq!(bridge.residue1, 5);
        assert_eq!(bridge.residue2, 15);
        assert_eq!(bridge.bridge_type, BetaType::Antiparallel);
    }

    /// `n` consecutive residues on chain A, numbered from 1.
    fn chain(n: usize) -> Vec<BackboneAtoms> {
        (0..n)
            .map(|i| BackboneAtoms::new("A".to_string(), i as i32 + 1, "ALA".to_string(), None))
            .collect()
    }

    /// Bridges found for H-bonds given as (donor, acceptor) residue indices.
    fn bridges_for(
        hbonds: &[(usize, usize)],
        residues: &[BackboneAtoms],
    ) -> Vec<(usize, usize, BetaType)> {
        let matrix: HashMap<(usize, usize), f64> = hbonds.iter().map(|&k| (k, -1.0)).collect();
        detect_beta_bridges(&matrix, residues.len(), residues)
            .into_iter()
            .map(|b| (b.residue1, b.residue2, b.bridge_type))
            .collect()
    }

    // Kabsch & Sander bridge patterns, with Hbond(x, y) meaning C=O of x
    // accepts from N-H of y, i.e. matrix key (donor y, acceptor x).

    #[test]
    fn test_antiparallel_bridge_mutual_hbonds() {
        // Hbond(i, j) and Hbond(j, i) for i = 2, j = 10
        let found = bridges_for(&[(10, 2), (2, 10)], &chain(14));
        assert_eq!(found, vec![(2, 10, BetaType::Antiparallel)]);
    }

    #[test]
    fn test_no_bridge_at_chain_ends() {
        // Like mkdssp, a bridge residue needs both sequence neighbours.
        assert!(bridges_for(&[(10, 0), (0, 10)], &chain(14)).is_empty());
        assert!(bridges_for(&[(13, 2), (2, 13)], &chain(14)).is_empty());
    }

    #[test]
    fn test_antiparallel_bridge_flanking_hbonds() {
        // Hbond(i-1, j+1) and Hbond(j-1, i+1) for i = 2, j = 10
        let found = bridges_for(&[(11, 1), (3, 9)], &chain(14));
        assert_eq!(found, vec![(2, 10, BetaType::Antiparallel)]);
    }

    #[test]
    fn test_parallel_bridge_first_pattern() {
        // Hbond(i-1, j) and Hbond(j, i+1) for i = 2, j = 10
        let found = bridges_for(&[(10, 1), (3, 10)], &chain(14));
        assert_eq!(found, vec![(2, 10, BetaType::Parallel)]);
    }

    #[test]
    fn test_parallel_bridge_second_pattern() {
        // Hbond(j-1, i) and Hbond(i, j+1) for i = 2, j = 10
        let found = bridges_for(&[(2, 9), (11, 2)], &chain(14));
        assert_eq!(found, vec![(2, 10, BetaType::Parallel)]);
    }

    #[test]
    fn test_parallel_pattern_with_reversed_hbonds_is_not_a_bridge() {
        // Hbond(j, i-1) and Hbond(i+1, j): the parallel pattern with every
        // H-bond pointing the wrong way.
        assert!(bridges_for(&[(1, 10), (10, 3)], &chain(14)).is_empty());
        assert!(bridges_for(&[(9, 2), (2, 11)], &chain(14)).is_empty());
    }

    #[test]
    fn test_bridge_requires_no_chain_break_around_residues() {
        let mut residues = chain(14);
        // Gap in numbering between index 1 and index 2.
        for r in residues.iter_mut().skip(2) {
            r.residue_seq += 5;
        }
        assert!(bridges_for(&[(11, 1), (3, 9)], &residues).is_empty());
    }

    fn ladders_for(bridges: &[(usize, usize)], ty: BetaType) -> Vec<(Vec<usize>, Vec<usize>)> {
        let bridges: Vec<_> = bridges
            .iter()
            .map(|&(i, j)| BetaBridge::new(i, j, ty))
            .collect();
        group_bridges_into_ladders(&bridges, &chain(40))
            .into_iter()
            .map(|l| (l.i, l.j))
            .collect()
    }

    #[test]
    fn test_ladders_from_consecutive_bridges() {
        let anti = ladders_for(&[(2, 20), (3, 19), (4, 18)], BetaType::Antiparallel);
        assert_eq!(anti, vec![(vec![2, 3, 4], vec![18, 19, 20])]);
        let par = ladders_for(&[(2, 20), (3, 21), (4, 22)], BetaType::Parallel);
        assert_eq!(par, vec![(vec![2, 3, 4], vec![20, 21, 22])]);
    }

    #[test]
    fn test_ladders_joined_across_beta_bulge() {
        // One extra residue (4) on the first strand, none on the partner.
        let found = ladders_for(
            &[(2, 20), (3, 19), (5, 18), (6, 17)],
            BetaType::Antiparallel,
        );
        assert_eq!(found, vec![(vec![2, 3, 5, 6], vec![17, 18, 19, 20])]);
    }

    #[test]
    fn test_ladders_not_joined_across_long_gap() {
        // Six-residue gap on the first strand is not a bulge.
        let found = ladders_for(
            &[(2, 30), (3, 29), (10, 28), (11, 27)],
            BetaType::Antiparallel,
        );
        assert_eq!(found.len(), 2);
    }

    #[test]
    fn test_bend_detection_small() {
        // Structure too small for bend detection
        let residues: Vec<BackboneAtoms> = vec![];
        let bends = detect_bends(&residues);
        assert!(bends.is_empty());
    }

    #[test]
    fn test_ppii_detection() {
        let mut dihedrals = vec![BackboneDihedrals::default(); 5];
        let assignments = vec![SecondaryStructure::Coil; 5];

        // No PPII angles set
        let ppii = detect_ppii_regions(&dihedrals, &assignments);
        assert!(ppii.iter().all(|&p| !p));

        // Set PPII-like angles for residues 1-3
        for d in dihedrals.iter_mut().skip(1).take(3) {
            d.phi = Some(-75.0);
            d.psi = Some(145.0);
        }

        let ppii = detect_ppii_regions(&dihedrals, &assignments);
        // Residues 1-3 should be PPII (3 consecutive)
        assert!(!ppii[0]);
        assert!(ppii[1]);
        assert!(ppii[2]);
        assert!(ppii[3]);
        assert!(!ppii[4]);
    }
}
