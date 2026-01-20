//! Secondary structure pattern recognition.
//!
//! This module implements pattern recognition for:
//! - Helices (α-helix, 3₁₀-helix, π-helix)
//! - Beta-sheets (parallel and antiparallel)
//! - Turns (H-bonded turns)
//! - Bends (high backbone curvature)
//! - PPII/κ-helix (dihedral-based)

use super::dihedral::BackboneDihedrals;
use super::hbond::BackboneAtoms;
use super::types::SecondaryStructure;
use std::collections::HashMap;

/// Minimum helix length for α-helix assignment.
pub const MIN_ALPHA_HELIX_LENGTH: usize = 4;

/// Minimum helix length for 3₁₀-helix assignment.
pub const MIN_310_HELIX_LENGTH: usize = 3;

/// Minimum helix length for π-helix assignment.
pub const MIN_PI_HELIX_LENGTH: usize = 5;

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

    /// Returns the minimum length for this helix type.
    pub fn min_length(&self) -> usize {
        match self {
            HelixType::Helix310 => MIN_310_HELIX_LENGTH,
            HelixType::Alpha => MIN_ALPHA_HELIX_LENGTH,
            HelixType::Pi => MIN_PI_HELIX_LENGTH,
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

/// Represents a helix segment.
#[derive(Debug, Clone)]
pub struct HelixSegment {
    /// Start residue index
    pub start: usize,
    /// End residue index (inclusive)
    pub end: usize,
    /// Helix type
    #[allow(dead_code)]
    pub helix_type: HelixType,
}

impl HelixSegment {
    /// Creates a new helix segment.
    pub fn new(start: usize, end: usize, helix_type: HelixType) -> Self {
        Self {
            start,
            end,
            helix_type,
        }
    }

    /// Returns the length of the helix segment.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.end - self.start + 1
    }

    /// Returns true if the segment is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
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

/// Represents a beta-sheet ladder (connected bridges).
#[derive(Debug, Clone)]
pub struct BetaLadder {
    /// All bridges in this ladder
    pub bridges: Vec<BetaBridge>,
    /// Ladder type (parallel or antiparallel)
    pub ladder_type: BetaType,
}

impl BetaLadder {
    /// Creates a new beta-ladder with one initial bridge.
    pub fn new(bridge: BetaBridge) -> Self {
        let ladder_type = bridge.bridge_type;
        Self {
            bridges: vec![bridge],
            ladder_type,
        }
    }

    /// Adds a bridge to the ladder.
    pub fn add_bridge(&mut self, bridge: BetaBridge) {
        self.bridges.push(bridge);
    }
}

/// Detects n-turns (H-bonded turns) of a specific length.
///
/// An n-turn at residue i exists if there's an H-bond from i+n to i.
#[allow(clippy::needless_range_loop)]
pub fn detect_n_turns(
    n: usize,
    hbond_matrix: &HashMap<(usize, usize), f64>,
    num_residues: usize,
) -> Vec<bool> {
    let mut turns = vec![false; num_residues];

    for i in 0..num_residues.saturating_sub(n) {
        // Check if there's an H-bond from i+n (donor) to i (acceptor)
        if hbond_matrix.contains_key(&(i + n, i)) {
            turns[i] = true;
        }
    }

    turns
}

/// Detects helix patterns based on consecutive n-turns.
///
/// A helix is assigned when there are consecutive n-turns.
/// For α-helix (n=4), we need consecutive 4-turns.
pub fn detect_helix_patterns(
    helix_type: HelixType,
    hbond_matrix: &HashMap<(usize, usize), f64>,
    num_residues: usize,
    residues: &[BackboneAtoms],
) -> Vec<HelixSegment> {
    let n = helix_type.offset();
    let min_len = helix_type.min_length();
    let mut segments = Vec::new();

    if num_residues < n + 1 {
        return segments;
    }

    // Detect n-turns
    let turns = detect_n_turns(n, hbond_matrix, num_residues);

    // Find consecutive turns that form helices
    let mut in_helix = false;
    let mut helix_start = 0;

    for i in 0..num_residues.saturating_sub(n - 1) {
        // Check for chain breaks
        let chain_ok = if i + 1 < num_residues {
            residues[i].chain_id == residues[i + 1].chain_id
        } else {
            true
        };

        if turns[i] && chain_ok {
            if !in_helix {
                helix_start = i;
                in_helix = true;
            }
        } else {
            if in_helix {
                // End of potential helix
                let helix_end = i + n - 2; // Extend to cover the last turn
                let segment_len = helix_end - helix_start + 1;
                if segment_len >= min_len && helix_end < num_residues {
                    segments.push(HelixSegment::new(helix_start, helix_end, helix_type));
                }
            }
            in_helix = false;
        }
    }

    // Handle helix at end of chain
    if in_helix {
        let helix_end = (num_residues - 1).min(
            turns
                .iter()
                .rposition(|&t| t)
                .map(|i| i + n - 1)
                .unwrap_or(helix_start + min_len - 1),
        );
        let segment_len = helix_end - helix_start + 1;
        if segment_len >= min_len {
            segments.push(HelixSegment::new(helix_start, helix_end, helix_type));
        }
    }

    segments
}

/// Detects beta-bridges between residue pairs.
///
/// A beta-bridge exists when there are H-bonds between residues i and j
/// in either parallel or antiparallel orientation.
///
/// Antiparallel: H-bonds i→j and j→i
/// Parallel: H-bonds i→j and j+1→i-1 (or i→j and j-1→i+1)
pub fn detect_beta_bridges(
    hbond_matrix: &HashMap<(usize, usize), f64>,
    num_residues: usize,
    _residues: &[BackboneAtoms],
) -> Vec<BetaBridge> {
    let mut bridges = Vec::new();

    for i in 0..num_residues {
        for j in (i + 3)..num_residues {
            // Check they're on the same chain or allow inter-chain
            // For simplicity, we currently only check intra-chain bridges

            // Antiparallel: i-1→j (NH of i-1 to CO of j) and j→i-1 (NH of j to CO of i-1)
            // Alternative: i→j-1 and j-1→i
            let antiparallel_1 =
                i > 0 && hbond_matrix.contains_key(&(i, j)) && hbond_matrix.contains_key(&(j, i));

            let antiparallel_2 = i > 0
                && j > 0
                && hbond_matrix.contains_key(&(i - 1, j))
                && hbond_matrix.contains_key(&(j, i - 1));

            let antiparallel_3 = j > 0
                && hbond_matrix.contains_key(&(i, j - 1))
                && hbond_matrix.contains_key(&(j - 1, i));

            if antiparallel_1 || antiparallel_2 || antiparallel_3 {
                bridges.push(BetaBridge::new(i, j, BetaType::Antiparallel));
                continue;
            }

            // Parallel: i-1→j (NH of i-1 to CO of j) and j→i+1 (NH of j to CO of i+1)
            // Alternative: j-1→i and i→j+1
            let parallel_1 = i > 0
                && i + 1 < num_residues
                && hbond_matrix.contains_key(&(i - 1, j))
                && hbond_matrix.contains_key(&(j, i + 1));

            let parallel_2 = j > 0
                && j + 1 < num_residues
                && hbond_matrix.contains_key(&(j - 1, i))
                && hbond_matrix.contains_key(&(i, j + 1));

            if parallel_1 || parallel_2 {
                bridges.push(BetaBridge::new(i, j, BetaType::Parallel));
            }
        }
    }

    bridges
}

/// Groups beta-bridges into ladders.
pub fn group_bridges_into_ladders(bridges: &[BetaBridge]) -> Vec<BetaLadder> {
    let mut ladders: Vec<BetaLadder> = Vec::new();

    for bridge in bridges {
        // Try to find an existing ladder this bridge belongs to
        let mut found = false;
        for ladder in &mut ladders {
            if ladder.ladder_type == bridge.bridge_type {
                // Check if this bridge is adjacent to any bridge in the ladder
                for existing in &ladder.bridges {
                    let adj1 = (existing.residue1 as isize - bridge.residue1 as isize).abs() <= 1
                        && (existing.residue2 as isize - bridge.residue2 as isize).abs() <= 1;
                    if adj1 {
                        ladder.add_bridge(bridge.clone());
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }

        if !found {
            ladders.push(BetaLadder::new(bridge.clone()));
        }
    }

    ladders
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
    fn test_helix_segment_length() {
        let segment = HelixSegment::new(0, 5, HelixType::Alpha);
        assert_eq!(segment.len(), 6);
        assert!(!segment.is_empty());
    }

    #[test]
    fn test_n_turns_detection() {
        let mut hbond_matrix = HashMap::new();
        // Add H-bonds for a 4-turn at position 0 (donor at 4, acceptor at 0)
        hbond_matrix.insert((4, 0), -1.0);
        hbond_matrix.insert((5, 1), -1.0);

        let turns = detect_n_turns(4, &hbond_matrix, 10);
        assert!(turns[0]);
        assert!(turns[1]);
        assert!(!turns[2]);
    }

    #[test]
    fn test_beta_bridge_creation() {
        let bridge = BetaBridge::new(5, 15, BetaType::Antiparallel);
        assert_eq!(bridge.residue1, 5);
        assert_eq!(bridge.residue2, 15);
        assert_eq!(bridge.bridge_type, BetaType::Antiparallel);
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
