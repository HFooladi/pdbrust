//! Backbone dihedral angle calculation for PPII helix detection.
//!
//! This module calculates the φ (phi) and ψ (psi) backbone dihedral angles
//! needed for polyproline II (PPII) helix detection in DSSP 4.
//!
//! PPII helix detection criteria (DSSP 4):
//! - φ = -75° ± 29°
//! - ψ = +145° ± 29°
//! - At least 2 consecutive residues matching these criteria
//! - Residues must be in coil conformation (not assigned to other SS types)

use super::hbond::BackboneAtoms;
use std::f64::consts::PI;

/// PPII φ angle center (degrees).
pub const PPII_PHI_CENTER: f64 = -75.0;

/// PPII ψ angle center (degrees).
pub const PPII_PSI_CENTER: f64 = 145.0;

/// PPII angle tolerance (degrees).
pub const PPII_ANGLE_TOLERANCE: f64 = 29.0;

/// Minimum consecutive residues for PPII helix assignment.
pub const PPII_MIN_CONSECUTIVE: usize = 2;

/// Backbone dihedral angles for a single residue.
#[derive(Debug, Clone, Default)]
pub struct BackboneDihedrals {
    /// φ (phi) angle: C(i-1)-N(i)-CA(i)-C(i) in degrees
    pub phi: Option<f64>,
    /// ψ (psi) angle: N(i)-CA(i)-C(i)-N(i+1) in degrees
    pub psi: Option<f64>,
    /// ω (omega) angle: CA(i-1)-C(i-1)-N(i)-CA(i) in degrees
    pub omega: Option<f64>,
}

impl BackboneDihedrals {
    /// Returns true if this residue has PPII-like dihedral angles.
    pub fn is_ppii_like(&self) -> bool {
        match (self.phi, self.psi) {
            (Some(phi), Some(psi)) => {
                let phi_diff = angle_difference(phi, PPII_PHI_CENTER);
                let psi_diff = angle_difference(psi, PPII_PSI_CENTER);
                phi_diff <= PPII_ANGLE_TOLERANCE && psi_diff <= PPII_ANGLE_TOLERANCE
            }
            _ => false,
        }
    }
}

/// Calculates the absolute difference between two angles in degrees,
/// accounting for periodicity (-180 to 180).
fn angle_difference(a1: f64, a2: f64) -> f64 {
    let mut diff = (a1 - a2).abs();
    if diff > 180.0 {
        diff = 360.0 - diff;
    }
    diff
}

/// Calculates the dihedral angle between four points.
///
/// The dihedral angle is the angle between the planes defined by
/// (p1, p2, p3) and (p2, p3, p4).
///
/// Returns the angle in degrees, in the range [-180, 180].
pub fn calculate_dihedral(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> f64 {
    // Vector b1 = p2 - p1
    let b1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    // Vector b2 = p3 - p2
    let b2 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];
    // Vector b3 = p4 - p3
    let b3 = [p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2]];

    // Normal to plane 1 (b1 x b2)
    let n1 = cross_product(&b1, &b2);
    // Normal to plane 2 (b2 x b3)
    let n2 = cross_product(&b2, &b3);

    // m1 = n1 x b2_normalized
    let b2_len = vector_length(&b2);
    if b2_len < 1e-10 {
        return 0.0;
    }
    let b2_normalized = [b2[0] / b2_len, b2[1] / b2_len, b2[2] / b2_len];
    let m1 = cross_product(&n1, &b2_normalized);

    // Calculate dihedral
    let x = dot_product(&n1, &n2);
    let y = dot_product(&m1, &n2);

    let angle_rad = y.atan2(x);
    angle_rad * 180.0 / PI
}

/// Calculates the cross product of two 3D vectors.
#[inline]
fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Calculates the dot product of two 3D vectors.
#[inline]
fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Calculates the length of a 3D vector.
#[inline]
fn vector_length(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Calculates φ (phi) angle for a residue.
///
/// φ = dihedral(C(i-1), N(i), CA(i), C(i))
pub fn calculate_phi(prev_residue: &BackboneAtoms, curr_residue: &BackboneAtoms) -> Option<f64> {
    let c_prev = prev_residue.c?;
    let n = curr_residue.n?;
    let ca = curr_residue.ca?;
    let c = curr_residue.c?;

    Some(calculate_dihedral(&c_prev, &n, &ca, &c))
}

/// Calculates ψ (psi) angle for a residue.
///
/// ψ = dihedral(N(i), CA(i), C(i), N(i+1))
pub fn calculate_psi(curr_residue: &BackboneAtoms, next_residue: &BackboneAtoms) -> Option<f64> {
    let n = curr_residue.n?;
    let ca = curr_residue.ca?;
    let c = curr_residue.c?;
    let n_next = next_residue.n?;

    Some(calculate_dihedral(&n, &ca, &c, &n_next))
}

/// Calculates ω (omega) angle for a peptide bond.
///
/// ω = dihedral(CA(i-1), C(i-1), N(i), CA(i))
pub fn calculate_omega(prev_residue: &BackboneAtoms, curr_residue: &BackboneAtoms) -> Option<f64> {
    let ca_prev = prev_residue.ca?;
    let c_prev = prev_residue.c?;
    let n = curr_residue.n?;
    let ca = curr_residue.ca?;

    Some(calculate_dihedral(&ca_prev, &c_prev, &n, &ca))
}

/// Calculates backbone dihedral angles for all residues.
///
/// Returns a vector of BackboneDihedrals, one per residue.
pub fn calculate_all_dihedrals(residues: &[BackboneAtoms]) -> Vec<BackboneDihedrals> {
    let n = residues.len();
    if n == 0 {
        return Vec::new();
    }

    let mut dihedrals = vec![BackboneDihedrals::default(); n];

    for i in 0..n {
        // Calculate φ (requires previous residue)
        if i > 0 && residues[i].chain_id == residues[i - 1].chain_id {
            // Check for sequence continuity
            let seq_diff = residues[i].residue_seq - residues[i - 1].residue_seq;
            if seq_diff == 1 || (seq_diff == 0 && residues[i].ins_code != residues[i - 1].ins_code)
            {
                dihedrals[i].phi = calculate_phi(&residues[i - 1], &residues[i]);
                dihedrals[i].omega = calculate_omega(&residues[i - 1], &residues[i]);
            }
        }

        // Calculate ψ (requires next residue)
        if i + 1 < n && residues[i].chain_id == residues[i + 1].chain_id {
            // Check for sequence continuity
            let seq_diff = residues[i + 1].residue_seq - residues[i].residue_seq;
            if seq_diff == 1 || (seq_diff == 0 && residues[i + 1].ins_code != residues[i].ins_code)
            {
                dihedrals[i].psi = calculate_psi(&residues[i], &residues[i + 1]);
            }
        }
    }

    dihedrals
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_angle_difference() {
        assert!((angle_difference(10.0, 20.0) - 10.0).abs() < 1e-10);
        assert!((angle_difference(-170.0, 170.0) - 20.0).abs() < 1e-10);
        assert!((angle_difference(0.0, 0.0) - 0.0).abs() < 1e-10);
        assert!((angle_difference(-180.0, 180.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross_product(&a, &b);
        assert!((c[0] - 0.0).abs() < 1e-10);
        assert!((c[1] - 0.0).abs() < 1e-10);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dot_product() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert!((dot_product(&a, &b) - 32.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_length() {
        let v = [3.0, 4.0, 0.0];
        assert!((vector_length(&v) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dihedral_calculation() {
        // Test with a simple trans configuration (should be ~180 or -180)
        let p1 = [0.0, 0.0, 0.0];
        let p2 = [1.0, 0.0, 0.0];
        let p3 = [2.0, 0.0, 0.0];
        let p4 = [3.0, 0.0, 0.0];

        // Collinear points, dihedral is undefined but should return 0
        let dihedral = calculate_dihedral(&p1, &p2, &p3, &p4);
        // For collinear points, result is undefined
        assert!(dihedral.is_finite());
    }

    #[test]
    fn test_dihedral_90_degrees() {
        // Set up points for a ~90 degree dihedral
        let p1 = [0.0, 1.0, 0.0];
        let p2 = [0.0, 0.0, 0.0];
        let p3 = [1.0, 0.0, 0.0];
        let p4 = [1.0, 0.0, 1.0];

        let dihedral = calculate_dihedral(&p1, &p2, &p3, &p4);
        // Should be approximately 90 degrees
        assert!((dihedral.abs() - 90.0).abs() < 1.0);
    }

    #[test]
    fn test_ppii_detection() {
        let mut dihedral = BackboneDihedrals::default();

        // Not PPII without phi/psi
        assert!(!dihedral.is_ppii_like());

        // Set to PPII-like values
        dihedral.phi = Some(-75.0);
        dihedral.psi = Some(145.0);
        assert!(dihedral.is_ppii_like());

        // Just outside tolerance
        dihedral.phi = Some(-75.0 + PPII_ANGLE_TOLERANCE + 5.0);
        assert!(!dihedral.is_ppii_like());

        // Within tolerance
        dihedral.phi = Some(-75.0 + PPII_ANGLE_TOLERANCE - 5.0);
        assert!(dihedral.is_ppii_like());
    }

    #[test]
    fn test_empty_dihedrals() {
        let residues: Vec<BackboneAtoms> = vec![];
        let dihedrals = calculate_all_dihedrals(&residues);
        assert!(dihedrals.is_empty());
    }
}
