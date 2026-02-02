//! Volume overlap calculation for protein-ligand complexes.
//!
//! This module implements grid-based volume overlap calculation following
//! the PoseBusters methodology. Volume overlap is defined as the percentage
//! of ligand volume that intersects with protein volume.

use super::clash::get_volume_radius;
use super::{VOLUME_GRID_SPACING, VOLUME_VDW_SCALING};
use crate::records::Atom;

/// Calculate the volume overlap between ligand and protein atoms.
///
/// Uses a grid-based approach where:
/// 1. A 3D grid is created around the ligand atoms
/// 2. Grid points inside ligand atom spheres are counted
/// 3. Grid points also inside protein atom spheres are counted
/// 4. Overlap percentage = shared_points / ligand_points × 100
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand
/// * `protein_atoms` - Atoms belonging to the protein (including nearby atoms)
///
/// # Returns
///
/// The volume overlap percentage (0.0 to 100.0).
/// Returns 0.0 if the ligand has no atoms or volume.
pub fn calculate_volume_overlap(ligand_atoms: &[&Atom], protein_atoms: &[&Atom]) -> f64 {
    if ligand_atoms.is_empty() {
        return 0.0;
    }

    // Get bounding box for the ligand with padding for atom radii
    let (min_coords, max_coords) = get_padded_bounding_box(ligand_atoms);

    // Count grid points
    let mut ligand_points = 0u64;
    let mut overlap_points = 0u64;

    // Iterate over grid
    let mut x = min_coords.0;
    while x <= max_coords.0 {
        let mut y = min_coords.1;
        while y <= max_coords.1 {
            let mut z = min_coords.2;
            while z <= max_coords.2 {
                let point = (x, y, z);

                // Check if point is inside any ligand atom
                if is_inside_any_atom(point, ligand_atoms) {
                    ligand_points += 1;

                    // Check if point is also inside any protein atom
                    if is_inside_any_atom(point, protein_atoms) {
                        overlap_points += 1;
                    }
                }

                z += VOLUME_GRID_SPACING;
            }
            y += VOLUME_GRID_SPACING;
        }
        x += VOLUME_GRID_SPACING;
    }

    if ligand_points == 0 {
        return 0.0;
    }

    (overlap_points as f64 / ligand_points as f64) * 100.0
}

/// Calculate an approximate ligand volume using grid counting.
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand
///
/// # Returns
///
/// The approximate volume in cubic Angstroms.
#[allow(dead_code)]
pub fn calculate_ligand_volume(ligand_atoms: &[&Atom]) -> f64 {
    if ligand_atoms.is_empty() {
        return 0.0;
    }

    let (min_coords, max_coords) = get_padded_bounding_box(ligand_atoms);

    let mut points = 0u64;

    let mut x = min_coords.0;
    while x <= max_coords.0 {
        let mut y = min_coords.1;
        while y <= max_coords.1 {
            let mut z = min_coords.2;
            while z <= max_coords.2 {
                if is_inside_any_atom((x, y, z), ligand_atoms) {
                    points += 1;
                }
                z += VOLUME_GRID_SPACING;
            }
            y += VOLUME_GRID_SPACING;
        }
        x += VOLUME_GRID_SPACING;
    }

    // Volume = number of points × volume per grid cell
    let voxel_volume = VOLUME_GRID_SPACING.powi(3);
    points as f64 * voxel_volume
}

/// Get the bounding box for a set of atoms, padded by the maximum atom radius.
///
/// # Arguments
///
/// * `atoms` - The atoms to compute the bounding box for
///
/// # Returns
///
/// A tuple of (min_coords, max_coords) where each is (x, y, z).
fn get_padded_bounding_box(atoms: &[&Atom]) -> ((f64, f64, f64), (f64, f64, f64)) {
    if atoms.is_empty() {
        return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0));
    }

    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    let mut z_min = f64::INFINITY;
    let mut z_max = f64::NEG_INFINITY;
    let mut max_radius = 0.0_f64;

    for atom in atoms {
        let radius = get_volume_radius(&atom.element, VOLUME_VDW_SCALING);
        max_radius = max_radius.max(radius);

        x_min = x_min.min(atom.x);
        x_max = x_max.max(atom.x);
        y_min = y_min.min(atom.y);
        y_max = y_max.max(atom.y);
        z_min = z_min.min(atom.z);
        z_max = z_max.max(atom.z);
    }

    // Pad by max radius + small buffer
    let padding = max_radius + VOLUME_GRID_SPACING;

    (
        (x_min - padding, y_min - padding, z_min - padding),
        (x_max + padding, y_max + padding, z_max + padding),
    )
}

/// Check if a point is inside any atom's sphere.
///
/// Uses scaled van der Waals radii (0.8 × vdW as per PoseBusters).
#[inline]
fn is_inside_any_atom(point: (f64, f64, f64), atoms: &[&Atom]) -> bool {
    for atom in atoms {
        let radius = get_volume_radius(&atom.element, VOLUME_VDW_SCALING);
        let radius_sq = radius * radius;

        let dx = point.0 - atom.x;
        let dy = point.1 - atom.y;
        let dz = point.2 - atom.z;
        let dist_sq = dx * dx + dy * dy + dz * dz;

        if dist_sq <= radius_sq {
            return true;
        }
    }
    false
}

/// Filter protein atoms to only those near the ligand for efficient overlap calculation.
///
/// # Arguments
///
/// * `ligand_atoms` - Atoms belonging to the ligand
/// * `protein_atoms` - All protein atoms
/// * `cutoff` - Distance cutoff for filtering
///
/// # Returns
///
/// Indices of protein atoms within cutoff distance of any ligand atom.
pub fn filter_nearby_protein_atoms<'a>(
    ligand_atoms: &[&Atom],
    protein_atoms: &[&'a Atom],
    cutoff: f64,
) -> Vec<&'a Atom> {
    let cutoff_sq = cutoff * cutoff;

    protein_atoms
        .iter()
        .filter(|prot_atom| {
            ligand_atoms.iter().any(|lig_atom| {
                let dx = prot_atom.x - lig_atom.x;
                let dy = prot_atom.y - lig_atom.y;
                let dz = prot_atom.z - lig_atom.z;
                dx * dx + dy * dy + dz * dz <= cutoff_sq
            })
        })
        .copied()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_atom(serial: i32, x: f64, y: f64, z: f64, element: &str) -> Atom {
        Atom {
            serial,
            name: "C".to_string(),
            alt_loc: None,
            residue_name: "LIG".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            ins_code: None,
            is_hetatm: true,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: element.to_string(),
        }
    }

    #[test]
    fn test_no_overlap() {
        // Ligand and protein far apart
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let prot_atom = create_test_atom(2, 20.0, 0.0, 0.0, "C");

        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let prot_atoms: Vec<&Atom> = vec![&prot_atom];

        let overlap = calculate_volume_overlap(&lig_atoms, &prot_atoms);

        assert!(overlap < 0.1, "Expected no overlap, got {}%", overlap);
    }

    #[test]
    fn test_complete_overlap() {
        // Ligand and protein at same position
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let prot_atom = create_test_atom(2, 0.0, 0.0, 0.0, "C");

        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let prot_atoms: Vec<&Atom> = vec![&prot_atom];

        let overlap = calculate_volume_overlap(&lig_atoms, &prot_atoms);

        // Should be close to 100% (exact value depends on grid resolution)
        assert!(overlap > 90.0, "Expected high overlap, got {}%", overlap);
    }

    #[test]
    fn test_partial_overlap() {
        // Ligand and protein partially overlapping
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let prot_atom = create_test_atom(2, 1.5, 0.0, 0.0, "C"); // Partial overlap

        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let prot_atoms: Vec<&Atom> = vec![&prot_atom];

        let overlap = calculate_volume_overlap(&lig_atoms, &prot_atoms);

        // Should be between 0 and 100
        assert!(
            overlap > 0.0 && overlap < 100.0,
            "Expected partial overlap, got {}%",
            overlap
        );
    }

    #[test]
    fn test_empty_ligand() {
        let lig_atoms: Vec<&Atom> = vec![];
        let prot_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let prot_atoms: Vec<&Atom> = vec![&prot_atom];

        let overlap = calculate_volume_overlap(&lig_atoms, &prot_atoms);

        assert_eq!(overlap, 0.0);
    }

    #[test]
    fn test_empty_protein() {
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let prot_atoms: Vec<&Atom> = vec![];

        let overlap = calculate_volume_overlap(&lig_atoms, &prot_atoms);

        assert_eq!(overlap, 0.0);
    }

    #[test]
    fn test_ligand_volume_single_atom() {
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let lig_atoms: Vec<&Atom> = vec![&lig_atom];

        let volume = calculate_ligand_volume(&lig_atoms);

        // Carbon vdW radius = 1.70, scaled = 1.36
        // Volume of sphere = 4/3 * pi * r^3 ≈ 10.5 Å³
        // Grid approximation will be close but not exact
        assert!(
            volume > 5.0 && volume < 20.0,
            "Unexpected volume: {}",
            volume
        );
    }

    #[test]
    fn test_filter_nearby_atoms() {
        let lig_atom = create_test_atom(1, 0.0, 0.0, 0.0, "C");
        let prot_near = create_test_atom(2, 3.0, 0.0, 0.0, "C");
        let prot_far = create_test_atom(3, 20.0, 0.0, 0.0, "C");

        let lig_atoms: Vec<&Atom> = vec![&lig_atom];
        let prot_atoms = vec![&prot_near, &prot_far];

        let nearby = filter_nearby_protein_atoms(&lig_atoms, &prot_atoms, 5.0);

        assert_eq!(nearby.len(), 1);
        assert_eq!(nearby[0].serial, 2);
    }

    #[test]
    fn test_bounding_box() {
        let atom1 = create_test_atom(1, -5.0, -3.0, -1.0, "C");
        let atom2 = create_test_atom(2, 5.0, 3.0, 1.0, "C");

        let atoms: Vec<&Atom> = vec![&atom1, &atom2];
        let (min, max) = get_padded_bounding_box(&atoms);

        // Should include padding for atom radii
        assert!(min.0 < -5.0);
        assert!(max.0 > 5.0);
        assert!(min.1 < -3.0);
        assert!(max.1 > 3.0);
    }
}
