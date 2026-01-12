//! Coordinate transformation utilities for structure superposition.
//!
//! This module provides functions for:
//! - Extracting coordinates from structures based on atom selection
//! - Computing centroids
//! - Applying rotation and translation transformations

use crate::core::PdbStructure;
use crate::error::PdbError;

/// A 3D coordinate point as (x, y, z).
pub type Point3D = (f64, f64, f64);

/// Residue identifier: (chain_id, residue_seq, residue_name).
pub type ResidueInfo = (String, i32, String);

/// Coordinate with associated residue information.
pub type CoordWithResidue = (ResidueInfo, Point3D);

/// Atom selection criteria for RMSD and alignment calculations.
#[derive(Debug, Clone, Default, PartialEq)]
pub enum AtomSelection {
    /// Use only CA (alpha-carbon) atoms. This is the default and most common
    /// choice for protein structure comparison.
    #[default]
    CaOnly,
    /// Use backbone atoms (N, CA, C, O).
    Backbone,
    /// Use all atoms. Requires exact atom correspondence between structures.
    AllAtoms,
    /// Custom selection based on atom name patterns.
    Custom(Vec<String>),
}

impl AtomSelection {
    /// Check if an atom name matches this selection.
    pub fn matches(&self, atom_name: &str) -> bool {
        let name = atom_name.trim();
        match self {
            AtomSelection::CaOnly => name == "CA",
            AtomSelection::Backbone => {
                matches!(name, "N" | "CA" | "C" | "O")
            }
            AtomSelection::AllAtoms => true,
            AtomSelection::Custom(names) => names.iter().any(|n| n.trim() == name),
        }
    }
}

/// Extract coordinates from a structure based on atom selection.
///
/// Returns coordinates as a vector of (x, y, z) tuples, along with
/// metadata about the atoms (residue info) for per-residue calculations.
///
/// # Arguments
/// * `structure` - The structure to extract coordinates from
/// * `selection` - The atom selection criteria
/// * `chain_id` - Optional chain ID to filter by
///
/// # Returns
/// Vector of (x, y, z) coordinates for matching atoms
pub fn extract_coords_by_selection(
    structure: &PdbStructure,
    selection: &AtomSelection,
    chain_id: Option<&str>,
) -> Vec<(f64, f64, f64)> {
    structure
        .atoms
        .iter()
        .filter(|atom| selection.matches(&atom.name) && chain_id.is_none_or(|c| atom.chain_id == c))
        .map(|atom| (atom.x, atom.y, atom.z))
        .collect()
}

/// Extract coordinates with residue information for per-residue RMSD.
///
/// Returns tuples of ((chain_id, residue_seq, residue_name), (x, y, z)).
pub fn extract_coords_with_residue_info(
    structure: &PdbStructure,
    selection: &AtomSelection,
    chain_id: Option<&str>,
) -> Vec<CoordWithResidue> {
    structure
        .atoms
        .iter()
        .filter(|atom| selection.matches(&atom.name) && chain_id.is_none_or(|c| atom.chain_id == c))
        .map(|atom| {
            (
                (
                    atom.chain_id.clone(),
                    atom.residue_seq,
                    atom.residue_name.clone(),
                ),
                (atom.x, atom.y, atom.z),
            )
        })
        .collect()
}

/// Compute the centroid (center of mass) of a set of coordinates.
///
/// # Arguments
/// * `coords` - Vector of (x, y, z) coordinates
///
/// # Returns
/// The centroid as (x, y, z), or error if coords is empty
pub fn compute_centroid(coords: &[(f64, f64, f64)]) -> Result<(f64, f64, f64), PdbError> {
    if coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "Cannot compute centroid of empty coordinate set".to_string(),
        ));
    }

    let n = coords.len() as f64;
    let cx = coords.iter().map(|c| c.0).sum::<f64>() / n;
    let cy = coords.iter().map(|c| c.1).sum::<f64>() / n;
    let cz = coords.iter().map(|c| c.2).sum::<f64>() / n;

    Ok((cx, cy, cz))
}

/// Center coordinates by subtracting the centroid.
///
/// # Arguments
/// * `coords` - Vector of (x, y, z) coordinates
///
/// # Returns
/// Tuple of (centered_coords, centroid)
pub fn center_coords(coords: &[Point3D]) -> Result<(Vec<Point3D>, Point3D), PdbError> {
    let centroid = compute_centroid(coords)?;

    let centered: Vec<_> = coords
        .iter()
        .map(|(x, y, z)| (x - centroid.0, y - centroid.1, z - centroid.2))
        .collect();

    Ok((centered, centroid))
}

/// Apply a rotation matrix and translation to a set of coordinates.
///
/// The transformation is applied as: new_coord = R * coord + t
///
/// # Arguments
/// * `coords` - Vector of (x, y, z) coordinates
/// * `rotation` - 3x3 rotation matrix as nested arrays
/// * `translation` - Translation vector [tx, ty, tz]
///
/// # Returns
/// Transformed coordinates
pub fn apply_transform_to_coords(
    coords: &[(f64, f64, f64)],
    rotation: &[[f64; 3]; 3],
    translation: &[f64; 3],
) -> Vec<(f64, f64, f64)> {
    coords
        .iter()
        .map(|(x, y, z)| {
            let rx = rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z;
            let ry = rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z;
            let rz = rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z;
            (
                rx + translation[0],
                ry + translation[1],
                rz + translation[2],
            )
        })
        .collect()
}

/// Apply a rotation matrix and translation to a structure, returning a new structure.
///
/// All atoms in the structure are transformed.
///
/// # Arguments
/// * `structure` - The structure to transform
/// * `rotation` - 3x3 rotation matrix as nested arrays
/// * `translation` - Translation vector [tx, ty, tz]
///
/// # Returns
/// New structure with transformed coordinates
pub fn apply_transform(
    structure: &PdbStructure,
    rotation: &[[f64; 3]; 3],
    translation: &[f64; 3],
) -> PdbStructure {
    let mut transformed = structure.clone();

    for atom in &mut transformed.atoms {
        let x = atom.x;
        let y = atom.y;
        let z = atom.z;

        atom.x = rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z + translation[0];
        atom.y = rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z + translation[1];
        atom.z = rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z + translation[2];
    }

    transformed
}

/// Translate a structure by a given vector.
///
/// # Arguments
/// * `structure` - The structure to translate
/// * `tx` - Translation in x
/// * `ty` - Translation in y
/// * `tz` - Translation in z
///
/// # Returns
/// New structure with translated coordinates
pub fn translate_structure(structure: &PdbStructure, tx: f64, ty: f64, tz: f64) -> PdbStructure {
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    apply_transform(structure, &identity, &[tx, ty, tz])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_selection_ca_only() {
        let sel = AtomSelection::CaOnly;
        assert!(sel.matches("CA"));
        assert!(sel.matches(" CA "));
        assert!(!sel.matches("N"));
        assert!(!sel.matches("C"));
    }

    #[test]
    fn test_atom_selection_backbone() {
        let sel = AtomSelection::Backbone;
        assert!(sel.matches("N"));
        assert!(sel.matches("CA"));
        assert!(sel.matches("C"));
        assert!(sel.matches("O"));
        assert!(!sel.matches("CB"));
        assert!(!sel.matches("OG"));
    }

    #[test]
    fn test_atom_selection_all() {
        let sel = AtomSelection::AllAtoms;
        assert!(sel.matches("N"));
        assert!(sel.matches("CA"));
        assert!(sel.matches("CB"));
        assert!(sel.matches("OG1"));
    }

    #[test]
    fn test_atom_selection_custom() {
        let sel = AtomSelection::Custom(vec!["CA".to_string(), "CB".to_string()]);
        assert!(sel.matches("CA"));
        assert!(sel.matches("CB"));
        assert!(!sel.matches("N"));
        assert!(!sel.matches("C"));
    }

    #[test]
    fn test_compute_centroid() {
        let coords = vec![(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 1.0, 0.0)];
        let centroid = compute_centroid(&coords).unwrap();
        assert!((centroid.0 - 1.0).abs() < 1e-10);
        assert!((centroid.1 - 1.0 / 3.0).abs() < 1e-10);
        assert!((centroid.2 - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_centroid_empty() {
        let coords: Vec<(f64, f64, f64)> = vec![];
        assert!(compute_centroid(&coords).is_err());
    }

    #[test]
    fn test_center_coords() {
        let coords = vec![(1.0, 2.0, 3.0), (3.0, 4.0, 5.0)];
        let (centered, centroid) = center_coords(&coords).unwrap();

        assert!((centroid.0 - 2.0).abs() < 1e-10);
        assert!((centroid.1 - 3.0).abs() < 1e-10);
        assert!((centroid.2 - 4.0).abs() < 1e-10);

        assert!((centered[0].0 - (-1.0)).abs() < 1e-10);
        assert!((centered[0].1 - (-1.0)).abs() < 1e-10);
        assert!((centered[0].2 - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_apply_transform_identity() {
        let coords = vec![(1.0, 2.0, 3.0)];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let translation = [0.0, 0.0, 0.0];

        let transformed = apply_transform_to_coords(&coords, &identity, &translation);
        assert!((transformed[0].0 - 1.0).abs() < 1e-10);
        assert!((transformed[0].1 - 2.0).abs() < 1e-10);
        assert!((transformed[0].2 - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_transform_translation() {
        let coords = vec![(0.0, 0.0, 0.0)];
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let translation = [1.0, 2.0, 3.0];

        let transformed = apply_transform_to_coords(&coords, &identity, &translation);
        assert!((transformed[0].0 - 1.0).abs() < 1e-10);
        assert!((transformed[0].1 - 2.0).abs() < 1e-10);
        assert!((transformed[0].2 - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_transform_rotation_90_z() {
        // 90 degree rotation around z-axis
        let coords = vec![(1.0, 0.0, 0.0)];
        let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let translation = [0.0, 0.0, 0.0];

        let transformed = apply_transform_to_coords(&coords, &rotation, &translation);
        assert!((transformed[0].0 - 0.0).abs() < 1e-10);
        assert!((transformed[0].1 - 1.0).abs() < 1e-10);
        assert!((transformed[0].2 - 0.0).abs() < 1e-10);
    }
}
