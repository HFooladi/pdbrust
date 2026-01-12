//! Structure superposition using the Kabsch algorithm.
//!
//! This module provides functions for optimal structural alignment using
//! the Kabsch algorithm, which finds the rotation matrix that minimizes
//! the RMSD between two sets of points.
//!
//! # Algorithm Overview
//!
//! The Kabsch algorithm works by:
//! 1. Centering both point sets at the origin
//! 2. Computing the cross-covariance matrix H = P^T * Q
//! 3. Using SVD to decompose H = U * S * V^T
//! 4. Computing the optimal rotation R = V * U^T
//! 5. Handling reflections (ensuring det(R) = +1)
//!
//! # Example
//!
//! ```rust,ignore
//! use pdbrust::geometry::{align_structures, AtomSelection};
//!
//! let mobile = parse_pdb_file("model1.pdb")?;
//! let target = parse_pdb_file("model2.pdb")?;
//!
//! let (aligned, result) = align_structures(&mobile, &target, AtomSelection::CaOnly)?;
//! println!("RMSD after alignment: {:.2} Angstroms", result.rmsd);
//! aligned.to_file("aligned.pdb")?;
//! ```

use std::collections::HashMap;

use nalgebra::{Matrix3, SVD, Vector3};

use crate::core::PdbStructure;
use crate::error::PdbError;

use super::rmsd::rmsd_from_coords;
use super::transform::{
    AtomSelection, CoordWithResidue, Point3D, apply_transform, center_coords,
    extract_coords_by_selection, extract_coords_with_residue_info,
};

/// Result type for superposition: (aligned_coords, AlignmentResult).
type SuperposeResult = Result<(Vec<Point3D>, AlignmentResult), PdbError>;

/// Map from residue key to list of coordinates with residue info.
type ResidueCoordMap = HashMap<(String, i32), Vec<CoordWithResidue>>;

/// Result of structure alignment containing RMSD and transformation parameters.
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Root mean square deviation after optimal alignment (Angstroms).
    pub rmsd: f64,
    /// 3x3 rotation matrix that aligns mobile to target.
    pub rotation: [[f64; 3]; 3],
    /// Translation vector applied after rotation.
    pub translation: [f64; 3],
    /// Number of atoms used in the alignment calculation.
    pub num_atoms: usize,
}

/// Per-residue RMSD information for flexibility analysis.
#[derive(Debug, Clone)]
pub struct PerResidueRmsd {
    /// Residue identifier: (chain_id, residue_seq).
    pub residue_id: (String, i32),
    /// Residue name (e.g., "ALA", "GLY").
    pub residue_name: String,
    /// RMSD for this residue after alignment (Angstroms).
    pub rmsd: f64,
    /// Number of atoms used for this residue's RMSD.
    pub num_atoms: usize,
}

/// Compute the optimal rotation matrix using the Kabsch algorithm.
///
/// Given two centered point sets P and Q, finds rotation R that minimizes
/// the sum of squared distances: ||Q - R*P||^2
///
/// # Arguments
/// * `p` - Mobile points (must be centered at origin)
/// * `q` - Target points (must be centered at origin)
///
/// # Returns
/// 3x3 rotation matrix as nalgebra Matrix3
fn kabsch_rotation(p: &[Vector3<f64>], q: &[Vector3<f64>]) -> Matrix3<f64> {
    // Build cross-covariance matrix H = P^T * Q
    let mut h = Matrix3::zeros();
    for (pi, qi) in p.iter().zip(q.iter()) {
        h += pi * qi.transpose();
    }

    // SVD decomposition: H = U * S * V^T
    let svd = SVD::new(h, true, true);
    let u = svd.u.expect("SVD should compute U matrix");
    let v_t = svd.v_t.expect("SVD should compute V^T matrix");

    // Compute rotation: R = V * U^T
    let mut rotation = v_t.transpose() * u.transpose();

    // Handle reflection case (det(R) = -1)
    // This ensures we get a proper rotation, not a reflection
    if rotation.determinant() < 0.0 {
        // Flip sign of last column of V
        let mut v = v_t.transpose();
        v.column_mut(2).scale_mut(-1.0);
        rotation = v * u.transpose();
    }

    rotation
}

/// Superpose mobile coordinates onto target coordinates.
///
/// Returns the transformed coordinates and alignment result.
///
/// # Arguments
/// * `mobile_coords` - Coordinates to be aligned
/// * `target_coords` - Reference coordinates
///
/// # Returns
/// Tuple of (aligned_coords, AlignmentResult)
fn superpose_coords(mobile_coords: &[Point3D], target_coords: &[Point3D]) -> SuperposeResult {
    // Validate input
    if mobile_coords.is_empty() || target_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "Cannot superpose empty coordinate sets".to_string(),
        ));
    }

    if mobile_coords.len() != target_coords.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: target_coords.len(),
            found: mobile_coords.len(),
        });
    }

    if mobile_coords.len() < 3 {
        return Err(PdbError::InsufficientAtoms(
            "Need at least 3 atoms for superposition".to_string(),
        ));
    }

    // Center both sets
    let (mobile_centered, mobile_centroid) = center_coords(mobile_coords)?;
    let (target_centered, target_centroid) = center_coords(target_coords)?;

    // Convert to nalgebra vectors
    let mobile_vec: Vec<Vector3<f64>> = mobile_centered
        .iter()
        .map(|(x, y, z)| Vector3::new(*x, *y, *z))
        .collect();
    let target_vec: Vec<Vector3<f64>> = target_centered
        .iter()
        .map(|(x, y, z)| Vector3::new(*x, *y, *z))
        .collect();

    // Compute optimal rotation using Kabsch algorithm
    let rotation = kabsch_rotation(&mobile_vec, &target_vec);

    // Apply transformation: rotate centered mobile, then translate to target centroid
    let mut aligned = Vec::with_capacity(mobile_coords.len());
    for m in &mobile_vec {
        let rotated = rotation * m;
        aligned.push((
            rotated.x + target_centroid.0,
            rotated.y + target_centroid.1,
            rotated.z + target_centroid.2,
        ));
    }

    // Compute RMSD after alignment
    let rmsd = rmsd_from_coords(&aligned, target_coords)?;

    // Convert rotation to array format
    let rotation_array = [
        [rotation[(0, 0)], rotation[(0, 1)], rotation[(0, 2)]],
        [rotation[(1, 0)], rotation[(1, 1)], rotation[(1, 2)]],
        [rotation[(2, 0)], rotation[(2, 1)], rotation[(2, 2)]],
    ];

    // The full transformation is:
    // 1. Translate by -mobile_centroid (center)
    // 2. Rotate by R
    // 3. Translate by +target_centroid
    // Combined translation: R * (-mobile_centroid) + target_centroid
    let neg_mobile = Vector3::new(-mobile_centroid.0, -mobile_centroid.1, -mobile_centroid.2);
    let rotated_neg_mobile = rotation * neg_mobile;
    let translation = [
        rotated_neg_mobile.x + target_centroid.0,
        rotated_neg_mobile.y + target_centroid.1,
        rotated_neg_mobile.z + target_centroid.2,
    ];

    let result = AlignmentResult {
        rmsd,
        rotation: rotation_array,
        translation,
        num_atoms: mobile_coords.len(),
    };

    Ok((aligned, result))
}

/// Align mobile structure onto target and return aligned structure + RMSD.
///
/// Uses the Kabsch algorithm to find optimal rotation and translation
/// that minimizes RMSD between the selected atoms.
///
/// # Arguments
/// * `mobile` - Structure to be aligned (will be transformed)
/// * `target` - Reference structure (unchanged)
/// * `selection` - Atom selection mode for calculating alignment
///
/// # Returns
/// Tuple of (aligned_structure, AlignmentResult)
///
/// # Errors
/// - `AtomCountMismatch` if structures have different numbers of selected atoms
/// - `NoAtomsSelected` if no atoms match the selection criteria
/// - `InsufficientAtoms` if fewer than 3 atoms are selected
///
/// # Example
///
/// ```rust,ignore
/// use pdbrust::geometry::{align_structures, AtomSelection};
///
/// let mobile = parse_pdb_file("model1.pdb")?;
/// let target = parse_pdb_file("reference.pdb")?;
///
/// let (aligned, result) = align_structures(&mobile, &target, AtomSelection::CaOnly)?;
/// println!("RMSD: {:.3} Angstroms ({} atoms)", result.rmsd, result.num_atoms);
/// ```
pub fn align_structures(
    mobile: &PdbStructure,
    target: &PdbStructure,
    selection: AtomSelection,
) -> Result<(PdbStructure, AlignmentResult), PdbError> {
    let mobile_coords = extract_coords_by_selection(mobile, &selection, None);
    let target_coords = extract_coords_by_selection(target, &selection, None);

    if mobile_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in mobile structure",
            selection
        )));
    }

    if target_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in target structure",
            selection
        )));
    }

    // Compute superposition on selected atoms
    let (_, result) = superpose_coords(&mobile_coords, &target_coords)?;

    // Apply transformation to ALL atoms in the mobile structure
    let aligned = apply_transform(mobile, &result.rotation, &result.translation);

    Ok((aligned, result))
}

/// Calculate optimal alignment without creating new structure.
///
/// Returns only the AlignmentResult with transformation parameters.
/// Useful when you only need the RMSD or transformation matrices.
///
/// # Arguments
/// * `mobile` - Structure to be aligned
/// * `target` - Reference structure
/// * `selection` - Atom selection mode for calculating alignment
///
/// # Returns
/// AlignmentResult containing RMSD and transformation parameters
pub fn calculate_alignment(
    mobile: &PdbStructure,
    target: &PdbStructure,
    selection: AtomSelection,
) -> Result<AlignmentResult, PdbError> {
    let mobile_coords = extract_coords_by_selection(mobile, &selection, None);
    let target_coords = extract_coords_by_selection(target, &selection, None);

    if mobile_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in mobile structure",
            selection
        )));
    }

    if target_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in target structure",
            selection
        )));
    }

    let (_, result) = superpose_coords(&mobile_coords, &target_coords)?;
    Ok(result)
}

/// Compute per-residue RMSD after alignment.
///
/// First aligns the structures using the specified atom selection,
/// then computes RMSD for each residue individually. This is useful
/// for identifying flexible regions in protein structures.
///
/// # Arguments
/// * `mobile` - Structure to be aligned
/// * `target` - Reference structure
/// * `selection` - Atom selection mode for alignment and RMSD calculation
///
/// # Returns
/// Vector of PerResidueRmsd containing RMSD for each residue
///
/// # Example
///
/// ```rust,ignore
/// use pdbrust::geometry::{per_residue_rmsd, AtomSelection};
///
/// let mobile = parse_pdb_file("model1.pdb")?;
/// let target = parse_pdb_file("model2.pdb")?;
///
/// let per_res = per_residue_rmsd(&mobile, &target, AtomSelection::CaOnly)?;
/// for r in &per_res {
///     if r.rmsd > 2.0 {  // Flexible residues
///         println!("{}{}: {:.2} Angstroms", r.residue_id.0, r.residue_id.1, r.rmsd);
///     }
/// }
/// ```
pub fn per_residue_rmsd(
    mobile: &PdbStructure,
    target: &PdbStructure,
    selection: AtomSelection,
) -> Result<Vec<PerResidueRmsd>, PdbError> {
    // First align the structures
    let (aligned, _) = align_structures(mobile, target, selection.clone())?;

    // Extract coordinates with residue info from aligned and target
    let aligned_with_res = extract_coords_with_residue_info(&aligned, &selection, None);
    let target_with_res = extract_coords_with_residue_info(target, &selection, None);

    // Group by residue
    let mut aligned_by_residue: ResidueCoordMap = HashMap::new();
    for item in aligned_with_res {
        let key = (item.0.0.clone(), item.0.1);
        aligned_by_residue.entry(key).or_default().push(item);
    }

    let mut target_by_residue: ResidueCoordMap = HashMap::new();
    for item in target_with_res {
        let key = (item.0.0.clone(), item.0.1);
        target_by_residue.entry(key).or_default().push(item);
    }

    // Compute per-residue RMSD
    let mut results = Vec::new();

    for (residue_key, aligned_atoms) in &aligned_by_residue {
        if let Some(target_atoms) = target_by_residue.get(residue_key) {
            if aligned_atoms.len() == target_atoms.len() && !aligned_atoms.is_empty() {
                let aligned_coords: Vec<_> = aligned_atoms.iter().map(|a| a.1).collect();
                let target_coords: Vec<_> = target_atoms.iter().map(|a| a.1).collect();

                if let Ok(rmsd) = rmsd_from_coords(&aligned_coords, &target_coords) {
                    results.push(PerResidueRmsd {
                        residue_id: residue_key.clone(),
                        residue_name: aligned_atoms[0].0.2.clone(),
                        rmsd,
                        num_atoms: aligned_atoms.len(),
                    });
                }
            }
        }
    }

    // Sort by chain and residue number
    results.sort_by(|a, b| {
        a.residue_id
            .0
            .cmp(&b.residue_id.0)
            .then(a.residue_id.1.cmp(&b.residue_id.1))
    });

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_atom(x: f64, y: f64, z: f64, name: &str, residue_seq: i32, chain_id: &str) -> Atom {
        Atom {
            serial: residue_seq,
            name: name.to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
            ins_code: None,
        }
    }

    fn create_linear_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(0.0, 0.0, 0.0, "CA", 1, "A"),
            create_atom(3.8, 0.0, 0.0, "CA", 2, "A"),
            create_atom(7.6, 0.0, 0.0, "CA", 3, "A"),
            create_atom(11.4, 0.0, 0.0, "CA", 4, "A"),
        ];
        structure
    }

    fn translate_structure(structure: &PdbStructure, tx: f64, ty: f64, tz: f64) -> PdbStructure {
        let mut translated = structure.clone();
        for atom in &mut translated.atoms {
            atom.x += tx;
            atom.y += ty;
            atom.z += tz;
        }
        translated
    }

    #[test]
    fn test_kabsch_identity() {
        let points = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];

        let rotation = kabsch_rotation(&points, &points);

        // Should be close to identity matrix
        assert!((rotation[(0, 0)] - 1.0).abs() < 1e-10);
        assert!((rotation[(1, 1)] - 1.0).abs() < 1e-10);
        assert!((rotation[(2, 2)] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_kabsch_rotation_is_proper() {
        // Create two point sets where one is a rotation of the other
        let p = vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        // 90 degree rotation around z-axis
        let q = vec![
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let rotation = kabsch_rotation(&p, &q);

        // Determinant should be +1 (proper rotation)
        assert!((rotation.determinant() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_superpose_identical() {
        let coords = vec![
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ];

        let (aligned, result) = superpose_coords(&coords, &coords).unwrap();

        assert!(result.rmsd < 1e-10);
        assert_eq!(result.num_atoms, 4);

        // Aligned coordinates should match original
        for (a, o) in aligned.iter().zip(coords.iter()) {
            assert!((a.0 - o.0).abs() < 1e-10);
            assert!((a.1 - o.1).abs() < 1e-10);
            assert!((a.2 - o.2).abs() < 1e-10);
        }
    }

    #[test]
    fn test_superpose_translated() {
        let coords1 = vec![
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ];

        // Translated by (10, 10, 10)
        let coords2: Vec<_> = coords1
            .iter()
            .map(|(x, y, z)| (x + 10.0, y + 10.0, z + 10.0))
            .collect();

        let (aligned, result) = superpose_coords(&coords2, &coords1).unwrap();

        // RMSD should be ~0 after alignment
        assert!(result.rmsd < 1e-10);

        // Aligned coordinates should match target
        for (a, t) in aligned.iter().zip(coords1.iter()) {
            assert!((a.0 - t.0).abs() < 1e-10);
            assert!((a.1 - t.1).abs() < 1e-10);
            assert!((a.2 - t.2).abs() < 1e-10);
        }
    }

    #[test]
    fn test_superpose_insufficient_atoms() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let result = superpose_coords(&coords, &coords);
        assert!(matches!(result, Err(PdbError::InsufficientAtoms(_))));
    }

    #[test]
    fn test_superpose_mismatched_length() {
        let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)];
        let coords2 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let result = superpose_coords(&coords1, &coords2);
        assert!(matches!(result, Err(PdbError::AtomCountMismatch { .. })));
    }

    #[test]
    fn test_align_structures_identical() {
        let structure = create_linear_structure();
        let (aligned, result) =
            align_structures(&structure, &structure, AtomSelection::CaOnly).unwrap();

        assert!(result.rmsd < 1e-10);
        assert_eq!(result.num_atoms, 4);
        assert_eq!(aligned.atoms.len(), structure.atoms.len());
    }

    #[test]
    fn test_align_structures_translated() {
        let target = create_linear_structure();
        let mobile = translate_structure(&target, 50.0, 50.0, 50.0);

        let (aligned, result) = align_structures(&mobile, &target, AtomSelection::CaOnly).unwrap();

        // RMSD should be ~0 after alignment
        assert!(result.rmsd < 1e-6);

        // Aligned atoms should be close to target
        for (aligned_atom, target_atom) in aligned.atoms.iter().zip(target.atoms.iter()) {
            assert!((aligned_atom.x - target_atom.x).abs() < 1e-6);
            assert!((aligned_atom.y - target_atom.y).abs() < 1e-6);
            assert!((aligned_atom.z - target_atom.z).abs() < 1e-6);
        }
    }

    #[test]
    fn test_calculate_alignment() {
        let target = create_linear_structure();
        let mobile = translate_structure(&target, 10.0, 20.0, 30.0);

        let result = calculate_alignment(&mobile, &target, AtomSelection::CaOnly).unwrap();

        assert!(result.rmsd < 1e-6);
        assert_eq!(result.num_atoms, 4);
    }

    #[test]
    fn test_per_residue_rmsd() {
        let target = create_linear_structure();
        let mobile = translate_structure(&target, 5.0, 5.0, 5.0);

        let per_res = per_residue_rmsd(&mobile, &target, AtomSelection::CaOnly).unwrap();

        // Should have 4 residues
        assert_eq!(per_res.len(), 4);

        // All per-residue RMSDs should be ~0 after alignment
        for r in &per_res {
            assert!(r.rmsd < 1e-6);
            assert_eq!(r.num_atoms, 1); // Only CA per residue
        }
    }
}
