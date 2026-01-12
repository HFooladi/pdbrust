//! RMSD (Root Mean Square Deviation) calculation.
//!
//! This module provides functions for computing RMSD between coordinate sets
//! or structures. RMSD measures the average distance between atoms after
//! optimal superposition.
//!
//! # Example
//!
//! ```rust,ignore
//! use pdbrust::geometry::{calculate_rmsd, AtomSelection};
//!
//! let structure1 = parse_pdb_file("protein1.pdb")?;
//! let structure2 = parse_pdb_file("protein2.pdb")?;
//!
//! // Calculate RMSD using CA atoms
//! let rmsd = calculate_rmsd(&structure1, &structure2, AtomSelection::CaOnly)?;
//! println!("RMSD: {:.2} Angstroms", rmsd);
//! ```

use crate::core::PdbStructure;
use crate::error::PdbError;

use super::transform::{extract_coords_by_selection, AtomSelection};

/// Calculate RMSD directly from two coordinate sets.
///
/// This computes RMSD without any alignment - coordinates are compared
/// directly in their current positions. For optimal superposition,
/// use the alignment functions in the `superposition` module first.
///
/// # Formula
///
/// RMSD = sqrt(1/N * sum((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2))
///
/// # Arguments
/// * `coords1` - First set of coordinates
/// * `coords2` - Second set of coordinates (must have same length)
///
/// # Returns
/// RMSD in the same units as the input coordinates (typically Angstroms)
///
/// # Errors
/// - `AtomCountMismatch` if coordinate sets have different lengths
/// - `NoAtomsSelected` if either coordinate set is empty
///
/// # Example
///
/// ```rust
/// use pdbrust::geometry::rmsd_from_coords;
///
/// let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
/// let coords2 = vec![(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]; // 1 Angstrom difference at atom 2
///
/// let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
/// // RMSD = sqrt((0 + 1) / 2) = sqrt(0.5) ≈ 0.707
/// ```
pub fn rmsd_from_coords(
    coords1: &[(f64, f64, f64)],
    coords2: &[(f64, f64, f64)],
) -> Result<f64, PdbError> {
    if coords1.is_empty() || coords2.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "Cannot compute RMSD with empty coordinate sets".to_string(),
        ));
    }

    if coords1.len() != coords2.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: coords1.len(),
            found: coords2.len(),
        });
    }

    let n = coords1.len() as f64;
    let sum_sq: f64 = coords1
        .iter()
        .zip(coords2.iter())
        .map(|((x1, y1, z1), (x2, y2, z2))| {
            let dx = x1 - x2;
            let dy = y1 - y2;
            let dz = z1 - z2;
            dx * dx + dy * dy + dz * dz
        })
        .sum();

    Ok((sum_sq / n).sqrt())
}

/// Calculate RMSD between two structures based on atom selection.
///
/// This function extracts coordinates based on the selection criteria
/// and computes RMSD without alignment. For aligned RMSD, use
/// `aligned_rmsd` or the alignment functions.
///
/// # Arguments
/// * `structure1` - First structure
/// * `structure2` - Second structure
/// * `selection` - Atom selection criteria (CA only, backbone, all, or custom)
///
/// # Returns
/// RMSD in Angstroms
///
/// # Errors
/// - `AtomCountMismatch` if structures have different numbers of selected atoms
/// - `NoAtomsSelected` if no atoms match the selection criteria
///
/// # Example
///
/// ```rust,ignore
/// use pdbrust::geometry::{calculate_rmsd, AtomSelection};
///
/// let s1 = parse_pdb_file("model1.pdb")?;
/// let s2 = parse_pdb_file("model2.pdb")?;
///
/// // RMSD using CA atoms only (default)
/// let rmsd_ca = calculate_rmsd(&s1, &s2, AtomSelection::CaOnly)?;
///
/// // RMSD using backbone atoms
/// let rmsd_bb = calculate_rmsd(&s1, &s2, AtomSelection::Backbone)?;
/// ```
pub fn calculate_rmsd(
    structure1: &PdbStructure,
    structure2: &PdbStructure,
    selection: AtomSelection,
) -> Result<f64, PdbError> {
    let coords1 = extract_coords_by_selection(structure1, &selection, None);
    let coords2 = extract_coords_by_selection(structure2, &selection, None);

    if coords1.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in first structure",
            selection
        )));
    }

    if coords2.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in second structure",
            selection
        )));
    }

    rmsd_from_coords(&coords1, &coords2)
}

/// Calculate RMSD between two structures for a specific chain.
///
/// # Arguments
/// * `structure1` - First structure
/// * `structure2` - Second structure
/// * `selection` - Atom selection criteria
/// * `chain_id` - Chain identifier to compare
///
/// # Returns
/// RMSD in Angstroms for the specified chain
pub fn calculate_rmsd_chain(
    structure1: &PdbStructure,
    structure2: &PdbStructure,
    selection: AtomSelection,
    chain_id: &str,
) -> Result<f64, PdbError> {
    let coords1 = extract_coords_by_selection(structure1, &selection, Some(chain_id));
    let coords2 = extract_coords_by_selection(structure2, &selection, Some(chain_id));

    if coords1.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in chain {} of first structure",
            selection, chain_id
        )));
    }

    if coords2.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in chain {} of second structure",
            selection, chain_id
        )));
    }

    rmsd_from_coords(&coords1, &coords2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_atom(x: f64, y: f64, z: f64, name: &str, residue_seq: i32) -> Atom {
        Atom {
            serial: residue_seq,
            name: name.to_string(),
            alt_loc: None,
            residue_name: "ALA".to_string(),
            chain_id: "A".to_string(),
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

    fn create_test_structure(offset_x: f64) -> PdbStructure {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(0.0 + offset_x, 0.0, 0.0, "CA", 1),
            create_atom(3.8 + offset_x, 0.0, 0.0, "CA", 2),
            create_atom(7.6 + offset_x, 0.0, 0.0, "CA", 3),
        ];
        structure
    }

    #[test]
    fn test_rmsd_identical_coords() {
        let coords = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)];
        let rmsd = rmsd_from_coords(&coords, &coords).unwrap();
        assert!(rmsd < 1e-10);
    }

    #[test]
    fn test_rmsd_known_displacement() {
        // All atoms displaced by 1.0 in x direction
        let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let coords2 = vec![(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)];
        let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
        assert!((rmsd - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_single_displaced() {
        // One atom at origin, one displaced by sqrt(3) ≈ 1.732
        let coords1 = vec![(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)];
        let coords2 = vec![(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)];
        let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
        // RMSD = sqrt((0 + 3) / 2) = sqrt(1.5) ≈ 1.2247
        assert!((rmsd - (1.5_f64).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_empty_coords() {
        let empty: Vec<(f64, f64, f64)> = vec![];
        let coords = vec![(0.0, 0.0, 0.0)];

        assert!(rmsd_from_coords(&empty, &coords).is_err());
        assert!(rmsd_from_coords(&coords, &empty).is_err());
    }

    #[test]
    fn test_rmsd_mismatched_length() {
        let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let coords2 = vec![(0.0, 0.0, 0.0)];

        let result = rmsd_from_coords(&coords1, &coords2);
        assert!(matches!(result, Err(PdbError::AtomCountMismatch { .. })));
    }

    #[test]
    fn test_calculate_rmsd_identical_structures() {
        let structure = create_test_structure(0.0);
        let rmsd = calculate_rmsd(&structure, &structure, AtomSelection::CaOnly).unwrap();
        assert!(rmsd < 1e-10);
    }

    #[test]
    fn test_calculate_rmsd_translated_structure() {
        let structure1 = create_test_structure(0.0);
        let structure2 = create_test_structure(1.0); // Shifted by 1.0 in x

        let rmsd = calculate_rmsd(&structure1, &structure2, AtomSelection::CaOnly).unwrap();
        assert!((rmsd - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_calculate_rmsd_no_matching_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![create_atom(0.0, 0.0, 0.0, "CB", 1)]; // CB, not CA

        let result = calculate_rmsd(&structure, &structure, AtomSelection::CaOnly);
        assert!(matches!(result, Err(PdbError::NoAtomsSelected(_))));
    }
}
