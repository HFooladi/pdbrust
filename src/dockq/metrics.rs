//! iRMSD, LRMSD, and DockQ score computation.
//!
//! Computes interface RMSD (iRMSD), ligand RMSD (LRMSD), and the combined
//! DockQ score for protein-protein interface quality assessment.

use std::collections::{HashMap, HashSet};

use crate::core::PdbStructure;
use crate::error::PdbError;
use crate::geometry::{
    AtomSelection, apply_transform_to_coords, extract_coords_by_selection, rmsd_from_coords,
    superpose_coords,
};

type CoordPair = (Vec<(f64, f64, f64)>, Vec<(f64, f64, f64)>);

use super::contacts::{ResidueKey, find_interface_residues};

/// Compute the DockQ score from its three components.
///
/// DockQ = (fnat + 1/(1+(iRMSD/1.5)^2) + 1/(1+(LRMSD/8.5)^2)) / 3
pub fn compute_dockq(fnat: f64, irmsd: f64, lrmsd: f64) -> f64 {
    let irmsd_score = 1.0 / (1.0 + (irmsd / 1.5_f64).powi(2));
    let lrmsd_score = 1.0 / (1.0 + (lrmsd / 8.5_f64).powi(2));
    (fnat + irmsd_score + lrmsd_score) / 3.0
}

/// Compute interface RMSD (iRMSD) between model and native.
///
/// 1. Identify interface residues in native (residues within threshold of the other chain)
/// 2. Collect backbone atoms for these interface residues from both structures
/// 3. Superpose and compute RMSD
#[allow(clippy::too_many_arguments)]
pub fn compute_irmsd(
    model: &PdbStructure,
    native: &PdbStructure,
    native_chain_a: &str,
    native_chain_b: &str,
    _model_chain_a: &str,
    _model_chain_b: &str,
    residue_mapping: &HashMap<ResidueKey, ResidueKey>,
    interface_threshold: f64,
) -> Result<f64, PdbError> {
    // Find interface residues in native
    let native_interface =
        find_interface_residues(native, native_chain_a, native_chain_b, interface_threshold);

    if native_interface.is_empty() {
        return Err(PdbError::NoInterfaceContacts(format!(
            "No interface residues between chains {} and {} within {} A",
            native_chain_a, native_chain_b, interface_threshold
        )));
    }

    // Build reverse mapping: native_key -> model_key
    let reverse_mapping: HashMap<&ResidueKey, &ResidueKey> = residue_mapping
        .iter()
        .map(|(model_key, native_key)| (native_key, model_key))
        .collect();

    // Collect backbone coords for interface residues
    let (native_coords, model_coords) =
        collect_matched_backbone_coords(native, model, &native_interface, &reverse_mapping)?;

    if native_coords.len() < 3 {
        return Err(PdbError::InsufficientAtoms(
            "Need at least 3 interface backbone atoms for iRMSD".to_string(),
        ));
    }

    // Superpose and compute RMSD
    let (_, result) = superpose_coords(&model_coords, &native_coords)?;
    Ok(result.rmsd)
}

/// Compute ligand RMSD (LRMSD) between model and native.
///
/// 1. Determine receptor (larger chain) and ligand (smaller chain)
/// 2. Superpose model receptor onto native receptor
/// 3. Apply the same transformation to model ligand
/// 4. Compute RMSD between transformed model ligand and native ligand
pub fn compute_lrmsd(
    model: &PdbStructure,
    native: &PdbStructure,
    native_chain_a: &str,
    native_chain_b: &str,
    model_chain_a: &str,
    model_chain_b: &str,
) -> Result<f64, PdbError> {
    // Determine receptor and ligand by size
    let native_a_count = count_residues(native, native_chain_a);
    let native_b_count = count_residues(native, native_chain_b);

    let (native_receptor, native_ligand, model_receptor, model_ligand) =
        if native_a_count >= native_b_count {
            (native_chain_a, native_chain_b, model_chain_a, model_chain_b)
        } else {
            (native_chain_b, native_chain_a, model_chain_b, model_chain_a)
        };

    let selection = AtomSelection::Backbone;

    // Extract receptor backbone coords
    let native_receptor_coords =
        extract_coords_by_selection(native, &selection, Some(native_receptor));
    let model_receptor_coords =
        extract_coords_by_selection(model, &selection, Some(model_receptor));

    if native_receptor_coords.is_empty() || model_receptor_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "No backbone atoms in receptor chain".to_string(),
        ));
    }

    if native_receptor_coords.len() != model_receptor_coords.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: native_receptor_coords.len(),
            found: model_receptor_coords.len(),
        });
    }

    if model_receptor_coords.len() < 3 {
        return Err(PdbError::InsufficientAtoms(
            "Need at least 3 receptor backbone atoms for LRMSD".to_string(),
        ));
    }

    // Superpose model receptor onto native receptor
    let (_, alignment) = superpose_coords(&model_receptor_coords, &native_receptor_coords)?;

    // Extract ligand backbone coords
    let model_ligand_coords = extract_coords_by_selection(model, &selection, Some(model_ligand));
    let native_ligand_coords = extract_coords_by_selection(native, &selection, Some(native_ligand));

    if native_ligand_coords.is_empty() || model_ligand_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "No backbone atoms in ligand chain".to_string(),
        ));
    }

    if native_ligand_coords.len() != model_ligand_coords.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: native_ligand_coords.len(),
            found: model_ligand_coords.len(),
        });
    }

    // Apply receptor transformation to ligand
    let transformed_ligand = apply_transform_to_coords(
        &model_ligand_coords,
        &alignment.rotation,
        &alignment.translation,
    );

    // Compute RMSD
    rmsd_from_coords(&transformed_ligand, &native_ligand_coords)
}

/// Collect matched backbone atom coordinates for a set of native interface residues.
///
/// Returns (native_coords, model_coords) for backbone atoms that exist in both structures.
fn collect_matched_backbone_coords(
    native: &PdbStructure,
    model: &PdbStructure,
    native_interface_residues: &HashSet<ResidueKey>,
    reverse_mapping: &HashMap<&ResidueKey, &ResidueKey>,
) -> Result<CoordPair, PdbError> {
    let backbone_names = ["N", "CA", "C", "O"];

    // Group native backbone atoms by (chain, resid, atom_name)
    let mut native_atom_map: HashMap<(&str, i32, &str), (f64, f64, f64)> = HashMap::new();
    for atom in &native.atoms {
        if !atom.is_hetatm {
            let key = (atom.chain_id.as_str(), atom.residue_seq);
            if native_interface_residues.contains(&(atom.chain_id.clone(), atom.residue_seq)) {
                let name = atom.name.trim();
                if backbone_names.contains(&name) {
                    native_atom_map.insert(
                        (atom.chain_id.as_str(), atom.residue_seq, name),
                        (atom.x, atom.y, atom.z),
                    );
                    let _ = key; // used in contains check above
                }
            }
        }
    }

    // Group model backbone atoms by (chain, resid, atom_name)
    let mut model_atom_map: HashMap<(&str, i32, &str), (f64, f64, f64)> = HashMap::new();
    for atom in &model.atoms {
        if !atom.is_hetatm {
            let name = atom.name.trim();
            if backbone_names.contains(&name) {
                model_atom_map.insert(
                    (atom.chain_id.as_str(), atom.residue_seq, name),
                    (atom.x, atom.y, atom.z),
                );
            }
        }
    }

    let mut native_coords = Vec::new();
    let mut model_coords = Vec::new();

    // For each native interface residue, find the corresponding model residue
    for native_res_key in native_interface_residues {
        let model_res_key = match reverse_mapping.get(native_res_key) {
            Some(k) => *k,
            None => continue,
        };

        for atom_name in &backbone_names {
            let native_atom_key = (native_res_key.0.as_str(), native_res_key.1, *atom_name);
            let model_atom_key = (model_res_key.0.as_str(), model_res_key.1, *atom_name);

            if let (Some(nc), Some(mc)) = (
                native_atom_map.get(&native_atom_key),
                model_atom_map.get(&model_atom_key),
            ) {
                native_coords.push(*nc);
                model_coords.push(*mc);
            }
        }
    }

    if native_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(
            "No matching backbone atoms found at interface".to_string(),
        ));
    }

    Ok((native_coords, model_coords))
}

fn count_residues(structure: &PdbStructure, chain_id: &str) -> usize {
    let mut seen = HashSet::new();
    for atom in &structure.atoms {
        if atom.chain_id == chain_id && !atom.is_hetatm {
            seen.insert(atom.residue_seq);
        }
    }
    seen.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dockq_perfect() {
        let score = compute_dockq(1.0, 0.0, 0.0);
        assert!((score - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dockq_zero_fnat() {
        let score = compute_dockq(0.0, 0.0, 0.0);
        // 0 + 1/(1+0) + 1/(1+0) = 2/3
        assert!((score - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_dockq_formula() {
        // Known values
        let fnat = 0.5;
        let irmsd = 1.5; // at d0, score = 0.5
        let lrmsd = 8.5; // at d0, score = 0.5
        let expected = (0.5 + 0.5 + 0.5) / 3.0;
        let score = compute_dockq(fnat, irmsd, lrmsd);
        assert!((score - expected).abs() < 1e-10);
    }

    #[test]
    fn test_dockq_quality_boundaries() {
        // Incorrect: < 0.23
        assert!(compute_dockq(0.0, 10.0, 50.0) < 0.23);
        // High: >= 0.80
        assert!(compute_dockq(1.0, 0.5, 1.0) > 0.80);
    }
}
