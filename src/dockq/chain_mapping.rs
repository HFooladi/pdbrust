//! Chain correspondence detection between model and native structures.
//!
//! Automatically determines which chains in the model correspond to which
//! chains in the native structure using sequence alignment.

use std::collections::HashMap;

use crate::core::PdbStructure;
use crate::error::PdbError;

use super::contacts::ResidueKey;
use super::sequence_align::{AlignmentParams, align_sequences, sequence_identity};

/// Extract the amino acid sequence for a chain from atom records.
///
/// Returns residue names in order of residue sequence number, deduplicating
/// by residue_seq.
pub fn extract_chain_sequence(structure: &PdbStructure, chain_id: &str) -> Vec<String> {
    let mut residues: Vec<(i32, String)> = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for atom in &structure.atoms {
        if atom.chain_id == chain_id
            && !atom.is_hetatm
            && atom.name.trim() == "CA"
            && seen.insert(atom.residue_seq)
        {
            residues.push((atom.residue_seq, atom.residue_name.clone()));
        }
    }

    residues.sort_by_key(|(seq, _)| *seq);
    residues.into_iter().map(|(_, name)| name).collect()
}

/// Find the best chain mapping between model and native structures.
///
/// Uses sequence identity to match chains. For small numbers of chains (<=8),
/// tries all permutations to find the optimal assignment. For larger complexes,
/// uses a greedy approach.
pub fn find_chain_mapping(
    model: &PdbStructure,
    native: &PdbStructure,
) -> Result<Vec<(String, String)>, PdbError> {
    let model_chains = model.get_chain_ids();
    let native_chains = native.get_chain_ids();

    if model_chains.is_empty() || native_chains.is_empty() {
        return Err(PdbError::NoChainMapping(
            "One or both structures have no chains".to_string(),
        ));
    }

    // Extract sequences for all chains
    let model_seqs: Vec<(String, Vec<String>)> = model_chains
        .iter()
        .map(|c| (c.clone(), extract_chain_sequence(model, c)))
        .filter(|(_, seq)| !seq.is_empty())
        .collect();

    let native_seqs: Vec<(String, Vec<String>)> = native_chains
        .iter()
        .map(|c| (c.clone(), extract_chain_sequence(native, c)))
        .filter(|(_, seq)| !seq.is_empty())
        .collect();

    if model_seqs.is_empty() || native_seqs.is_empty() {
        return Err(PdbError::NoChainMapping(
            "No protein chains with CA atoms found".to_string(),
        ));
    }

    // Compute pairwise sequence identity
    let mut similarity_matrix: Vec<Vec<f64>> = Vec::new();
    for (_, model_seq) in &model_seqs {
        let mut row = Vec::new();
        for (_, native_seq) in &native_seqs {
            row.push(sequence_identity(model_seq, native_seq));
        }
        similarity_matrix.push(row);
    }

    // Find optimal assignment
    let assignment = if model_seqs.len() <= 8 && native_seqs.len() <= 8 {
        optimal_assignment(&similarity_matrix, model_seqs.len(), native_seqs.len())
    } else {
        greedy_assignment(&similarity_matrix, model_seqs.len(), native_seqs.len())
    };

    let mapping: Vec<(String, String)> = assignment
        .into_iter()
        .map(|(mi, ni)| (model_seqs[mi].0.clone(), native_seqs[ni].0.clone()))
        .collect();

    if mapping.is_empty() {
        return Err(PdbError::NoChainMapping(
            "Could not find any chain correspondence".to_string(),
        ));
    }

    Ok(mapping)
}

/// Build a residue-level mapping between model and native for a chain pair.
///
/// Maps (model_chain, model_resid) -> (native_chain, native_resid) using
/// sequence alignment to establish residue correspondence.
pub fn build_residue_mapping(
    model: &PdbStructure,
    native: &PdbStructure,
    model_chain: &str,
    native_chain: &str,
) -> HashMap<ResidueKey, ResidueKey> {
    let model_residues = extract_residues_with_seq(model, model_chain);
    let native_residues = extract_residues_with_seq(native, native_chain);

    let model_names: Vec<String> = model_residues
        .iter()
        .map(|(_, name)| name.clone())
        .collect();
    let native_names: Vec<String> = native_residues
        .iter()
        .map(|(_, name)| name.clone())
        .collect();

    let alignment = align_sequences(&model_names, &native_names, &AlignmentParams::default());

    let mut mapping = HashMap::new();
    for (i, opt_j) in alignment.mapping.iter().enumerate() {
        if let Some(j) = opt_j {
            let model_key = (model_chain.to_string(), model_residues[i].0);
            let native_key = (native_chain.to_string(), native_residues[*j].0);
            mapping.insert(model_key, native_key);
        }
    }

    mapping
}

/// Extract residues with their sequence numbers for a chain.
fn extract_residues_with_seq(structure: &PdbStructure, chain_id: &str) -> Vec<(i32, String)> {
    let mut residues = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for atom in &structure.atoms {
        if atom.chain_id == chain_id
            && !atom.is_hetatm
            && atom.name.trim() == "CA"
            && seen.insert(atom.residue_seq)
        {
            residues.push((atom.residue_seq, atom.residue_name.clone()));
        }
    }

    residues.sort_by_key(|(seq, _)| *seq);
    residues
}

/// Find optimal chain assignment using brute-force permutation search.
///
/// Returns vec of (model_index, native_index) pairs.
fn optimal_assignment(
    similarity: &[Vec<f64>],
    n_model: usize,
    n_native: usize,
) -> Vec<(usize, usize)> {
    let n_assign = n_model.min(n_native);

    if n_assign == 0 {
        return Vec::new();
    }

    // Generate permutations of native indices, pick n_assign
    // For efficiency, use itertools permutations
    use itertools::Itertools;

    let native_indices: Vec<usize> = (0..n_native).collect();

    let mut best_score = f64::NEG_INFINITY;
    let mut best_assignment = Vec::new();

    for perm in native_indices.iter().permutations(n_assign) {
        // Try assigning model[0..n_assign] to native[perm[0..n_assign]]
        let score: f64 = (0..n_assign).map(|i| similarity[i][*perm[i]]).sum();

        if score > best_score {
            best_score = score;
            best_assignment = (0..n_assign).map(|i| (i, *perm[i])).collect();
        }
    }

    // If n_model > n_native, also try assigning different model chains
    if n_model > n_native {
        let model_indices: Vec<usize> = (0..n_model).collect();
        for model_perm in model_indices.iter().permutations(n_assign) {
            for perm in native_indices.iter().permutations(n_assign) {
                let score: f64 = (0..n_assign)
                    .map(|i| similarity[*model_perm[i]][*perm[i]])
                    .sum();

                if score > best_score {
                    best_score = score;
                    best_assignment = (0..n_assign).map(|i| (*model_perm[i], *perm[i])).collect();
                }
            }
        }
    }

    best_assignment
}

/// Greedy chain assignment for large complexes.
fn greedy_assignment(
    similarity: &[Vec<f64>],
    n_model: usize,
    n_native: usize,
) -> Vec<(usize, usize)> {
    let mut used_model = vec![false; n_model];
    let mut used_native = vec![false; n_native];
    let mut assignment = Vec::new();

    let n_assign = n_model.min(n_native);

    for _ in 0..n_assign {
        let mut best_score = f64::NEG_INFINITY;
        let mut best_pair = (0, 0);

        for i in 0..n_model {
            if used_model[i] {
                continue;
            }
            for j in 0..n_native {
                if used_native[j] {
                    continue;
                }
                if similarity[i][j] > best_score {
                    best_score = similarity[i][j];
                    best_pair = (i, j);
                }
            }
        }

        if best_score > f64::NEG_INFINITY {
            used_model[best_pair.0] = true;
            used_native[best_pair.1] = true;
            assignment.push(best_pair);
        }
    }

    assignment
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn make_ca_atom(chain_id: &str, residue_seq: i32, residue_name: &str, x: f64) -> Atom {
        Atom {
            serial: residue_seq,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_extract_chain_sequence() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_ca_atom("A", 1, "ALA", 0.0),
            make_ca_atom("A", 2, "GLY", 3.8),
            make_ca_atom("A", 3, "VAL", 7.6),
        ];

        let seq = extract_chain_sequence(&structure, "A");
        assert_eq!(seq, vec!["ALA", "GLY", "VAL"]);
    }

    #[test]
    fn test_chain_mapping_identical() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_ca_atom("A", 1, "ALA", 0.0),
            make_ca_atom("A", 2, "GLY", 3.8),
            make_ca_atom("A", 3, "VAL", 7.6),
            make_ca_atom("B", 1, "LEU", 20.0),
            make_ca_atom("B", 2, "ILE", 23.8),
            make_ca_atom("B", 3, "PHE", 27.6),
        ];

        let mapping = find_chain_mapping(&structure, &structure).unwrap();
        assert_eq!(mapping.len(), 2);
        // A->A and B->B
        assert!(mapping.contains(&("A".to_string(), "A".to_string())));
        assert!(mapping.contains(&("B".to_string(), "B".to_string())));
    }

    #[test]
    fn test_chain_mapping_swapped() {
        let mut native = PdbStructure::new();
        native.atoms = vec![
            make_ca_atom("A", 1, "ALA", 0.0),
            make_ca_atom("A", 2, "GLY", 3.8),
            make_ca_atom("A", 3, "VAL", 7.6),
            make_ca_atom("B", 1, "LEU", 20.0),
            make_ca_atom("B", 2, "ILE", 23.8),
            make_ca_atom("B", 3, "PHE", 27.6),
        ];

        // Model has chains swapped: what was A is now B and vice versa
        let mut model = PdbStructure::new();
        model.atoms = vec![
            make_ca_atom("B", 1, "ALA", 0.0),
            make_ca_atom("B", 2, "GLY", 3.8),
            make_ca_atom("B", 3, "VAL", 7.6),
            make_ca_atom("A", 1, "LEU", 20.0),
            make_ca_atom("A", 2, "ILE", 23.8),
            make_ca_atom("A", 3, "PHE", 27.6),
        ];

        let mapping = find_chain_mapping(&model, &native).unwrap();
        assert_eq!(mapping.len(), 2);
        // Model B -> Native A (same sequence)
        assert!(mapping.contains(&("B".to_string(), "A".to_string())));
        // Model A -> Native B (same sequence)
        assert!(mapping.contains(&("A".to_string(), "B".to_string())));
    }

    #[test]
    fn test_residue_mapping() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            make_ca_atom("A", 1, "ALA", 0.0),
            make_ca_atom("A", 2, "GLY", 3.8),
            make_ca_atom("A", 3, "VAL", 7.6),
        ];

        let mapping = build_residue_mapping(&structure, &structure, "A", "A");
        assert_eq!(mapping.len(), 3);
        assert_eq!(
            mapping.get(&("A".to_string(), 1)),
            Some(&("A".to_string(), 1))
        );
    }

    #[test]
    fn test_greedy_assignment() {
        let similarity = vec![vec![1.0, 0.3], vec![0.2, 0.9]];
        let result = greedy_assignment(&similarity, 2, 2);
        assert_eq!(result.len(), 2);
        // Should assign 0->0 (1.0) and 1->1 (0.9)
        assert!(result.contains(&(0, 0)));
        assert!(result.contains(&(1, 1)));
    }
}
