//! LDDT (Local Distance Difference Test) calculation.
//!
//! LDDT is a superposition-free metric for comparing protein structures,
//! widely used in AlphaFold (pLDDT) and CASP evaluations. It measures the
//! fraction of inter-atomic distances that are preserved within specified
//! thresholds.
//!
//! # Algorithm
//!
//! 1. For each atom in the reference structure, find all atoms within the
//!    inclusion radius (default: 15Å)
//! 2. For each pair (i,j) in reference, compute distance d_ref
//! 3. Find the same pair in the model and compute d_model
//! 4. Check if |d_ref - d_model| < threshold for each threshold
//! 5. Score = average fraction of preserved distances across all thresholds
//!
//! The default thresholds are 0.5Å, 1.0Å, 2.0Å, and 4.0Å.
//!
//! # Key Properties
//!
//! - **Superposition-free**: LDDT is invariant to rotation and translation
//! - **Local**: Focuses on local structure quality, not global fold
//! - **Range**: 0.0 (poor) to 1.0 (perfect)
//!
//! # Example
//!
//! ```rust,ignore
//! use pdbrust::geometry::{calculate_lddt, LddtOptions, AtomSelection};
//!
//! let model = parse_pdb_file("model.pdb")?;
//! let reference = parse_pdb_file("reference.pdb")?;
//!
//! // Calculate LDDT with default options (CA atoms, 15Å radius)
//! let result = model.lddt_to(&reference)?;
//! println!("LDDT: {:.4}", result.score);
//!
//! // With custom options
//! let options = LddtOptions::default().with_inclusion_radius(10.0);
//! let result = model.lddt_to_with_options(&reference, AtomSelection::Backbone, options)?;
//! ```
//!
//! # References
//!
//! - Mariani V, et al. (2013) "lDDT: a local superposition-free score for
//!   comparing protein structures and models using distance difference tests"
//!   Bioinformatics 29(21):2722-2728

use crate::core::PdbStructure;
use crate::error::PdbError;

use super::transform::{AtomSelection, CoordWithResidue, extract_coords_with_residue_info};

/// Options for LDDT calculation.
///
/// # Example
///
/// ```rust
/// use pdbrust::geometry::LddtOptions;
///
/// // Default options
/// let options = LddtOptions::default();
/// assert_eq!(options.inclusion_radius, 15.0);
///
/// // Custom options
/// let options = LddtOptions::default()
///     .with_inclusion_radius(10.0)
///     .with_thresholds(vec![0.5, 1.0, 2.0, 4.0, 8.0]);
/// ```
#[derive(Debug, Clone)]
pub struct LddtOptions {
    /// Maximum distance to consider for local comparisons (default: 15.0 Å).
    ///
    /// Pairs of atoms with reference distance > inclusion_radius are ignored.
    pub inclusion_radius: f64,

    /// Distance difference thresholds (default: [0.5, 1.0, 2.0, 4.0] Å).
    ///
    /// For each threshold, count how many distance pairs are preserved
    /// (i.e., |d_ref - d_model| < threshold). The final score is the
    /// average across all thresholds.
    pub thresholds: Vec<f64>,
}

impl Default for LddtOptions {
    fn default() -> Self {
        Self {
            inclusion_radius: 15.0,
            thresholds: vec![0.5, 1.0, 2.0, 4.0],
        }
    }
}

impl LddtOptions {
    /// Create new options with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the inclusion radius.
    ///
    /// Only pairs of atoms with distance <= inclusion_radius in the
    /// reference structure are considered.
    pub fn with_inclusion_radius(mut self, radius: f64) -> Self {
        self.inclusion_radius = radius;
        self
    }

    /// Set custom thresholds.
    ///
    /// The LDDT score is the average fraction of preserved distances
    /// across all thresholds.
    pub fn with_thresholds(mut self, thresholds: Vec<f64>) -> Self {
        self.thresholds = thresholds;
        self
    }
}

/// Result of LDDT calculation.
#[derive(Debug, Clone)]
pub struct LddtResult {
    /// Global LDDT score (0.0 to 1.0).
    ///
    /// This is the average fraction of preserved distances across all
    /// thresholds, averaged over all residues.
    pub score: f64,

    /// Number of distance pairs evaluated.
    pub num_pairs: usize,

    /// Score for each threshold.
    ///
    /// Each value is the fraction of pairs preserved at that threshold.
    pub per_threshold_scores: Vec<f64>,

    /// Number of residues evaluated.
    pub num_residues: usize,
}

/// Per-residue LDDT information.
#[derive(Debug, Clone)]
pub struct PerResidueLddt {
    /// Residue identifier as (chain_id, residue_seq).
    pub residue_id: (String, i32),

    /// Residue name (e.g., "ALA", "GLY").
    pub residue_name: String,

    /// LDDT score for this residue (0.0 to 1.0).
    pub score: f64,

    /// Number of distance pairs involving this residue.
    pub num_pairs: usize,
}

/// Calculate LDDT between two structures.
///
/// This function computes the Local Distance Difference Test score between
/// a model structure and a reference structure. The score measures how well
/// local inter-atomic distances are preserved.
///
/// # Arguments
///
/// * `model` - The model structure to evaluate
/// * `reference` - The reference structure (ground truth)
/// * `selection` - Atom selection criteria (e.g., CA only, backbone)
/// * `options` - LDDT calculation options (thresholds, inclusion radius)
///
/// # Returns
///
/// `LddtResult` containing the global score, per-threshold scores, and counts.
///
/// # Errors
///
/// - `NoAtomsSelected` if no atoms match the selection
/// - `AtomCountMismatch` if structures have different numbers of selected atoms
///
/// # Example
///
/// ```rust,ignore
/// use pdbrust::geometry::{calculate_lddt, LddtOptions, AtomSelection};
///
/// let result = calculate_lddt(&model, &reference, AtomSelection::CaOnly, LddtOptions::default())?;
/// println!("LDDT: {:.4}", result.score);
/// ```
pub fn calculate_lddt(
    model: &PdbStructure,
    reference: &PdbStructure,
    selection: AtomSelection,
    options: LddtOptions,
) -> Result<LddtResult, PdbError> {
    // Extract coordinates with residue info
    let model_coords = extract_coords_with_residue_info(model, &selection, None);
    let ref_coords = extract_coords_with_residue_info(reference, &selection, None);

    if model_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in model structure",
            selection
        )));
    }

    if ref_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in reference structure",
            selection
        )));
    }

    if model_coords.len() != ref_coords.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: ref_coords.len(),
            found: model_coords.len(),
        });
    }

    calculate_lddt_from_coords(&model_coords, &ref_coords, &options)
}

/// Calculate per-residue LDDT scores.
///
/// Returns LDDT scores for each residue individually, useful for
/// identifying poorly modeled regions.
///
/// # Arguments
///
/// * `model` - The model structure to evaluate
/// * `reference` - The reference structure (ground truth)
/// * `selection` - Atom selection criteria
/// * `options` - LDDT calculation options
///
/// # Returns
///
/// Vector of `PerResidueLddt` for each residue.
///
/// # Example
///
/// ```rust,ignore
/// let per_res = per_residue_lddt(&model, &reference, AtomSelection::CaOnly, LddtOptions::default())?;
/// for r in per_res.iter().filter(|r| r.score < 0.7) {
///     println!("Low LDDT: {}{} = {:.2}", r.residue_id.0, r.residue_id.1, r.score);
/// }
/// ```
pub fn per_residue_lddt(
    model: &PdbStructure,
    reference: &PdbStructure,
    selection: AtomSelection,
    options: LddtOptions,
) -> Result<Vec<PerResidueLddt>, PdbError> {
    // Extract coordinates with residue info
    let model_coords = extract_coords_with_residue_info(model, &selection, None);
    let ref_coords = extract_coords_with_residue_info(reference, &selection, None);

    if model_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in model structure",
            selection
        )));
    }

    if ref_coords.is_empty() {
        return Err(PdbError::NoAtomsSelected(format!(
            "No atoms matching {:?} selection in reference structure",
            selection
        )));
    }

    if model_coords.len() != ref_coords.len() {
        return Err(PdbError::AtomCountMismatch {
            expected: ref_coords.len(),
            found: model_coords.len(),
        });
    }

    per_residue_lddt_from_coords(&model_coords, &ref_coords, &options)
}

/// Calculate Euclidean distance between two 3D points.
#[inline]
fn distance(p1: &(f64, f64, f64), p2: &(f64, f64, f64)) -> f64 {
    let dx = p1.0 - p2.0;
    let dy = p1.1 - p2.1;
    let dz = p1.2 - p2.2;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Calculate LDDT from coordinate arrays with residue information.
fn calculate_lddt_from_coords(
    model_coords: &[CoordWithResidue],
    ref_coords: &[CoordWithResidue],
    options: &LddtOptions,
) -> Result<LddtResult, PdbError> {
    let n = ref_coords.len();
    let inclusion_radius_sq = options.inclusion_radius * options.inclusion_radius;

    // Count preserved distances for each threshold
    let mut preserved_counts: Vec<usize> = vec![0; options.thresholds.len()];
    let mut total_pairs: usize = 0;

    // Count unique residues
    let mut unique_residues = std::collections::HashSet::new();

    // For each pair of atoms (i, j) where i < j
    for i in 0..n {
        let ref_i = &ref_coords[i];
        let model_i = &model_coords[i];

        unique_residues.insert((ref_i.0.0.clone(), ref_i.0.1));

        for j in (i + 1)..n {
            let ref_j = &ref_coords[j];
            let model_j = &model_coords[j];

            // Calculate reference distance
            let d_ref = distance(&ref_i.1, &ref_j.1);

            // Skip pairs outside inclusion radius
            if d_ref * d_ref > inclusion_radius_sq {
                continue;
            }

            // Calculate model distance
            let d_model = distance(&model_i.1, &model_j.1);

            // Check each threshold
            let diff = (d_ref - d_model).abs();
            for (k, threshold) in options.thresholds.iter().enumerate() {
                if diff < *threshold {
                    preserved_counts[k] += 1;
                }
            }

            total_pairs += 1;
        }
    }

    // Calculate per-threshold scores
    let per_threshold_scores: Vec<f64> = if total_pairs > 0 {
        preserved_counts
            .iter()
            .map(|&count| count as f64 / total_pairs as f64)
            .collect()
    } else {
        vec![1.0; options.thresholds.len()] // Perfect score if no pairs
    };

    // Global LDDT is average of per-threshold scores
    let score = if per_threshold_scores.is_empty() {
        1.0
    } else {
        per_threshold_scores.iter().sum::<f64>() / per_threshold_scores.len() as f64
    };

    Ok(LddtResult {
        score,
        num_pairs: total_pairs,
        per_threshold_scores,
        num_residues: unique_residues.len(),
    })
}

/// Calculate per-residue LDDT from coordinate arrays.
fn per_residue_lddt_from_coords(
    model_coords: &[CoordWithResidue],
    ref_coords: &[CoordWithResidue],
    options: &LddtOptions,
) -> Result<Vec<PerResidueLddt>, PdbError> {
    use std::collections::HashMap;

    let n = ref_coords.len();
    let inclusion_radius_sq = options.inclusion_radius * options.inclusion_radius;
    let num_thresholds = options.thresholds.len();

    // Track per-residue statistics: (chain_id, residue_seq) -> (preserved_counts, total_pairs, residue_name)
    let mut residue_stats: HashMap<(String, i32), (Vec<usize>, usize, String)> = HashMap::new();

    // Initialize residue entries
    for ref_coord in ref_coords {
        let key = (ref_coord.0.0.clone(), ref_coord.0.1);
        residue_stats
            .entry(key)
            .or_insert_with(|| (vec![0; num_thresholds], 0, ref_coord.0.2.clone()));
    }

    // For each pair of atoms (i, j)
    for i in 0..n {
        let ref_i = &ref_coords[i];
        let model_i = &model_coords[i];

        for j in (i + 1)..n {
            let ref_j = &ref_coords[j];
            let model_j = &model_coords[j];

            // Calculate reference distance
            let d_ref = distance(&ref_i.1, &ref_j.1);

            // Skip pairs outside inclusion radius
            if d_ref * d_ref > inclusion_radius_sq {
                continue;
            }

            // Calculate model distance
            let d_model = distance(&model_i.1, &model_j.1);
            let diff = (d_ref - d_model).abs();

            // Update statistics for residue i
            let key_i = (ref_i.0.0.clone(), ref_i.0.1);
            if let Some((preserved, total, _)) = residue_stats.get_mut(&key_i) {
                for (k, threshold) in options.thresholds.iter().enumerate() {
                    if diff < *threshold {
                        preserved[k] += 1;
                    }
                }
                *total += 1;
            }

            // Update statistics for residue j
            let key_j = (ref_j.0.0.clone(), ref_j.0.1);
            if let Some((preserved, total, _)) = residue_stats.get_mut(&key_j) {
                for (k, threshold) in options.thresholds.iter().enumerate() {
                    if diff < *threshold {
                        preserved[k] += 1;
                    }
                }
                *total += 1;
            }
        }
    }

    // Convert to PerResidueLddt
    let mut results: Vec<PerResidueLddt> = residue_stats
        .into_iter()
        .map(|(key, (preserved, total, residue_name))| {
            let score = if total > 0 {
                let per_threshold: Vec<f64> =
                    preserved.iter().map(|&p| p as f64 / total as f64).collect();
                per_threshold.iter().sum::<f64>() / num_thresholds as f64
            } else {
                1.0 // Perfect score if no pairs (isolated atom)
            };

            PerResidueLddt {
                residue_id: key,
                residue_name,
                score,
                num_pairs: total,
            }
        })
        .collect();

    // Sort by chain_id then residue_seq
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

    fn create_atom(x: f64, y: f64, z: f64, residue_seq: i32, chain_id: &str) -> Atom {
        Atom {
            serial: residue_seq,
            name: "CA".to_string(),
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
            is_hetatm: false,
        }
    }

    fn create_linear_structure(spacing: f64) -> PdbStructure {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(0.0, 0.0, 0.0, 1, "A"),
            create_atom(spacing, 0.0, 0.0, 2, "A"),
            create_atom(spacing * 2.0, 0.0, 0.0, 3, "A"),
            create_atom(spacing * 3.0, 0.0, 0.0, 4, "A"),
            create_atom(spacing * 4.0, 0.0, 0.0, 5, "A"),
        ];
        structure
    }

    #[test]
    fn test_lddt_self_comparison() {
        let structure = create_linear_structure(3.8);
        let result = calculate_lddt(
            &structure,
            &structure,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
        .unwrap();

        // Self-comparison should have perfect LDDT
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "Self-LDDT should be 1.0, got {}",
            result.score
        );
    }

    #[test]
    fn test_lddt_translation_invariance() {
        let reference = create_linear_structure(3.8);
        let mut model = create_linear_structure(3.8);

        // Translate the model
        for atom in &mut model.atoms {
            atom.x += 100.0;
            atom.y += 50.0;
            atom.z += 25.0;
        }

        let result = calculate_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
        .unwrap();

        // LDDT should be invariant to translation
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "LDDT should be translation invariant, got {}",
            result.score
        );
    }

    #[test]
    fn test_lddt_rotation_invariance() {
        let reference = create_linear_structure(3.8);
        let mut model = create_linear_structure(3.8);

        // Rotate 90 degrees around z-axis
        for atom in &mut model.atoms {
            let x = atom.x;
            let y = atom.y;
            atom.x = -y;
            atom.y = x;
        }

        let result = calculate_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
        .unwrap();

        // LDDT should be invariant to rotation
        assert!(
            (result.score - 1.0).abs() < 1e-10,
            "LDDT should be rotation invariant, got {}",
            result.score
        );
    }

    #[test]
    fn test_lddt_perturbed_structure() {
        let reference = create_linear_structure(3.8);
        let mut model = create_linear_structure(3.8);

        // Perturb one atom significantly
        model.atoms[2].y += 5.0; // Move middle atom by 5 Angstroms

        let result = calculate_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
        .unwrap();

        // LDDT should be less than 1.0 for perturbed structure
        assert!(
            result.score < 1.0,
            "LDDT should be < 1.0 for perturbed structure, got {}",
            result.score
        );
        assert!(
            result.score > 0.0,
            "LDDT should be > 0.0 for perturbed structure"
        );
    }

    #[test]
    fn test_lddt_custom_options() {
        let reference = create_linear_structure(3.8);
        let mut model = create_linear_structure(3.8);

        // Larger perturbation that will fail strict thresholds but pass lenient ones
        model.atoms[2].y += 1.5; // 1.5 Angstrom perturbation

        let options_lenient = LddtOptions::default().with_thresholds(vec![2.0, 4.0]);
        let options_strict = LddtOptions::default().with_thresholds(vec![0.5, 1.0]);

        let result_lenient =
            calculate_lddt(&model, &reference, AtomSelection::CaOnly, options_lenient).unwrap();
        let result_strict =
            calculate_lddt(&model, &reference, AtomSelection::CaOnly, options_strict).unwrap();

        // Stricter thresholds should give lower or equal score
        assert!(
            result_strict.score <= result_lenient.score,
            "Stricter thresholds should give lower or equal LDDT: {} vs {}",
            result_strict.score,
            result_lenient.score
        );
    }

    #[test]
    fn test_lddt_inclusion_radius() {
        let reference = create_linear_structure(10.0); // Larger spacing
        let mut model = create_linear_structure(10.0);

        // Perturb distant atom
        model.atoms[4].y += 2.0;

        let options_small_radius = LddtOptions::default().with_inclusion_radius(5.0);
        let options_large_radius = LddtOptions::default().with_inclusion_radius(50.0);

        let result_small = calculate_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            options_small_radius,
        )
        .unwrap();
        let result_large = calculate_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            options_large_radius,
        )
        .unwrap();

        // Smaller radius should have fewer pairs and potentially different score
        assert!(result_small.num_pairs <= result_large.num_pairs);
    }

    #[test]
    fn test_per_residue_lddt() {
        let reference = create_linear_structure(3.8);
        let mut model = create_linear_structure(3.8);

        // Perturb middle residue
        model.atoms[2].y += 5.0;

        let per_res = per_residue_lddt(
            &model,
            &reference,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        )
        .unwrap();

        assert_eq!(per_res.len(), 5, "Should have 5 residues");

        // Find the perturbed residue (residue 3)
        let perturbed = per_res.iter().find(|r| r.residue_id.1 == 3).unwrap();

        // The perturbed residue should have lower LDDT
        let others: Vec<_> = per_res.iter().filter(|r| r.residue_id.1 != 3).collect();
        let avg_others = others.iter().map(|r| r.score).sum::<f64>() / others.len() as f64;

        assert!(
            perturbed.score < avg_others,
            "Perturbed residue should have lower LDDT: {} vs {}",
            perturbed.score,
            avg_others
        );
    }

    #[test]
    fn test_lddt_empty_structure() {
        let structure = PdbStructure::new();

        let result = calculate_lddt(
            &structure,
            &structure,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        );

        assert!(matches!(result, Err(PdbError::NoAtomsSelected(_))));
    }

    #[test]
    fn test_lddt_mismatched_structures() {
        let structure1 = create_linear_structure(3.8);
        let mut structure2 = create_linear_structure(3.8);
        structure2.atoms.pop(); // Remove one atom

        let result = calculate_lddt(
            &structure1,
            &structure2,
            AtomSelection::CaOnly,
            LddtOptions::default(),
        );

        assert!(matches!(result, Err(PdbError::AtomCountMismatch { .. })));
    }

    #[test]
    fn test_lddt_options_builder() {
        let options = LddtOptions::new()
            .with_inclusion_radius(10.0)
            .with_thresholds(vec![0.25, 0.5, 1.0]);

        assert_eq!(options.inclusion_radius, 10.0);
        assert_eq!(options.thresholds, vec![0.25, 0.5, 1.0]);
    }

    #[test]
    fn test_distance_function() {
        let p1 = (0.0, 0.0, 0.0);
        let p2 = (3.0, 4.0, 0.0);
        let d = distance(&p1, &p2);
        assert!((d - 5.0).abs() < 1e-10, "Distance should be 5.0, got {}", d);
    }
}
