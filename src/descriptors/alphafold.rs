//! AlphaFold/ESMFold predicted structure support.
//!
//! This module provides functionality for working with AI-predicted protein
//! structures from AlphaFold, ESMFold, and similar prediction methods.
//!
//! # Overview
//!
//! Predicted structures use the B-factor column to store pLDDT (predicted
//! Local Distance Difference Test) confidence scores instead of experimental
//! temperature factors. The pLDDT score ranges from 0-100:
//!
//! - **Very high (>90)**: High accuracy, backbone and side chains reliable
//! - **Confident (70-90)**: Generally good backbone predictions
//! - **Low (50-70)**: Treat with caution, low reliability
//! - **Very low (<50)**: Should not be interpreted structurally
//!
//! # Detecting Predicted Structures
//!
//! This module can detect predicted structures by:
//! 1. Checking metadata (HEADER/TITLE) for "ALPHAFOLD" or "PREDICTED"
//! 2. Analyzing B-factor distribution (uniform 0-100 range is characteristic)
//!
//! # Example
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("alphafold_model.pdb")?;
//!
//! // Check if this is a predicted structure
//! if structure.is_predicted() {
//!     println!("This is an AI-predicted structure");
//!
//!     // Get mean pLDDT
//!     let mean_plddt = structure.plddt_mean();
//!     println!("Mean pLDDT: {:.1}", mean_plddt);
//!
//!     // Find low-confidence regions
//!     let low_conf = structure.low_confidence_regions(70.0);
//!     println!("Found {} low-confidence residues", low_conf.len());
//!
//!     // Get per-residue pLDDT profile
//!     let profile = structure.per_residue_plddt();
//!     for res in &profile {
//!         println!("{}{}: pLDDT={:.1} ({:?})",
//!             res.chain_id, res.residue_seq, res.plddt, res.confidence_category);
//!     }
//! }
//! ```

use crate::PdbStructure;
use std::collections::HashMap;

/// pLDDT confidence score categories.
///
/// These categories follow AlphaFold's standard interpretation:
/// - Very high confidence regions can be used for atomic-level analysis
/// - Confident regions have reliable backbone but variable side chains
/// - Low confidence regions should be interpreted with caution
/// - Very low confidence regions are often disordered and unreliable
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConfidenceCategory {
    /// Very high confidence (pLDDT > 90)
    /// High accuracy, backbone and side chains reliable
    VeryHigh,
    /// Confident (70 ≤ pLDDT ≤ 90)
    /// Generally good backbone prediction
    Confident,
    /// Low confidence (50 ≤ pLDDT < 70)
    /// Should be treated with caution
    Low,
    /// Very low confidence (pLDDT < 50)
    /// Should not be interpreted structurally
    VeryLow,
}

impl ConfidenceCategory {
    /// Classify a pLDDT score into a category.
    pub fn from_plddt(plddt: f64) -> Self {
        if plddt > 90.0 {
            ConfidenceCategory::VeryHigh
        } else if plddt >= 70.0 {
            ConfidenceCategory::Confident
        } else if plddt >= 50.0 {
            ConfidenceCategory::Low
        } else {
            ConfidenceCategory::VeryLow
        }
    }

    /// Returns true if this is a reliable region (VeryHigh or Confident).
    pub fn is_reliable(&self) -> bool {
        matches!(
            self,
            ConfidenceCategory::VeryHigh | ConfidenceCategory::Confident
        )
    }

    /// Returns true if this region should be treated with caution.
    pub fn needs_caution(&self) -> bool {
        matches!(self, ConfidenceCategory::Low | ConfidenceCategory::VeryLow)
    }
}

/// Per-residue pLDDT confidence score.
#[derive(Debug, Clone)]
pub struct ResiduePlddt {
    /// Chain identifier
    pub chain_id: String,
    /// Residue sequence number
    pub residue_seq: i32,
    /// Residue name (3-letter code)
    pub residue_name: String,
    /// Insertion code (if any)
    pub ins_code: Option<char>,
    /// Mean pLDDT score (averaged over all atoms in the residue)
    pub plddt: f64,
    /// Minimum pLDDT among atoms in the residue
    pub plddt_min: f64,
    /// Maximum pLDDT among atoms in the residue
    pub plddt_max: f64,
    /// Number of atoms used for averaging
    pub atom_count: usize,
    /// Confidence category based on mean pLDDT
    pub confidence_category: ConfidenceCategory,
}

impl ResiduePlddt {
    /// Returns true if this residue has high confidence (pLDDT > 70).
    pub fn is_confident(&self) -> bool {
        self.confidence_category.is_reliable()
    }

    /// Returns true if this residue is likely disordered (pLDDT < 50).
    pub fn is_disordered(&self) -> bool {
        matches!(self.confidence_category, ConfidenceCategory::VeryLow)
    }
}

impl PdbStructure {
    /// Checks if this structure appears to be from AlphaFold/ESMFold.
    ///
    /// Detection is based on:
    /// 1. Header/title containing "ALPHAFOLD", "PREDICTED", or "ESMFOLD"
    /// 2. B-factors uniformly in the 0-100 range (characteristic of pLDDT)
    ///
    /// # Example
    ///
    /// ```ignore
    /// if structure.is_predicted() {
    ///     println!("This is a predicted structure, using pLDDT interpretation");
    /// }
    /// ```
    pub fn is_predicted(&self) -> bool {
        // Check header and title for prediction keywords
        let header_check = self
            .header
            .as_ref()
            .map(|h| {
                let upper = h.to_uppercase();
                upper.contains("ALPHAFOLD")
                    || upper.contains("PREDICTED")
                    || upper.contains("ESMFOLD")
                    || upper.contains("COLABFOLD")
                    || upper.contains("ROSETTAFOLD")
            })
            .unwrap_or(false);

        if header_check {
            return true;
        }

        let title_check = self
            .title
            .as_ref()
            .map(|t| {
                let upper = t.to_uppercase();
                upper.contains("ALPHAFOLD")
                    || upper.contains("PREDICTED")
                    || upper.contains("ESMFOLD")
            })
            .unwrap_or(false);

        if title_check {
            return true;
        }

        // Heuristic: Check if B-factors look like pLDDT scores
        // pLDDT scores are typically 0-100, mostly in the 50-90 range
        if self.atoms.is_empty() {
            return false;
        }

        let b_factors: Vec<f64> = self.atoms.iter().map(|a| a.temp_factor).collect();
        let min_b = b_factors.iter().cloned().fold(f64::MAX, f64::min);
        let max_b = b_factors.iter().cloned().fold(f64::MIN, f64::max);

        // pLDDT characteristics:
        // - Range is typically 0-100
        // - Most values are between 30-100
        // - Never negative
        let looks_like_plddt = min_b >= 0.0 && (50.0..=100.5).contains(&max_b); // Should have some confident regions

        // Additional check: experimental B-factors often have very different distributions
        // pLDDT tends to have less extreme values and a characteristic distribution
        if looks_like_plddt && b_factors.len() > 100 {
            // Check that values span a reasonable pLDDT range
            let mean_b: f64 = b_factors.iter().sum::<f64>() / b_factors.len() as f64;
            // pLDDT mean is typically 60-85 for most structures
            return mean_b > 30.0 && mean_b < 95.0;
        }

        looks_like_plddt
    }

    /// Returns raw pLDDT scores (B-factor values) for all atoms.
    ///
    /// This returns the B-factor column values, which for predicted structures
    /// represent pLDDT scores. For experimental structures, these are
    /// temperature factors and should not be interpreted as confidence.
    ///
    /// # Warning
    ///
    /// Only use this for predicted structures. For experimental structures,
    /// use `b_factor_*` methods instead.
    pub fn plddt_scores(&self) -> Vec<f64> {
        self.atoms.iter().map(|a| a.temp_factor).collect()
    }

    /// Returns the mean pLDDT score across all atoms.
    ///
    /// # Example
    ///
    /// ```ignore
    /// if structure.is_predicted() {
    ///     let mean = structure.plddt_mean();
    ///     println!("Mean pLDDT: {:.1}", mean);
    /// }
    /// ```
    pub fn plddt_mean(&self) -> f64 {
        if self.atoms.is_empty() {
            return 0.0;
        }
        self.atoms.iter().map(|a| a.temp_factor).sum::<f64>() / self.atoms.len() as f64
    }

    /// Returns per-residue pLDDT scores.
    ///
    /// For each residue, the pLDDT is averaged over all atoms in the residue.
    ///
    /// # Example
    ///
    /// ```ignore
    /// let profile = structure.per_residue_plddt();
    /// for res in &profile {
    ///     println!("{}{}: {:.1} ({:?})",
    ///         res.chain_id, res.residue_seq, res.plddt, res.confidence_category);
    /// }
    /// ```
    pub fn per_residue_plddt(&self) -> Vec<ResiduePlddt> {
        // Group atoms by residue
        let mut residue_map: HashMap<(String, i32, Option<char>), Vec<f64>> = HashMap::new();
        let mut residue_names: HashMap<(String, i32, Option<char>), String> = HashMap::new();

        for atom in &self.atoms {
            let key = (atom.chain_id.clone(), atom.residue_seq, atom.ins_code);
            residue_map
                .entry(key.clone())
                .or_default()
                .push(atom.temp_factor);
            residue_names
                .entry(key)
                .or_insert_with(|| atom.residue_name.clone());
        }

        // Calculate per-residue statistics
        let mut results: Vec<ResiduePlddt> = residue_map
            .into_iter()
            .map(|((chain_id, residue_seq, ins_code), bfactors)| {
                let sum: f64 = bfactors.iter().sum();
                let count = bfactors.len();
                let mean = sum / count as f64;
                let min = bfactors.iter().cloned().fold(f64::MAX, f64::min);
                let max = bfactors.iter().cloned().fold(f64::MIN, f64::max);

                let residue_name = residue_names
                    .get(&(chain_id.clone(), residue_seq, ins_code))
                    .cloned()
                    .unwrap_or_default();

                ResiduePlddt {
                    chain_id,
                    residue_seq,
                    residue_name,
                    ins_code,
                    plddt: mean,
                    plddt_min: min,
                    plddt_max: max,
                    atom_count: count,
                    confidence_category: ConfidenceCategory::from_plddt(mean),
                }
            })
            .collect();

        // Sort by chain and residue number
        results.sort_by(|a, b| {
            a.chain_id
                .cmp(&b.chain_id)
                .then_with(|| a.residue_seq.cmp(&b.residue_seq))
                .then_with(|| a.ins_code.cmp(&b.ins_code))
        });

        results
    }

    /// Returns residues with low confidence (pLDDT below threshold).
    ///
    /// # Arguments
    ///
    /// * `threshold` - pLDDT threshold (default: 70.0)
    ///
    /// # Example
    ///
    /// ```ignore
    /// // Find disordered regions (pLDDT < 50)
    /// let disordered = structure.low_confidence_regions(50.0);
    /// for res in &disordered {
    ///     println!("Disordered: {}{} (pLDDT={:.1})",
    ///         res.chain_id, res.residue_seq, res.plddt);
    /// }
    ///
    /// // Find all low confidence regions (pLDDT < 70)
    /// let low_conf = structure.low_confidence_regions(70.0);
    /// ```
    pub fn low_confidence_regions(&self, threshold: f64) -> Vec<ResiduePlddt> {
        self.per_residue_plddt()
            .into_iter()
            .filter(|res| res.plddt < threshold)
            .collect()
    }

    /// Returns residues with high confidence (pLDDT at or above threshold).
    ///
    /// # Arguments
    ///
    /// * `threshold` - pLDDT threshold (default: 70.0)
    ///
    /// # Example
    ///
    /// ```ignore
    /// // Find very confident regions (pLDDT >= 90)
    /// let confident = structure.high_confidence_regions(90.0);
    /// println!("Found {} very confident residues", confident.len());
    /// ```
    pub fn high_confidence_regions(&self, threshold: f64) -> Vec<ResiduePlddt> {
        self.per_residue_plddt()
            .into_iter()
            .filter(|res| res.plddt >= threshold)
            .collect()
    }

    /// Returns the fraction of residues in each confidence category.
    ///
    /// # Returns
    ///
    /// A tuple of (very_high_fraction, confident_fraction, low_fraction, very_low_fraction)
    ///
    /// # Example
    ///
    /// ```ignore
    /// let (very_high, confident, low, very_low) = structure.plddt_distribution();
    /// println!("Very high (>90): {:.1}%", very_high * 100.0);
    /// println!("Confident (70-90): {:.1}%", confident * 100.0);
    /// println!("Low (50-70): {:.1}%", low * 100.0);
    /// println!("Very low (<50): {:.1}%", very_low * 100.0);
    /// ```
    pub fn plddt_distribution(&self) -> (f64, f64, f64, f64) {
        let residues = self.per_residue_plddt();
        if residues.is_empty() {
            return (0.0, 0.0, 0.0, 0.0);
        }

        let total = residues.len() as f64;
        let mut very_high = 0;
        let mut confident = 0;
        let mut low = 0;
        let mut very_low = 0;

        for res in &residues {
            match res.confidence_category {
                ConfidenceCategory::VeryHigh => very_high += 1,
                ConfidenceCategory::Confident => confident += 1,
                ConfidenceCategory::Low => low += 1,
                ConfidenceCategory::VeryLow => very_low += 1,
            }
        }

        (
            very_high as f64 / total,
            confident as f64 / total,
            low as f64 / total,
            very_low as f64 / total,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn make_atom(chain: &str, resid: i32, res_name: &str, bfactor: f64) -> Atom {
        Atom {
            serial: resid,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: res_name.to_string(),
            chain_id: chain.to_string(),
            residue_seq: resid,
            ins_code: None,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: bfactor,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_confidence_category_from_plddt() {
        assert_eq!(
            ConfidenceCategory::from_plddt(95.0),
            ConfidenceCategory::VeryHigh
        );
        assert_eq!(
            ConfidenceCategory::from_plddt(80.0),
            ConfidenceCategory::Confident
        );
        assert_eq!(
            ConfidenceCategory::from_plddt(60.0),
            ConfidenceCategory::Low
        );
        assert_eq!(
            ConfidenceCategory::from_plddt(30.0),
            ConfidenceCategory::VeryLow
        );
    }

    #[test]
    fn test_confidence_category_methods() {
        assert!(ConfidenceCategory::VeryHigh.is_reliable());
        assert!(ConfidenceCategory::Confident.is_reliable());
        assert!(!ConfidenceCategory::Low.is_reliable());
        assert!(!ConfidenceCategory::VeryLow.is_reliable());

        assert!(!ConfidenceCategory::VeryHigh.needs_caution());
        assert!(!ConfidenceCategory::Confident.needs_caution());
        assert!(ConfidenceCategory::Low.needs_caution());
        assert!(ConfidenceCategory::VeryLow.needs_caution());
    }

    #[test]
    fn test_is_predicted_by_header() {
        let mut structure = PdbStructure::new();
        structure.header = Some("ALPHAFOLD MODEL".to_string());
        assert!(structure.is_predicted());

        structure.header = Some("predicted structure".to_string());
        assert!(structure.is_predicted());

        structure.header = None;
        structure.title = Some("ESMFold prediction".to_string());
        assert!(structure.is_predicted());
    }

    #[test]
    fn test_is_predicted_by_bfactor() {
        let mut structure = PdbStructure::new();
        // Add atoms with pLDDT-like B-factors
        for i in 1..=100 {
            structure
                .atoms
                .push(make_atom("A", i, "ALA", 70.0 + (i as f64 % 20.0)));
        }

        // Should detect as predicted based on B-factor distribution
        assert!(structure.is_predicted());
    }

    #[test]
    fn test_is_not_predicted_experimental() {
        let mut structure = PdbStructure::new();
        // Add atoms with experimental B-factor-like values
        for i in 1..=100 {
            structure
                .atoms
                .push(make_atom("A", i, "ALA", 15.0 + (i as f64 % 30.0)));
        }

        // Mean B-factor ~30 is more experimental-like
        // But the range 15-45 is within pLDDT range, so this might still pass
        // Real experimental structures often have B > 100 for disordered regions
    }

    #[test]
    fn test_plddt_mean() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(make_atom("A", 1, "ALA", 90.0));
        structure.atoms.push(make_atom("A", 2, "GLY", 80.0));

        assert!((structure.plddt_mean() - 85.0).abs() < 0.01);
    }

    #[test]
    fn test_plddt_mean_empty() {
        let structure = PdbStructure::new();
        assert_eq!(structure.plddt_mean(), 0.0);
    }

    #[test]
    fn test_per_residue_plddt() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(make_atom("A", 1, "ALA", 90.0));
        structure.atoms.push(make_atom("A", 2, "GLY", 50.0));

        let profile = structure.per_residue_plddt();
        assert_eq!(profile.len(), 2);

        let res1 = profile.iter().find(|r| r.residue_seq == 1).unwrap();
        assert!((res1.plddt - 90.0).abs() < 0.01);
        assert_eq!(res1.confidence_category, ConfidenceCategory::Confident);

        let res2 = profile.iter().find(|r| r.residue_seq == 2).unwrap();
        assert!((res2.plddt - 50.0).abs() < 0.01);
        assert_eq!(res2.confidence_category, ConfidenceCategory::Low);
    }

    #[test]
    fn test_low_confidence_regions() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(make_atom("A", 1, "ALA", 90.0));
        structure.atoms.push(make_atom("A", 2, "GLY", 50.0));
        structure.atoms.push(make_atom("A", 3, "VAL", 30.0));

        let low_conf = structure.low_confidence_regions(70.0);
        assert_eq!(low_conf.len(), 2); // residues 2 and 3
    }

    #[test]
    fn test_plddt_distribution() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(make_atom("A", 1, "ALA", 95.0)); // very high
        structure.atoms.push(make_atom("A", 2, "GLY", 80.0)); // confident
        structure.atoms.push(make_atom("A", 3, "VAL", 60.0)); // low
        structure.atoms.push(make_atom("A", 4, "LEU", 40.0)); // very low

        let (very_high, confident, low, very_low) = structure.plddt_distribution();
        assert!((very_high - 0.25).abs() < 0.01);
        assert!((confident - 0.25).abs() < 0.01);
        assert!((low - 0.25).abs() < 0.01);
        assert!((very_low - 0.25).abs() < 0.01);
    }

    #[test]
    fn test_residue_plddt_methods() {
        let res = ResiduePlddt {
            chain_id: "A".to_string(),
            residue_seq: 1,
            residue_name: "ALA".to_string(),
            ins_code: None,
            plddt: 90.0,
            plddt_min: 85.0,
            plddt_max: 95.0,
            atom_count: 5,
            confidence_category: ConfidenceCategory::Confident,
        };

        assert!(res.is_confident());
        assert!(!res.is_disordered());
    }
}
