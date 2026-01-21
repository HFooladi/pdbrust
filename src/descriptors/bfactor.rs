//! B-factor (temperature factor) analysis for protein structures.
//!
//! B-factors represent the displacement of atoms from their mean positions
//! and are commonly used to identify flexible regions in protein structures.
//!
//! # Examples
//!
//! ```ignore
//! use pdbrust::PdbStructure;
//!
//! let structure = PdbStructure::from_file("protein.pdb")?;
//!
//! // Basic statistics
//! let mean = structure.b_factor_mean();
//! let mean_ca = structure.b_factor_mean_ca();
//!
//! // Per-residue profile
//! let profile = structure.b_factor_profile();
//! for res in &profile {
//!     println!("{}{}: {:.2} Å²", res.chain_id, res.residue_seq, res.b_factor_mean);
//! }
//!
//! // Identify flexible regions
//! let flexible = structure.flexible_residues(50.0);  // B > 50 Å²
//! ```

use super::ResidueBFactor;
use crate::core::PdbStructure;
use std::collections::HashMap;

impl PdbStructure {
    /// Calculate the mean B-factor across all atoms.
    ///
    /// # Returns
    ///
    /// The arithmetic mean of B-factors (temp_factor) in Å².
    /// Returns 0.0 for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let mean_b = structure.b_factor_mean();
    /// println!("Mean B-factor: {:.2} Å²", mean_b);
    /// ```
    pub fn b_factor_mean(&self) -> f64 {
        if self.atoms.is_empty() {
            return 0.0;
        }

        let sum: f64 = self.atoms.iter().map(|a| a.temp_factor).sum();
        sum / self.atoms.len() as f64
    }

    /// Calculate the mean B-factor for CA atoms only.
    ///
    /// This is often more representative of residue-level mobility
    /// as it avoids the influence of side chain dynamics.
    ///
    /// # Returns
    ///
    /// The arithmetic mean of B-factors for CA atoms in Å².
    /// Returns 0.0 if no CA atoms are present.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let mean_ca_b = structure.b_factor_mean_ca();
    /// println!("Mean CA B-factor: {:.2} Å²", mean_ca_b);
    /// ```
    pub fn b_factor_mean_ca(&self) -> f64 {
        let ca_bfactors: Vec<f64> = self
            .atoms
            .iter()
            .filter(|a| a.name.trim() == "CA")
            .map(|a| a.temp_factor)
            .collect();

        if ca_bfactors.is_empty() {
            return 0.0;
        }

        ca_bfactors.iter().sum::<f64>() / ca_bfactors.len() as f64
    }

    /// Get the minimum B-factor value in the structure.
    ///
    /// # Returns
    ///
    /// The minimum B-factor value in Å².
    /// Returns 0.0 for empty structures.
    pub fn b_factor_min(&self) -> f64 {
        if self.atoms.is_empty() {
            return 0.0;
        }
        self.atoms
            .iter()
            .map(|a| a.temp_factor)
            .fold(f64::INFINITY, f64::min)
    }

    /// Get the maximum B-factor value in the structure.
    ///
    /// # Returns
    ///
    /// The maximum B-factor value in Å².
    /// Returns 0.0 for empty structures.
    pub fn b_factor_max(&self) -> f64 {
        if self.atoms.is_empty() {
            return 0.0;
        }
        self.atoms
            .iter()
            .map(|a| a.temp_factor)
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Calculate the standard deviation of B-factors.
    ///
    /// # Returns
    ///
    /// The standard deviation of B-factors in Å².
    /// Returns 0.0 for structures with fewer than 2 atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let std_b = structure.b_factor_std();
    /// println!("B-factor std dev: {:.2} Å²", std_b);
    /// ```
    pub fn b_factor_std(&self) -> f64 {
        if self.atoms.len() < 2 {
            return 0.0;
        }

        let mean = self.b_factor_mean();
        let variance: f64 = self
            .atoms
            .iter()
            .map(|a| (a.temp_factor - mean).powi(2))
            .sum::<f64>()
            / self.atoms.len() as f64;

        variance.sqrt()
    }

    /// Compute per-residue B-factor statistics.
    ///
    /// For each unique residue (identified by chain_id, residue_seq, ins_code),
    /// computes the mean, min, and max B-factors across all atoms in that residue.
    ///
    /// # Returns
    ///
    /// A vector of `ResidueBFactor` structs, sorted by chain_id and residue_seq.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let profile = structure.b_factor_profile();
    ///
    /// for res in &profile {
    ///     println!("{} {}{}: mean={:.2}, min={:.2}, max={:.2}",
    ///         res.residue_name, res.chain_id, res.residue_seq,
    ///         res.b_factor_mean, res.b_factor_min, res.b_factor_max);
    /// }
    /// ```
    pub fn b_factor_profile(&self) -> Vec<ResidueBFactor> {
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

        // Compute statistics for each residue
        let mut profile: Vec<ResidueBFactor> = residue_map
            .into_iter()
            .map(|((chain_id, residue_seq, ins_code), bfactors)| {
                let n = bfactors.len();
                let sum: f64 = bfactors.iter().sum();
                let min = bfactors.iter().cloned().fold(f64::INFINITY, f64::min);
                let max = bfactors.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

                ResidueBFactor {
                    chain_id: chain_id.clone(),
                    residue_seq,
                    ins_code,
                    residue_name: residue_names
                        .get(&(chain_id, residue_seq, ins_code))
                        .cloned()
                        .unwrap_or_default(),
                    b_factor_mean: sum / n as f64,
                    b_factor_min: min,
                    b_factor_max: max,
                    atom_count: n,
                }
            })
            .collect();

        // Sort by chain_id, then by residue_seq
        profile.sort_by(|a, b| {
            a.chain_id
                .cmp(&b.chain_id)
                .then_with(|| a.residue_seq.cmp(&b.residue_seq))
        });

        profile
    }

    /// Identify residues with B-factors above a threshold (flexible regions).
    ///
    /// Returns residues where the mean B-factor exceeds the specified threshold.
    /// High B-factors typically indicate mobile or disordered regions.
    ///
    /// # Arguments
    ///
    /// * `threshold` - B-factor cutoff in Å². Residues with mean B > threshold
    ///   are returned. Common values: 40-60 Å² for identifying flexible loops.
    ///
    /// # Returns
    ///
    /// A vector of `ResidueBFactor` for residues above the threshold,
    /// sorted by chain_id and residue_seq.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let flexible = structure.flexible_residues(50.0);
    ///
    /// println!("Flexible regions (B > 50 Å²):");
    /// for res in &flexible {
    ///     println!("  {}{}: {:.2} Å²", res.chain_id, res.residue_seq, res.b_factor_mean);
    /// }
    /// ```
    pub fn flexible_residues(&self, threshold: f64) -> Vec<ResidueBFactor> {
        self.b_factor_profile()
            .into_iter()
            .filter(|res| res.b_factor_mean > threshold)
            .collect()
    }

    /// Identify residues with B-factors below a threshold (rigid regions).
    ///
    /// Returns residues where the mean B-factor is below the specified threshold.
    /// Low B-factors typically indicate well-ordered, rigid regions.
    ///
    /// # Arguments
    ///
    /// * `threshold` - B-factor cutoff in Å². Residues with mean B < threshold
    ///   are returned. Common values: 15-25 Å² for identifying rigid core.
    ///
    /// # Returns
    ///
    /// A vector of `ResidueBFactor` for residues below the threshold,
    /// sorted by chain_id and residue_seq.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let rigid = structure.rigid_residues(20.0);
    ///
    /// println!("Rigid core residues (B < 20 Å²):");
    /// for res in &rigid {
    ///     println!("  {}{}: {:.2} Å²", res.chain_id, res.residue_seq, res.b_factor_mean);
    /// }
    /// ```
    pub fn rigid_residues(&self, threshold: f64) -> Vec<ResidueBFactor> {
        self.b_factor_profile()
            .into_iter()
            .filter(|res| res.b_factor_mean < threshold)
            .collect()
    }

    /// Create a new structure with Z-score normalized B-factors.
    ///
    /// Normalization transforms B-factors to have mean 0 and standard deviation 1,
    /// enabling comparison of flexibility between different structures.
    ///
    /// # Formula
    ///
    /// B_normalized = (B - mean) / std
    ///
    /// # Returns
    ///
    /// A new `PdbStructure` with normalized B-factors. If the structure has
    /// fewer than 2 atoms or zero standard deviation, B-factors are set to 0.0.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let normalized = structure.normalize_b_factors();
    ///
    /// // Normalized structure has mean ~0 and std ~1
    /// println!("Normalized mean: {:.4}", normalized.b_factor_mean());
    /// println!("Normalized std:  {:.4}", normalized.b_factor_std());
    /// ```
    pub fn normalize_b_factors(&self) -> PdbStructure {
        let mut normalized = self.clone();

        if self.atoms.len() < 2 {
            // Set all to 0.0 for degenerate cases
            for atom in &mut normalized.atoms {
                atom.temp_factor = 0.0;
            }
            return normalized;
        }

        let mean = self.b_factor_mean();
        let std = self.b_factor_std();

        if std < 1e-10 {
            // All B-factors are identical, normalize to 0.0
            for atom in &mut normalized.atoms {
                atom.temp_factor = 0.0;
            }
        } else {
            for atom in &mut normalized.atoms {
                atom.temp_factor = (atom.temp_factor - mean) / std;
            }
        }

        normalized
    }

    /// Get the percentile rank of an atom's B-factor within the structure.
    ///
    /// # Arguments
    ///
    /// * `atom_serial` - The serial number of the atom to query.
    ///
    /// # Returns
    ///
    /// The percentile (0.0 to 100.0) of the atom's B-factor, or None
    /// if the atom is not found.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// if let Some(percentile) = structure.b_factor_percentile(42) {
    ///     println!("Atom 42 B-factor is at the {:.1}th percentile", percentile);
    /// }
    /// ```
    pub fn b_factor_percentile(&self, atom_serial: i32) -> Option<f64> {
        let target_atom = self.atoms.iter().find(|a| a.serial == atom_serial)?;
        let target_b = target_atom.temp_factor;

        if self.atoms.is_empty() {
            return Some(0.0);
        }

        let count_below = self
            .atoms
            .iter()
            .filter(|a| a.temp_factor < target_b)
            .count();

        Some(100.0 * count_below as f64 / self.atoms.len() as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        temp_factor: f64,
    ) -> Atom {
        Atom {
            serial,
            name: name.to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor,
            element: "C".to_string(),
        }
    }

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", 1, 10.0),
            create_atom(2, " CB ", "ALA", "A", 1, 15.0),
            create_atom(3, " CA ", "GLY", "A", 2, 20.0),
            create_atom(4, " CA ", "VAL", "A", 3, 30.0),
            create_atom(5, " CB ", "VAL", "A", 3, 40.0),
            create_atom(6, " CG1", "VAL", "A", 3, 50.0),
        ];
        structure
    }

    #[test]
    fn test_b_factor_mean() {
        let structure = create_test_structure();
        // Mean of [10, 15, 20, 30, 40, 50] = 165 / 6 = 27.5
        let mean = structure.b_factor_mean();
        assert!((mean - 27.5).abs() < 1e-10);
    }

    #[test]
    fn test_b_factor_mean_empty() {
        let structure = PdbStructure::new();
        assert_eq!(structure.b_factor_mean(), 0.0);
    }

    #[test]
    fn test_b_factor_mean_ca() {
        let structure = create_test_structure();
        // Mean of CA only [10, 20, 30] = 60 / 3 = 20.0
        let mean_ca = structure.b_factor_mean_ca();
        assert!((mean_ca - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_b_factor_mean_ca_no_ca() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CB ", "ALA", "A", 1, 10.0),
            create_atom(2, " CB ", "GLY", "A", 2, 20.0),
        ];
        assert_eq!(structure.b_factor_mean_ca(), 0.0);
    }

    #[test]
    fn test_b_factor_min_max() {
        let structure = create_test_structure();
        assert!((structure.b_factor_min() - 10.0).abs() < 1e-10);
        assert!((structure.b_factor_max() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_b_factor_min_max_empty() {
        let structure = PdbStructure::new();
        assert_eq!(structure.b_factor_min(), 0.0);
        assert_eq!(structure.b_factor_max(), 0.0);
    }

    #[test]
    fn test_b_factor_std() {
        let structure = create_test_structure();
        // Variance of [10, 15, 20, 30, 40, 50] with mean 27.5
        // = (306.25 + 156.25 + 56.25 + 6.25 + 156.25 + 506.25) / 6
        // = 1187.5 / 6 = 197.916...
        // Std = sqrt(197.916...) ≈ 14.068
        let std = structure.b_factor_std();
        assert!((std - 14.068).abs() < 0.01);
    }

    #[test]
    fn test_b_factor_std_single_atom() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![create_atom(1, " CA ", "ALA", "A", 1, 10.0)];
        assert_eq!(structure.b_factor_std(), 0.0);
    }

    #[test]
    fn test_b_factor_profile() {
        let structure = create_test_structure();
        let profile = structure.b_factor_profile();

        assert_eq!(profile.len(), 3);

        // First residue: ALA A1, atoms with B=[10, 15], mean=12.5
        assert_eq!(profile[0].chain_id, "A");
        assert_eq!(profile[0].residue_seq, 1);
        assert_eq!(profile[0].residue_name, "ALA");
        assert!((profile[0].b_factor_mean - 12.5).abs() < 1e-10);
        assert!((profile[0].b_factor_min - 10.0).abs() < 1e-10);
        assert!((profile[0].b_factor_max - 15.0).abs() < 1e-10);
        assert_eq!(profile[0].atom_count, 2);

        // Second residue: GLY A2, atoms with B=[20], mean=20.0
        assert_eq!(profile[1].chain_id, "A");
        assert_eq!(profile[1].residue_seq, 2);
        assert_eq!(profile[1].residue_name, "GLY");
        assert!((profile[1].b_factor_mean - 20.0).abs() < 1e-10);
        assert_eq!(profile[1].atom_count, 1);

        // Third residue: VAL A3, atoms with B=[30, 40, 50], mean=40.0
        assert_eq!(profile[2].chain_id, "A");
        assert_eq!(profile[2].residue_seq, 3);
        assert_eq!(profile[2].residue_name, "VAL");
        assert!((profile[2].b_factor_mean - 40.0).abs() < 1e-10);
        assert!((profile[2].b_factor_min - 30.0).abs() < 1e-10);
        assert!((profile[2].b_factor_max - 50.0).abs() < 1e-10);
        assert_eq!(profile[2].atom_count, 3);
    }

    #[test]
    fn test_flexible_residues() {
        let structure = create_test_structure();
        let flexible = structure.flexible_residues(25.0);

        // Only VAL A3 has mean B > 25 (mean=40.0)
        assert_eq!(flexible.len(), 1);
        assert_eq!(flexible[0].residue_seq, 3);
    }

    #[test]
    fn test_rigid_residues() {
        let structure = create_test_structure();
        let rigid = structure.rigid_residues(15.0);

        // Only ALA A1 has mean B < 15 (mean=12.5)
        assert_eq!(rigid.len(), 1);
        assert_eq!(rigid[0].residue_seq, 1);
    }

    #[test]
    fn test_normalize_b_factors() {
        let structure = create_test_structure();
        let normalized = structure.normalize_b_factors();

        // After normalization, mean should be ~0 and std should be ~1
        let norm_mean = normalized.b_factor_mean();
        let norm_std = normalized.b_factor_std();

        assert!(norm_mean.abs() < 1e-10, "Normalized mean should be ~0");
        assert!(
            (norm_std - 1.0).abs() < 1e-10,
            "Normalized std should be ~1"
        );
    }

    #[test]
    fn test_normalize_b_factors_uniform() {
        // All B-factors are identical
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", 1, 20.0),
            create_atom(2, " CA ", "GLY", "A", 2, 20.0),
            create_atom(3, " CA ", "VAL", "A", 3, 20.0),
        ];

        let normalized = structure.normalize_b_factors();

        // All should be 0.0 when std is 0
        for atom in &normalized.atoms {
            assert!((atom.temp_factor - 0.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_b_factor_percentile() {
        let structure = create_test_structure();

        // Atom 1 has B=10 (minimum), should be 0th percentile
        assert!((structure.b_factor_percentile(1).unwrap() - 0.0).abs() < 1e-10);

        // Atom 6 has B=50 (maximum), should be at 5/6 * 100 = 83.33... percentile
        let p = structure.b_factor_percentile(6).unwrap();
        assert!((p - 83.333).abs() < 0.01);
    }

    #[test]
    fn test_b_factor_percentile_not_found() {
        let structure = create_test_structure();
        assert!(structure.b_factor_percentile(999).is_none());
    }

    #[test]
    fn test_multi_chain_profile() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, " CA ", "ALA", "A", 1, 10.0),
            create_atom(2, " CA ", "GLY", "B", 1, 20.0),
            create_atom(3, " CA ", "VAL", "A", 2, 30.0),
        ];

        let profile = structure.b_factor_profile();
        assert_eq!(profile.len(), 3);

        // Should be sorted by chain, then by residue
        assert_eq!(profile[0].chain_id, "A");
        assert_eq!(profile[0].residue_seq, 1);
        assert_eq!(profile[1].chain_id, "A");
        assert_eq!(profile[1].residue_seq, 2);
        assert_eq!(profile[2].chain_id, "B");
        assert_eq!(profile[2].residue_seq, 1);
    }
}
