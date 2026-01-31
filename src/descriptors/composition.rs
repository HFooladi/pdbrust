//! Amino acid composition analysis.
//!
//! Functions for analyzing the amino acid composition of protein structures.

use crate::core::PdbStructure;
use std::collections::HashMap;

/// Hydrophobic amino acid residue names.
///
/// Based on the Kyte-Doolittle hydropathy scale, these residues
/// have positive hydropathy values.
pub const HYDROPHOBIC_RESIDUES: &[&str] = &[
    "ALA", // Alanine
    "VAL", // Valine
    "ILE", // Isoleucine
    "LEU", // Leucine
    "MET", // Methionine
    "PHE", // Phenylalanine
    "TRP", // Tryptophan
    "PRO", // Proline
];

/// Polar (hydrophilic) amino acid residue names.
pub const POLAR_RESIDUES: &[&str] = &[
    "SER", // Serine
    "THR", // Threonine
    "CYS", // Cysteine
    "TYR", // Tyrosine
    "ASN", // Asparagine
    "GLN", // Glutamine
];

/// Charged amino acid residue names.
pub const CHARGED_RESIDUES: &[&str] = &[
    "ASP", // Aspartic acid (negative)
    "GLU", // Glutamic acid (negative)
    "LYS", // Lysine (positive)
    "ARG", // Arginine (positive)
    "HIS", // Histidine (positive at low pH)
];

/// Aromatic amino acid residue names.
pub const AROMATIC_RESIDUES: &[&str] = &[
    "PHE", // Phenylalanine
    "TYR", // Tyrosine
    "TRP", // Tryptophan
];

/// Small amino acid residue names.
pub const SMALL_RESIDUES: &[&str] = &[
    "GLY", // Glycine
    "ALA", // Alanine
    "SER", // Serine
    "PRO", // Proline
];

impl PdbStructure {
    /// Calculate amino acid composition as fractions.
    ///
    /// Returns a HashMap mapping 3-letter amino acid codes to their
    /// frequency in the structure (values between 0.0 and 1.0).
    /// Composition is based on Cα atom counts.
    ///
    /// # Returns
    ///
    /// A HashMap where keys are 3-letter amino acid codes and values
    /// are the fraction of that amino acid in the structure.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let composition = structure.aa_composition();
    ///
    /// // Check fraction of alanine
    /// if let Some(&ala_fraction) = composition.get("ALA") {
    ///     println!("Alanine: {:.1}%", ala_fraction * 100.0);
    /// }
    /// ```
    pub fn aa_composition(&self) -> HashMap<String, f64> {
        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut total = 0usize;

        // Count residues based on CA atoms
        for atom in &self.atoms {
            if atom.name.trim() == "CA" {
                *counts.entry(atom.residue_name.clone()).or_insert(0) += 1;
                total += 1;
            }
        }

        // Convert counts to fractions
        let mut composition = HashMap::new();
        if total > 0 {
            for (residue, count) in counts {
                composition.insert(residue, count as f64 / total as f64);
            }
        }

        composition
    }

    /// Calculate the fraction of glycine residues.
    ///
    /// Glycine is the smallest amino acid and is often found in
    /// flexible regions of proteins.
    ///
    /// # Returns
    ///
    /// The fraction of glycine residues (0.0 to 1.0).
    pub fn glycine_ratio(&self) -> f64 {
        self.aa_composition().get("GLY").copied().unwrap_or(0.0)
    }

    /// Calculate the fraction of hydrophobic residues.
    ///
    /// Hydrophobic residues (ALA, VAL, ILE, LEU, MET, PHE, TRP, PRO)
    /// tend to be buried in the protein core.
    ///
    /// # Returns
    ///
    /// The fraction of hydrophobic residues (0.0 to 1.0).
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let hydrophobic = structure.hydrophobic_ratio();
    ///
    /// println!("Hydrophobic content: {:.1}%", hydrophobic * 100.0);
    /// ```
    pub fn hydrophobic_ratio(&self) -> f64 {
        let composition = self.aa_composition();
        HYDROPHOBIC_RESIDUES
            .iter()
            .filter_map(|&aa| composition.get(aa))
            .sum()
    }

    /// Calculate the fraction of polar (hydrophilic) residues.
    ///
    /// Polar residues (SER, THR, CYS, TYR, ASN, GLN) can form
    /// hydrogen bonds and are often found on the protein surface.
    ///
    /// # Returns
    ///
    /// The fraction of polar residues (0.0 to 1.0).
    pub fn polar_ratio(&self) -> f64 {
        let composition = self.aa_composition();
        POLAR_RESIDUES
            .iter()
            .filter_map(|&aa| composition.get(aa))
            .sum()
    }

    /// Calculate the fraction of charged residues.
    ///
    /// Charged residues (ASP, GLU, LYS, ARG, HIS) are typically
    /// found on the protein surface and are important for
    /// protein-protein interactions.
    ///
    /// # Returns
    ///
    /// The fraction of charged residues (0.0 to 1.0).
    pub fn charged_ratio(&self) -> f64 {
        let composition = self.aa_composition();
        CHARGED_RESIDUES
            .iter()
            .filter_map(|&aa| composition.get(aa))
            .sum()
    }

    /// Calculate the fraction of aromatic residues.
    ///
    /// Aromatic residues (PHE, TYR, TRP) have large side chains
    /// and are often involved in protein-ligand interactions.
    ///
    /// # Returns
    ///
    /// The fraction of aromatic residues (0.0 to 1.0).
    pub fn aromatic_ratio(&self) -> f64 {
        let composition = self.aa_composition();
        AROMATIC_RESIDUES
            .iter()
            .filter_map(|&aa| composition.get(aa))
            .sum()
    }

    /// Calculate the fraction of small residues.
    ///
    /// Small residues (GLY, ALA, SER, PRO) have compact side chains.
    ///
    /// # Returns
    ///
    /// The fraction of small residues (0.0 to 1.0).
    pub fn small_ratio(&self) -> f64 {
        let composition = self.aa_composition();
        SMALL_RESIDUES
            .iter()
            .filter_map(|&aa| composition.get(aa))
            .sum()
    }

    /// Count the number of residues based on Cα atoms.
    ///
    /// This provides a more accurate residue count than checking
    /// unique residue numbers, as it's based on actual Cα atoms present.
    ///
    /// # Returns
    ///
    /// The number of Cα atoms (equivalent to residue count for complete structures).
    pub fn count_ca_residues(&self) -> usize {
        self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .count()
    }

    /// Calculate the missing residue ratio based on sequence gaps.
    ///
    /// Compares the expected number of residues (based on the range
    /// from first to last residue number) to the actual number of
    /// residues present.
    ///
    /// # Returns
    ///
    /// The fraction of missing residues (0.0 = complete, approaching 1.0 = many gaps).
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let missing = structure.missing_residue_ratio();
    ///
    /// if missing > 0.1 {
    ///     println!("Warning: {:.1}% of residues are missing", missing * 100.0);
    /// }
    /// ```
    pub fn missing_residue_ratio(&self) -> f64 {
        // Get all CA atoms and their residue numbers per chain
        let mut chain_residues: HashMap<String, Vec<i32>> = HashMap::new();

        for atom in &self.atoms {
            if atom.name.trim() == "CA" {
                chain_residues
                    .entry(atom.chain_id.clone())
                    .or_default()
                    .push(atom.residue_seq);
            }
        }

        if chain_residues.is_empty() {
            return 0.0;
        }

        let mut total_expected = 0i32;
        let mut total_actual = 0usize;

        for residues in chain_residues.values() {
            if residues.is_empty() {
                continue;
            }

            let min_res = *residues.iter().min().unwrap();
            let max_res = *residues.iter().max().unwrap();
            let expected = max_res - min_res + 1;
            let actual = residues.len();

            total_expected += expected;
            total_actual += actual;
        }

        if total_expected <= 0 {
            return 0.0;
        }

        // Handle case where actual > expected (shouldn't happen but can with unusual numbering)
        if total_actual >= total_expected as usize {
            return 0.0;
        }

        let missing = total_expected as usize - total_actual;
        missing as f64 / total_expected as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_test_structure() -> PdbStructure {
        let mut structure = PdbStructure::new();

        // Create a simple structure with known composition
        // 2 ALA, 2 GLY, 1 VAL = 5 residues
        // Hydrophobic: ALA (2) + VAL (1) = 3/5 = 0.6
        // Glycine: 2/5 = 0.4
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1),
            create_ca_atom(2, "ALA", "A", 2),
            create_ca_atom(3, "GLY", "A", 3),
            create_ca_atom(4, "GLY", "A", 4),
            create_ca_atom(5, "VAL", "A", 5),
        ];

        structure
    }

    fn create_ca_atom(serial: i32, residue_name: &str, chain_id: &str, residue_seq: i32) -> Atom {
        Atom {
            serial,
            name: " CA ".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x: serial as f64,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[test]
    fn test_aa_composition() {
        let structure = create_test_structure();
        let composition = structure.aa_composition();

        assert_eq!(composition.len(), 3); // ALA, GLY, VAL

        let ala = composition.get("ALA").unwrap();
        let gly = composition.get("GLY").unwrap();
        let val = composition.get("VAL").unwrap();

        assert!((ala - 0.4).abs() < 1e-10); // 2/5
        assert!((gly - 0.4).abs() < 1e-10); // 2/5
        assert!((val - 0.2).abs() < 1e-10); // 1/5
    }

    #[test]
    fn test_glycine_ratio() {
        let structure = create_test_structure();
        let gly_ratio = structure.glycine_ratio();

        assert!((gly_ratio - 0.4).abs() < 1e-10); // 2/5
    }

    #[test]
    fn test_hydrophobic_ratio() {
        let structure = create_test_structure();
        let hydro_ratio = structure.hydrophobic_ratio();

        // ALA (2) + VAL (1) = 3/5 = 0.6
        assert!((hydro_ratio - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_count_ca_residues() {
        let structure = create_test_structure();
        assert_eq!(structure.count_ca_residues(), 5);
    }

    #[test]
    fn test_missing_residue_ratio_no_gaps() {
        let structure = create_test_structure();
        let missing = structure.missing_residue_ratio();

        // Residues 1-5 are all present, no gaps
        assert!((missing - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_missing_residue_ratio_with_gaps() {
        let mut structure = PdbStructure::new();

        // Residues 1, 3, 5 present (2 and 4 missing)
        // Expected: 5, Actual: 3, Missing: 2/5 = 0.4
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1),
            create_ca_atom(2, "ALA", "A", 3),
            create_ca_atom(3, "ALA", "A", 5),
        ];

        let missing = structure.missing_residue_ratio();
        assert!((missing - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_empty_structure_composition() {
        let structure = PdbStructure::new();

        let composition = structure.aa_composition();
        assert!(composition.is_empty());

        assert_eq!(structure.glycine_ratio(), 0.0);
        assert_eq!(structure.hydrophobic_ratio(), 0.0);
        assert_eq!(structure.count_ca_residues(), 0);
        assert_eq!(structure.missing_residue_ratio(), 0.0);
    }

    #[test]
    fn test_composition_ignores_non_ca_atoms() {
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1),
            Atom {
                serial: 2,
                name: " N  ".to_string(), // Not CA
                alt_loc: None,
                residue_name: "ALA".to_string(),
                chain_id: "A".to_string(),
                residue_seq: 1,
                ins_code: None,
                is_hetatm: false,
                x: 0.0,
                y: 0.0,
                z: 0.0,
                occupancy: 1.0,
                temp_factor: 20.0,
                element: "N".to_string(),
            },
            create_ca_atom(3, "GLY", "A", 2),
        ];

        let composition = structure.aa_composition();

        // Should only count 2 residues (ALA and GLY), not the N atom
        assert_eq!(structure.count_ca_residues(), 2);
        assert!((composition.get("ALA").unwrap() - 0.5).abs() < 1e-10);
        assert!((composition.get("GLY").unwrap() - 0.5).abs() < 1e-10);
    }
}
