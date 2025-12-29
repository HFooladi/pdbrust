//! Geometric descriptors for protein structures.
//!
//! Functions for computing geometric properties of protein structures
//! based on Cα atom positions.

use crate::core::PdbStructure;
use super::StructureDescriptors;

/// Threshold distance for consecutive Cα-Cα atoms in ordered structures.
/// Typical Cα-Cα distance in regular secondary structures is ~3.8 Å.
const ORDERED_CA_DISTANCE_THRESHOLD: f64 = 4.0;

impl PdbStructure {
    /// Calculate the radius of gyration.
    ///
    /// The radius of gyration (Rg) is a measure of the overall size and
    /// compactness of a protein structure. It represents the root-mean-square
    /// distance of all Cα atoms from the structure's centroid.
    ///
    /// # Formula
    ///
    /// Rg = sqrt(Σ(r_i - r_centroid)² / N)
    ///
    /// # Returns
    ///
    /// The radius of gyration in Angstroms. Returns 0.0 for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let rg = structure.radius_of_gyration();
    ///
    /// println!("Radius of gyration: {:.2} Å", rg);
    /// ```
    pub fn radius_of_gyration(&self) -> f64 {
        let ca_coords: Vec<(f64, f64, f64)> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        if ca_coords.is_empty() {
            return 0.0;
        }

        let n = ca_coords.len() as f64;

        // Calculate centroid
        let cx = ca_coords.iter().map(|c| c.0).sum::<f64>() / n;
        let cy = ca_coords.iter().map(|c| c.1).sum::<f64>() / n;
        let cz = ca_coords.iter().map(|c| c.2).sum::<f64>() / n;

        // Calculate sum of squared distances from centroid
        let rg_squared: f64 = ca_coords
            .iter()
            .map(|(x, y, z)| {
                (x - cx).powi(2) + (y - cy).powi(2) + (z - cz).powi(2)
            })
            .sum::<f64>() / n;

        rg_squared.sqrt()
    }

    /// Calculate the maximum Cα-Cα distance.
    ///
    /// This is the largest pairwise distance between any two Cα atoms
    /// in the structure, representing the maximum spatial extent.
    ///
    /// # Complexity
    ///
    /// O(n²) where n is the number of Cα atoms.
    ///
    /// # Returns
    ///
    /// The maximum Cα-Cα distance in Angstroms. Returns 0.0 for structures
    /// with fewer than 2 Cα atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let max_dist = structure.max_ca_distance();
    ///
    /// println!("Maximum extent: {:.2} Å", max_dist);
    /// ```
    pub fn max_ca_distance(&self) -> f64 {
        let ca_coords: Vec<(f64, f64, f64)> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        if ca_coords.len() < 2 {
            return 0.0;
        }

        let mut max_dist_sq = 0.0f64;

        for i in 0..ca_coords.len() {
            for j in (i + 1)..ca_coords.len() {
                let (x1, y1, z1) = ca_coords[i];
                let (x2, y2, z2) = ca_coords[j];

                let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2);

                if dist_sq > max_dist_sq {
                    max_dist_sq = dist_sq;
                }
            }
        }

        max_dist_sq.sqrt()
    }

    /// Estimate secondary structure content based on Cα-Cα distances.
    ///
    /// This is a heuristic measure based on the observation that
    /// consecutive Cα atoms in regular secondary structures (α-helices
    /// and β-sheets) are typically ~3.8 Å apart, while irregular
    /// regions may have larger distances.
    ///
    /// # Algorithm
    ///
    /// Counts the fraction of consecutive Cα-Cα pairs with distances
    /// below 4.0 Å (indicating ordered structure).
    ///
    /// # Returns
    ///
    /// The fraction of consecutive Cα pairs in ordered structure (0.0 to 1.0).
    /// Returns 0.0 for structures with fewer than 2 Cα atoms.
    ///
    /// # Note
    ///
    /// This is a rough heuristic, not a true DSSP-style secondary structure
    /// assignment. It's useful for rapid screening but not for detailed
    /// structural analysis.
    pub fn secondary_structure_ratio(&self) -> f64 {
        // Get CA atoms sorted by chain and residue number
        let mut ca_atoms: Vec<_> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .collect();

        if ca_atoms.len() < 2 {
            return 0.0;
        }

        // Sort by chain, then by residue number
        ca_atoms.sort_by(|a, b| {
            match a.chain_id.cmp(&b.chain_id) {
                std::cmp::Ordering::Equal => a.residue_seq.cmp(&b.residue_seq),
                other => other,
            }
        });

        let mut ordered_count = 0usize;
        let mut total_pairs = 0usize;

        for window in ca_atoms.windows(2) {
            let atom1 = window[0];
            let atom2 = window[1];

            // Only consider consecutive residues in the same chain
            if atom1.chain_id == atom2.chain_id
                && (atom2.residue_seq - atom1.residue_seq).abs() == 1
            {
                let dist_sq = (atom2.x - atom1.x).powi(2)
                    + (atom2.y - atom1.y).powi(2)
                    + (atom2.z - atom1.z).powi(2);
                let dist = dist_sq.sqrt();

                total_pairs += 1;
                if dist < ORDERED_CA_DISTANCE_THRESHOLD {
                    ordered_count += 1;
                }
            }
        }

        if total_pairs == 0 {
            return 0.0;
        }

        ordered_count as f64 / total_pairs as f64
    }

    /// Calculate the compactness index.
    ///
    /// The compactness index normalizes the radius of gyration by the
    /// number of residues to allow comparison between proteins of
    /// different sizes.
    ///
    /// # Formula
    ///
    /// Compactness = Rg / n^(1/3)
    ///
    /// Lower values indicate more compact (globular) structures,
    /// higher values indicate more extended structures.
    ///
    /// # Returns
    ///
    /// The compactness index. Returns 0.0 for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let compactness = structure.compactness_index();
    ///
    /// if compactness < 2.0 {
    ///     println!("Highly compact globular structure");
    /// } else if compactness > 3.0 {
    ///     println!("Extended or elongated structure");
    /// }
    /// ```
    pub fn compactness_index(&self) -> f64 {
        let n = self.count_ca_residues();
        if n < 2 {
            return 0.0;
        }

        let rg = self.radius_of_gyration();
        rg / (n as f64).powf(1.0 / 3.0)
    }

    /// Calculate the Cα atom density.
    ///
    /// The density is computed as the number of Cα atoms divided by
    /// the volume of the bounding box containing all Cα atoms.
    ///
    /// # Returns
    ///
    /// The Cα density in atoms per cubic Angstrom. Returns 0.0 if
    /// the bounding box has zero volume or no Cα atoms.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let density = structure.ca_density();
    ///
    /// println!("Cα density: {:.4} atoms/Å³", density);
    /// ```
    pub fn ca_density(&self) -> f64 {
        let ca_coords: Vec<(f64, f64, f64)> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        if ca_coords.is_empty() {
            return 0.0;
        }

        // Calculate bounding box
        let xs: Vec<f64> = ca_coords.iter().map(|c| c.0).collect();
        let ys: Vec<f64> = ca_coords.iter().map(|c| c.1).collect();
        let zs: Vec<f64> = ca_coords.iter().map(|c| c.2).collect();

        let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let z_min = zs.iter().cloned().fold(f64::INFINITY, f64::min);
        let z_max = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        let dx = x_max - x_min;
        let dy = y_max - y_min;
        let dz = z_max - z_min;

        // Handle degenerate cases (linear or planar arrangements)
        let volume = dx.max(0.001) * dy.max(0.001) * dz.max(0.001);

        ca_coords.len() as f64 / volume
    }

    /// Get the bounding box dimensions for Cα atoms.
    ///
    /// # Returns
    ///
    /// A tuple of ((x_min, x_max), (y_min, y_max), (z_min, z_max)).
    /// Returns ((0,0), (0,0), (0,0)) for empty structures.
    pub fn ca_bounding_box(&self) -> ((f64, f64), (f64, f64), (f64, f64)) {
        let ca_coords: Vec<(f64, f64, f64)> = self.atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        if ca_coords.is_empty() {
            return ((0.0, 0.0), (0.0, 0.0), (0.0, 0.0));
        }

        let xs: Vec<f64> = ca_coords.iter().map(|c| c.0).collect();
        let ys: Vec<f64> = ca_coords.iter().map(|c| c.1).collect();
        let zs: Vec<f64> = ca_coords.iter().map(|c| c.2).collect();

        let x_min = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let z_min = zs.iter().cloned().fold(f64::INFINITY, f64::min);
        let z_max = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        ((x_min, x_max), (y_min, y_max), (z_min, z_max))
    }

    /// Compute all structure descriptors at once.
    ///
    /// This is more efficient than calling individual descriptor methods
    /// when you need multiple descriptors, as some intermediate calculations
    /// can be shared.
    ///
    /// # Returns
    ///
    /// A `StructureDescriptors` struct containing all computed descriptors.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let descriptors = structure.structure_descriptors();
    ///
    /// println!("Residues: {}", descriptors.num_residues);
    /// println!("Rg: {:.2} Å", descriptors.radius_of_gyration);
    /// println!("Compactness: {:.2}", descriptors.compactness_index);
    /// ```
    pub fn structure_descriptors(&self) -> StructureDescriptors {
        StructureDescriptors {
            num_residues: self.count_ca_residues(),
            num_atoms: self.get_num_atoms(),
            aa_composition: self.aa_composition(),
            glycine_ratio: self.glycine_ratio(),
            hydrophobic_ratio: self.hydrophobic_ratio(),
            radius_of_gyration: self.radius_of_gyration(),
            max_ca_distance: self.max_ca_distance(),
            missing_residue_ratio: self.missing_residue_ratio(),
            secondary_structure_ratio: self.secondary_structure_ratio(),
            compactness_index: self.compactness_index(),
            ca_density: self.ca_density(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_ca_atom(serial: i32, residue_name: &str, chain_id: &str, residue_seq: i32, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial,
            name: " CA ".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    fn create_linear_structure() -> PdbStructure {
        // Create a linear structure along the x-axis
        // Points at (0,0,0), (3.8,0,0), (7.6,0,0), (11.4,0,0), (15.2,0,0)
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "ALA", "A", 2, 3.8, 0.0, 0.0),
            create_ca_atom(3, "ALA", "A", 3, 7.6, 0.0, 0.0),
            create_ca_atom(4, "ALA", "A", 4, 11.4, 0.0, 0.0),
            create_ca_atom(5, "ALA", "A", 5, 15.2, 0.0, 0.0),
        ];

        structure
    }

    fn create_cubic_structure() -> PdbStructure {
        // Create 8 atoms at corners of a cube with side length 10
        let mut structure = PdbStructure::new();

        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "ALA", "A", 2, 10.0, 0.0, 0.0),
            create_ca_atom(3, "ALA", "A", 3, 0.0, 10.0, 0.0),
            create_ca_atom(4, "ALA", "A", 4, 10.0, 10.0, 0.0),
            create_ca_atom(5, "ALA", "A", 5, 0.0, 0.0, 10.0),
            create_ca_atom(6, "ALA", "A", 6, 10.0, 0.0, 10.0),
            create_ca_atom(7, "ALA", "A", 7, 0.0, 10.0, 10.0),
            create_ca_atom(8, "ALA", "A", 8, 10.0, 10.0, 10.0),
        ];

        structure
    }

    #[test]
    fn test_radius_of_gyration_linear() {
        let structure = create_linear_structure();
        let rg = structure.radius_of_gyration();

        // For 5 points at 0, 3.8, 7.6, 11.4, 15.2
        // Centroid at 7.6
        // Distances: 7.6, 3.8, 0, 3.8, 7.6
        // Sum of squares: 7.6² + 3.8² + 0 + 3.8² + 7.6² = 57.76 + 14.44 + 14.44 + 57.76 = 144.4
        // Mean = 144.4 / 5 = 28.88
        // Rg = sqrt(28.88) ≈ 5.374
        assert!((rg - 5.374).abs() < 0.01, "Rg = {}", rg);
    }

    #[test]
    fn test_radius_of_gyration_empty() {
        let structure = PdbStructure::new();
        assert_eq!(structure.radius_of_gyration(), 0.0);
    }

    #[test]
    fn test_max_ca_distance_linear() {
        let structure = create_linear_structure();
        let max_dist = structure.max_ca_distance();

        // Maximum distance is from first to last point: 15.2 - 0 = 15.2
        assert!((max_dist - 15.2).abs() < 1e-10);
    }

    #[test]
    fn test_max_ca_distance_cubic() {
        let structure = create_cubic_structure();
        let max_dist = structure.max_ca_distance();

        // Maximum distance is the space diagonal of the cube: sqrt(10² + 10² + 10²) = sqrt(300) ≈ 17.32
        let expected = (300.0f64).sqrt();
        assert!((max_dist - expected).abs() < 0.01, "max_dist = {}", max_dist);
    }

    #[test]
    fn test_max_ca_distance_single_atom() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
        ];

        assert_eq!(structure.max_ca_distance(), 0.0);
    }

    #[test]
    fn test_secondary_structure_ratio() {
        let structure = create_linear_structure();
        let ss_ratio = structure.secondary_structure_ratio();

        // All consecutive pairs are 3.8 Å apart, which is < 4.0 threshold
        // So ratio should be 1.0 (all ordered)
        assert!((ss_ratio - 1.0).abs() < 1e-10, "ss_ratio = {}", ss_ratio);
    }

    #[test]
    fn test_secondary_structure_ratio_disordered() {
        let mut structure = PdbStructure::new();

        // Create structure with large gaps between consecutive residues
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "ALA", "A", 2, 10.0, 0.0, 0.0), // 10 Å > 4 Å threshold
            create_ca_atom(3, "ALA", "A", 3, 20.0, 0.0, 0.0),
            create_ca_atom(4, "ALA", "A", 4, 30.0, 0.0, 0.0),
        ];

        let ss_ratio = structure.secondary_structure_ratio();

        // All pairs are 10 Å apart, which is > 4.0 threshold
        // So ratio should be 0.0 (all disordered)
        assert!((ss_ratio - 0.0).abs() < 1e-10, "ss_ratio = {}", ss_ratio);
    }

    #[test]
    fn test_compactness_index() {
        let structure = create_linear_structure();
        let compactness = structure.compactness_index();

        // Rg ≈ 5.374, n = 5
        // Compactness = 5.374 / 5^(1/3) = 5.374 / 1.71 ≈ 3.14
        let n = 5.0f64;
        let expected = 5.374 / n.powf(1.0 / 3.0);
        assert!((compactness - expected).abs() < 0.1, "compactness = {}", compactness);
    }

    #[test]
    fn test_ca_density_cubic() {
        let structure = create_cubic_structure();
        let density = structure.ca_density();

        // 8 atoms in 10x10x10 = 1000 Å³ volume
        // Density = 8 / 1000 = 0.008
        assert!((density - 0.008).abs() < 0.001, "density = {}", density);
    }

    #[test]
    fn test_ca_bounding_box() {
        let structure = create_cubic_structure();
        let ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = structure.ca_bounding_box();

        assert!((x_min - 0.0).abs() < 1e-10);
        assert!((x_max - 10.0).abs() < 1e-10);
        assert!((y_min - 0.0).abs() < 1e-10);
        assert!((y_max - 10.0).abs() < 1e-10);
        assert!((z_min - 0.0).abs() < 1e-10);
        assert!((z_max - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_structure_descriptors() {
        let structure = create_linear_structure();
        let desc = structure.structure_descriptors();

        assert_eq!(desc.num_residues, 5);
        assert_eq!(desc.num_atoms, 5);
        assert!(desc.radius_of_gyration > 0.0);
        assert!(desc.max_ca_distance > 0.0);
        assert!(desc.compactness_index > 0.0);
    }

    #[test]
    fn test_empty_structure_geometry() {
        let structure = PdbStructure::new();

        assert_eq!(structure.radius_of_gyration(), 0.0);
        assert_eq!(structure.max_ca_distance(), 0.0);
        assert_eq!(structure.secondary_structure_ratio(), 0.0);
        assert_eq!(structure.compactness_index(), 0.0);
        assert_eq!(structure.ca_density(), 0.0);
    }
}
