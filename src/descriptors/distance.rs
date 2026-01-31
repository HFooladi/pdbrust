//! Distance matrix and contact map calculations for protein structures.
//!
//! Functions for computing pairwise distances and contact maps
//! based on atom positions. These are essential for machine learning
//! applications (GNNs, protein transformers) and structural analysis.

use crate::core::PdbStructure;

/// Default contact threshold for Cα-Cα contacts in Angstroms.
/// 8.0 Å is the standard threshold used in most structural biology tools.
pub const DEFAULT_CA_CONTACT_THRESHOLD: f64 = 8.0;

/// Default contact threshold for all-atom contacts in Angstroms.
/// 4.5 Å is typical for atom-atom contacts.
pub const DEFAULT_ATOM_CONTACT_THRESHOLD: f64 = 4.5;

impl PdbStructure {
    /// Compute distance matrix for Cα atoms.
    ///
    /// Returns an N×N symmetric matrix where `matrix[i][j]` is the Euclidean
    /// distance in Angstroms between Cα atoms i and j.
    ///
    /// # Complexity
    ///
    /// O(n²) time and space where n is the number of Cα atoms.
    ///
    /// # Returns
    ///
    /// A 2D vector where `result[i][j]` is the distance between Cα atoms i and j.
    /// Returns an empty vector for structures with no Cα atoms.
    /// The diagonal elements are always 0.0.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let dist_matrix = structure.distance_matrix_ca();
    ///
    /// // Distance between first and second CA atoms
    /// println!("Distance: {:.2} Å", dist_matrix[0][1]);
    /// ```
    pub fn distance_matrix_ca(&self) -> Vec<Vec<f64>> {
        let ca_coords: Vec<(f64, f64, f64)> = self
            .atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        compute_distance_matrix(&ca_coords)
    }

    /// Compute contact map for Cα atoms.
    ///
    /// Returns an N×N boolean matrix where `matrix[i][j]` is true if the
    /// distance between Cα atoms i and j is less than or equal to the threshold.
    ///
    /// # Arguments
    ///
    /// * `threshold` - Distance threshold in Angstroms. Use 8.0 Å as a standard
    ///   value for Cα-Cα contacts.
    ///
    /// # Returns
    ///
    /// A 2D boolean vector where `result[i][j]` is true if atoms i and j are
    /// within the threshold distance. The diagonal elements are always true.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let contacts = structure.contact_map_ca(8.0);
    ///
    /// // Check if residues 0 and 5 are in contact
    /// if contacts[0][5] {
    ///     println!("Residues 0 and 5 are in contact");
    /// }
    /// ```
    pub fn contact_map_ca(&self, threshold: f64) -> Vec<Vec<bool>> {
        let ca_coords: Vec<(f64, f64, f64)> = self
            .atoms
            .iter()
            .filter(|atom| atom.name.trim() == "CA")
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        compute_contact_map(&ca_coords, threshold)
    }

    /// Compute distance matrix for all atoms.
    ///
    /// Returns an N×N symmetric matrix where `matrix[i][j]` is the Euclidean
    /// distance in Angstroms between atoms i and j.
    ///
    /// # Warning
    ///
    /// For large structures, this can be very memory-intensive. A structure
    /// with 10,000 atoms will require ~800 MB for the distance matrix.
    /// Consider using `distance_matrix_ca()` for most applications.
    ///
    /// # Complexity
    ///
    /// O(n²) time and space where n is the total number of atoms.
    ///
    /// # Returns
    ///
    /// A 2D vector where `result[i][j]` is the distance between atoms i and j.
    /// Returns an empty vector for empty structures.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let dist_matrix = structure.distance_matrix();
    ///
    /// println!("Matrix size: {} x {}", dist_matrix.len(), dist_matrix.len());
    /// ```
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        let coords: Vec<(f64, f64, f64)> = self
            .atoms
            .iter()
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        compute_distance_matrix(&coords)
    }

    /// Compute contact map for all atoms.
    ///
    /// Returns an N×N boolean matrix where `matrix[i][j]` is true if the
    /// distance between atoms i and j is less than or equal to the threshold.
    ///
    /// # Arguments
    ///
    /// * `threshold` - Distance threshold in Angstroms. Use 4.5 Å as a typical
    ///   value for atom-atom contacts.
    ///
    /// # Warning
    ///
    /// For large structures, this can be memory-intensive. Consider using
    /// `contact_map_ca()` for most applications.
    ///
    /// # Returns
    ///
    /// A 2D boolean vector where `result[i][j]` is true if atoms i and j are
    /// within the threshold distance.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use pdbrust::PdbStructure;
    ///
    /// let structure = PdbStructure::from_file("protein.pdb")?;
    /// let contacts = structure.contact_map(4.5);
    ///
    /// // Count total contacts
    /// let total: usize = contacts.iter().map(|row| row.iter().filter(|&&x| x).count()).sum();
    /// println!("Total atom-atom contacts: {}", total);
    /// ```
    pub fn contact_map(&self, threshold: f64) -> Vec<Vec<bool>> {
        let coords: Vec<(f64, f64, f64)> = self
            .atoms
            .iter()
            .map(|atom| (atom.x, atom.y, atom.z))
            .collect();

        compute_contact_map(&coords, threshold)
    }
}

/// Compute the distance matrix for a set of 3D coordinates.
///
/// This is an internal helper function used by both CA and all-atom variants.
fn compute_distance_matrix(coords: &[(f64, f64, f64)]) -> Vec<Vec<f64>> {
    if coords.is_empty() {
        return vec![];
    }

    let n = coords.len();
    let mut matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let (x1, y1, z1) = coords[i];
            let (x2, y2, z2) = coords[j];

            let dist = ((x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2)).sqrt();

            matrix[i][j] = dist;
            matrix[j][i] = dist;
        }
    }

    matrix
}

/// Compute the contact map for a set of 3D coordinates with a distance threshold.
///
/// This is an internal helper function used by both CA and all-atom variants.
fn compute_contact_map(coords: &[(f64, f64, f64)], threshold: f64) -> Vec<Vec<bool>> {
    if coords.is_empty() {
        return vec![];
    }

    let n = coords.len();
    let threshold_sq = threshold * threshold;
    let mut matrix = vec![vec![false; n]; n];

    // Diagonal is always true (distance to self is 0)
    for (i, row) in matrix.iter_mut().enumerate() {
        row[i] = true;
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let (x1, y1, z1) = coords[i];
            let (x2, y2, z2) = coords[j];

            let dist_sq = (x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2);

            let in_contact = dist_sq <= threshold_sq;
            matrix[i][j] = in_contact;
            matrix[j][i] = in_contact;
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::records::Atom;

    fn create_ca_atom(
        serial: i32,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
    ) -> Atom {
        Atom {
            serial,
            name: " CA ".to_string(),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: "C".to_string(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn create_atom(
        serial: i32,
        name: &str,
        residue_name: &str,
        chain_id: &str,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
    ) -> Atom {
        Atom {
            serial,
            name: format!("{:>4}", name),
            alt_loc: None,
            residue_name: residue_name.to_string(),
            chain_id: chain_id.to_string(),
            residue_seq,
            ins_code: None,
            is_hetatm: false,
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 20.0,
            element: name.chars().next().unwrap().to_string(),
        }
    }

    #[test]
    fn test_distance_matrix_ca_empty() {
        let structure = PdbStructure::new();
        let matrix = structure.distance_matrix_ca();
        assert!(matrix.is_empty());
    }

    #[test]
    fn test_distance_matrix_ca_single_atom() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0)];

        let matrix = structure.distance_matrix_ca();
        assert_eq!(matrix.len(), 1);
        assert_eq!(matrix[0].len(), 1);
        assert_eq!(matrix[0][0], 0.0);
    }

    #[test]
    fn test_distance_matrix_ca_two_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "ALA", "A", 2, 3.0, 4.0, 0.0), // distance = 5.0
        ];

        let matrix = structure.distance_matrix_ca();
        assert_eq!(matrix.len(), 2);

        // Diagonal should be 0
        assert_eq!(matrix[0][0], 0.0);
        assert_eq!(matrix[1][1], 0.0);

        // Off-diagonal should be 5.0
        assert!((matrix[0][1] - 5.0).abs() < 1e-10);
        assert!((matrix[1][0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_ca_symmetry() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "GLY", "A", 2, 1.0, 2.0, 3.0),
            create_ca_atom(3, "VAL", "A", 3, 4.0, 5.0, 6.0),
        ];

        let matrix = structure.distance_matrix_ca();

        // Verify symmetry
        for (i, row) in matrix.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - matrix[j][i]).abs() < 1e-10,
                    "Matrix should be symmetric at [{},{}]",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_distance_matrix_ca_known_distances() {
        let mut structure = PdbStructure::new();
        // Create atoms at known positions
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "GLY", "A", 2, 1.0, 0.0, 0.0), // distance = 1.0
            create_ca_atom(3, "VAL", "A", 3, 0.0, 1.0, 0.0), // distance from 1 = 1.0, from 2 = sqrt(2)
        ];

        let matrix = structure.distance_matrix_ca();

        assert!((matrix[0][1] - 1.0).abs() < 1e-10);
        assert!((matrix[0][2] - 1.0).abs() < 1e-10);
        assert!((matrix[1][2] - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn test_contact_map_ca_empty() {
        let structure = PdbStructure::new();
        let matrix = structure.contact_map_ca(8.0);
        assert!(matrix.is_empty());
    }

    #[test]
    fn test_contact_map_ca_single_atom() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0)];

        let matrix = structure.contact_map_ca(8.0);
        assert_eq!(matrix.len(), 1);
        assert!(matrix[0][0]); // Self-contact is always true
    }

    #[test]
    fn test_contact_map_ca_threshold() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "GLY", "A", 2, 5.0, 0.0, 0.0), // distance = 5.0
            create_ca_atom(3, "VAL", "A", 3, 10.0, 0.0, 0.0), // distance from 1 = 10.0
        ];

        // With threshold 6.0, only adjacent atoms should be in contact
        let contacts_6 = structure.contact_map_ca(6.0);
        assert!(contacts_6[0][1]); // 5.0 <= 6.0
        assert!(contacts_6[1][2]); // 5.0 <= 6.0
        assert!(!contacts_6[0][2]); // 10.0 > 6.0

        // With threshold 4.0, no off-diagonal contacts
        let contacts_4 = structure.contact_map_ca(4.0);
        assert!(!contacts_4[0][1]); // 5.0 > 4.0
        assert!(!contacts_4[1][2]); // 5.0 > 4.0
        assert!(!contacts_4[0][2]); // 10.0 > 4.0

        // Diagonal should always be true
        assert!(contacts_4[0][0]);
        assert!(contacts_4[1][1]);
        assert!(contacts_4[2][2]);
    }

    #[test]
    fn test_contact_map_ca_symmetry() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "GLY", "A", 2, 3.0, 4.0, 0.0),
            create_ca_atom(3, "VAL", "A", 3, 6.0, 8.0, 0.0),
        ];

        let contacts = structure.contact_map_ca(8.0);

        for (i, row) in contacts.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert_eq!(
                    val, contacts[j][i],
                    "Contact map should be symmetric at [{},{}]",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_distance_matrix_all_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_atom(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0),
            create_atom(3, "C", "ALA", "A", 1, 3.0, 0.0, 0.0),
        ];

        let matrix = structure.distance_matrix();
        assert_eq!(matrix.len(), 3);

        // Check known distances
        assert!((matrix[0][1] - 1.5).abs() < 1e-10);
        assert!((matrix[1][2] - 1.5).abs() < 1e-10);
        assert!((matrix[0][2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_contact_map_all_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_atom(2, "CA", "ALA", "A", 1, 2.0, 0.0, 0.0),
            create_atom(3, "C", "ALA", "A", 1, 5.0, 0.0, 0.0),
        ];

        let contacts = structure.contact_map(3.0);
        assert!(contacts[0][1]); // 2.0 <= 3.0
        assert!(!contacts[0][2]); // 5.0 > 3.0
        assert!(contacts[1][2]); // 3.0 <= 3.0 (exactly on threshold)
    }

    #[test]
    fn test_contact_map_boundary_threshold() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_ca_atom(1, "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "GLY", "A", 2, 8.0, 0.0, 0.0), // exactly at threshold
        ];

        // At exactly threshold distance, should be in contact (<=)
        let contacts = structure.contact_map_ca(8.0);
        assert!(contacts[0][1]);

        // Just below threshold, should NOT be in contact
        let contacts_below = structure.contact_map_ca(7.99);
        assert!(!contacts_below[0][1]);
    }

    #[test]
    fn test_distance_matrix_ignores_non_ca_atoms() {
        let mut structure = PdbStructure::new();
        structure.atoms = vec![
            create_atom(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
            create_ca_atom(2, "ALA", "A", 1, 1.0, 0.0, 0.0),
            create_atom(3, "C", "ALA", "A", 1, 2.0, 0.0, 0.0),
            create_ca_atom(4, "GLY", "A", 2, 4.0, 0.0, 0.0),
            create_atom(5, "O", "GLY", "A", 2, 5.0, 0.0, 0.0),
        ];

        let matrix = structure.distance_matrix_ca();

        // Only 2 CA atoms, so 2x2 matrix
        assert_eq!(matrix.len(), 2);
        assert!((matrix[0][1] - 3.0).abs() < 1e-10); // Distance between CA atoms
    }
}
