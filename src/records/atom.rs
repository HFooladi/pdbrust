//! ATOM record structure and implementations

/// Represents an atom record from a PDB file.
///
/// Contains all standard PDB ATOM record fields including position,
/// identification, and thermal factor information.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atom serial number.
    pub serial: i32,
    /// Atom name.
    pub name: String,
    /// Alternate location indicator (if any).
    pub alt_loc: Option<char>,
    /// Residue name.
    pub residue_name: String,
    /// Chain identifier.
    pub chain_id: String,
    /// Residue sequence number.
    pub residue_seq: i32,
    /// X coordinate in Angstroms.
    pub x: f64,
    /// Y coordinate in Angstroms.
    pub y: f64,
    /// Z coordinate in Angstroms.
    pub z: f64,
    /// Occupancy.
    pub occupancy: f64,
    /// Temperature factor.
    pub temp_factor: f64,
    /// Element symbol.
    pub element: String,
    /// Insertion code.
    pub ins_code: Option<char>,
}

impl Atom {
    /// Creates a new Atom with the given parameters.
    pub fn new(
        serial: i32,
        name: String,
        alt_loc: Option<char>,
        residue_name: String,
        chain_id: String,
        residue_seq: i32,
        x: f64,
        y: f64,
        z: f64,
        occupancy: f64,
        temp_factor: f64,
        element: String,
        ins_code: Option<char>,
    ) -> Self {
        Self {
            serial,
            name,
            alt_loc,
            residue_name,
            chain_id,
            residue_seq,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
            ins_code,
        }
    }

    /// Creates a new Atom from a PDB ATOM/HETATM record line.
    pub fn from_pdb_line(line: &str) -> Result<Self, crate::error::PdbError> {
        if line.len() < 80 {
            return Err(crate::error::PdbError::InvalidRecord(
                "Line too short for ATOM/HETATM record".to_string(),
            ));
        }

        let serial = line[6..11].trim().parse().unwrap_or(0);
        let name = line[12..16].trim().to_string();
        let alt_loc = if line[16..17].trim().is_empty() {
            None
        } else {
            Some(line[16..17].chars().next().unwrap())
        };
        let residue_name = line[17..20].trim().to_string();
        let chain_id = line[21..22].trim().to_string();
        let residue_seq = line[22..26].trim().parse().unwrap_or(0);
        let ins_code = if line[26..27].trim().is_empty() {
            None
        } else {
            Some(line[26..27].chars().next().unwrap())
        };
        let x = line[30..38].trim().parse().unwrap_or(0.0);
        let y = line[38..46].trim().parse().unwrap_or(0.0);
        let z = line[46..54].trim().parse().unwrap_or(0.0);
        let occupancy = line[54..60].trim().parse().unwrap_or(1.0);
        let temp_factor = line[60..66].trim().parse().unwrap_or(0.0);
        let element = line[76..78].trim().to_string();

        Ok(Self {
            serial,
            name,
            alt_loc,
            residue_name,
            chain_id,
            residue_seq,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element,
            ins_code,
        })
    }

    /// Returns the 3D coordinates of the atom as a tuple.
    pub fn coordinates(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }

    /// Returns the distance between this atom and another atom.
    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Returns the squared distance between this atom and another atom.
    pub fn distance_squared_to(&self, other: &Atom) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        dx * dx + dy * dy + dz * dz
    }

    /// Gives the angle between the centers of three atoms in degrees.
    /// The angle is calculated as the angle between the two lines that include
    /// atoms [1, 2] and [2, 3]
    pub fn angle_between(&self, atom2: &Atom, atom3: &Atom) -> f64 {
        // Get the coordinates as individual values
        let (x1, y1, z1) = self.coordinates();
        let (x2, y2, z2) = atom2.coordinates();
        let (x3, y3, z3) = atom3.coordinates();
        
        // Calculate vectors between atoms
        let v1x = x2 - x1;
        let v1y = y2 - y1;
        let v1z = z2 - z1;
        
        let v2x = x3 - x2;
        let v2y = y3 - y2;
        let v2z = z3 - z2;
        
        // Calculate dot product
        let dot = v1x * v2x + v1y * v2y + v1z * v2z;
        
        // Calculate magnitudes
        let mag1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
        let mag2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt();
        
        // Calculate angle and convert to degrees
        (dot / (mag1 * mag2)).acos().to_degrees()
    }

    /// Returns true if this atom is a backbone atom (N, CA, C, O).
    pub fn is_backbone(&self) -> bool {
        matches!(
            self.name.trim(),
            "N" | "CA" | "C" | "O" | "N1" | "C1" | "O1" | "CA1"
        )
    }

    /// Returns true if this atom is a hydrogen atom.
    pub fn is_hydrogen(&self) -> bool {
        self.name.starts_with('H')
    }

    /// Returns true if this atom is a heavy atom (not hydrogen).
    pub fn is_heavy_atom(&self) -> bool {
        !self.is_hydrogen()
    }

    /// Returns the residue identifier as a tuple of (chain_id, residue_seq, ins_code).
    pub fn residue_id(&self) -> (&str, i32, Option<char>) {
        (&self.chain_id, self.residue_seq, self.ins_code)
    }

    /// Returns a formatted string representation of the atom's position.
    pub fn position_string(&self) -> String {
        format!("{:.3} {:.3} {:.3}", self.x, self.y, self.z)
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
}

impl Eq for Atom {}

impl std::hash::Hash for Atom {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.serial.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_creation() {
        let atom = Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            1.0,
            2.0,
            3.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "CA");
        assert_eq!(atom.residue_name, "ALA");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.residue_seq, 1);
        assert_eq!(atom.coordinates(), (1.0, 2.0, 3.0));
    }

    #[test]
    fn test_atom_distance() {
        let atom1 = Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        let atom2 = Atom::new(
            2,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            2,
            1.0,
            1.0,
            1.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        assert!((atom1.distance_to(&atom2) - 1.732050808).abs() < 1e-6);
    }

    #[test]
    fn test_atom_angle() {
        let atom1 = Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        let atom2 = Atom::new(
            2,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            2,
            1.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        let atom3 = Atom::new(
            3,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            3,
            1.0,
            1.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        
        // The expected angle should be 90 degrees (right angle) because:
        // Vector from atom1 to atom2 is (1,0,0) and vector from atom2 to atom3 is (0,1,0)
        let angle = atom1.angle_between(&atom2, &atom3);
        
        // Use appropriate precision comparison for floating-point values
        assert!((angle - 90.0).abs() < 1e-6, "Expected angle to be 90°, got {}°", angle);
        
        // Test another configuration for thoroughness
        let atom4 = Atom::new(
            4,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            4,
            2.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        
        // Vector from atom1 to atom2 is (1,0,0) and vector from atom2 to atom4 is (1,0,0)
        // These are parallel so the angle should be 0
        let angle2 = atom1.angle_between(&atom2, &atom4);
        assert!(angle2.abs() < 1e-6, "Expected angle to be 0°, got {}°", angle2);
        
        // Additional test for 180 degree angle
        let atom5 = Atom::new(
            5,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            5,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );
        
        // Vector from atom2 to atom5 is (-1,0,0) and vector from atom2 to atom4 is (1,0,0)
        // These are in opposite directions so the angle should be 180
        let angle3 = atom2.angle_between(&atom5, &atom4);
        assert!((angle3 - 180.0).abs() < 1e-6, "Expected angle to be 180°, got {}°", angle3);
    }

    #[test]
    fn test_atom_backbone() {
        let backbone_atom = Atom::new(
            1,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        let sidechain_atom = Atom::new(
            2,
            "CB".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        assert!(backbone_atom.is_backbone());
        assert!(!sidechain_atom.is_backbone());
    }

    #[test]
    fn test_atom_hydrogen() {
        let hydrogen_atom = Atom::new(
            1,
            "H".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "H".to_string(),
            None,
        );

        let heavy_atom = Atom::new(
            2,
            "CA".to_string(),
            None,
            "ALA".to_string(),
            "A".to_string(),
            1,
            0.0,
            0.0,
            0.0,
            1.0,
            20.0,
            "C".to_string(),
            None,
        );

        assert!(hydrogen_atom.is_hydrogen());
        assert!(!heavy_atom.is_hydrogen());
        assert!(heavy_atom.is_heavy_atom());
        assert!(!hydrogen_atom.is_heavy_atom());
    }
}
