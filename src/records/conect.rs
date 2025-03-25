//! CONECT record structure and implementations

/// Represents a CONECT record from a PDB file.
///
/// Specifies connectivity between atoms in the structure.
#[derive(Debug, Clone)]
pub struct Conect {
    /// Serial number of atom 1.
    pub atom1: i32,
    /// Serial number of atom 2.
    pub atom2: i32,
    /// Serial number of atom 3 (optional).
    pub atom3: Option<i32>,
    /// Serial number of atom 4 (optional).
    pub atom4: Option<i32>,
}
