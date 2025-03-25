//! CONECT record structure and implementations

/// Represents a CONECT record from a PDB file.
///
/// Specifies connectivity between atoms in the structure.
#[derive(Debug, Clone)]
pub struct Conect {
    /// Serial number of the central atom.
    pub atom_serial: i32,
    /// Serial numbers of atoms bonded to the central atom.
    pub bonded_atoms: Vec<i32>,
}