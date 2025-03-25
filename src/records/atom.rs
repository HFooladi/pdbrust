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