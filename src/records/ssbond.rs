//! SSBOND record structure and implementations

/// Represents an SSBOND record from a PDB file.
///
/// Describes disulfide bonds between pairs of cysteine residues.
#[derive(Debug, Clone)]
pub struct SSBond {
    /// Serial number of the SSBOND record.
    pub serial: i32,
    /// Name of the first residue.
    pub residue1_name: String,
    /// Chain identifier for the first residue.
    pub chain1_id: String,
    /// Sequence number of the first residue.
    pub residue1_seq: i32,
    /// Insertion code for the first residue (optional).
    pub icode1: Option<char>,
    /// Name of the second residue.
    pub residue2_name: String,
    /// Chain identifier for the second residue.
    pub chain2_id: String,
    /// Sequence number of the second residue.
    pub residue2_seq: i32,
    /// Insertion code for the second residue (optional).
    pub icode2: Option<char>,
    /// Symmetry operator for the first residue.
    pub sym1: i32,
    /// Symmetry operator for the second residue.
    pub sym2: i32,
    /// Distance between the two sulfur atoms.
    pub length: f64,
}