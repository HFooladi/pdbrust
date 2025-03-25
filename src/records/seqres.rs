//! SEQRES record structure and implementations


/// Represents a SEQRES record from a PDB file.
///
/// Contains sequence information for a specific chain in the structure.
#[derive(Debug, Clone)]
pub struct SeqRes {
    /// Serial number of the SEQRES record.
    pub serial: i32,
    /// Chain identifier.
    pub chain_id: String,
    /// Number of residues in the sequence.
    pub num_residues: i32,
    /// List of residue names in the sequence.
    pub residues: Vec<String>,
}