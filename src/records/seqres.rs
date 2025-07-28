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

impl SeqRes {
    /// Creates a new SEQRES record.
    pub fn new(serial: i32, chain_id: String, num_residues: i32, residues: Vec<String>) -> Self {
        Self {
            serial,
            chain_id,
            num_residues,
            residues,
        }
    }

    /// Creates a new SEQRES record from a PDB SEQRES line.
    pub fn from_pdb_line(line: &str) -> Result<Self, crate::error::PdbError> {
        if line.len() < 19 {
            return Err(crate::error::PdbError::InvalidRecord(
                "Line too short for SEQRES record".to_string(),
            ));
        }

        let serial = line[8..11].trim().parse().unwrap_or(0);
        let chain_id = line[11..12].trim().to_string();
        let num_residues = line[13..17].trim().parse().unwrap_or(0);

        // Parse residue names (13 residues per line, 4 characters each)
        let mut residues = Vec::new();
        for i in 0..13 {
            let start = 19 + (i * 4);
            let end = start + 3;
            if end <= line.len() {
                let residue = line[start..end].trim().to_string();
                if !residue.is_empty() {
                    residues.push(residue);
                }
            }
        }

        Ok(Self {
            serial,
            chain_id,
            num_residues,
            residues,
        })
    }

    /// Returns the sequence as a single string.
    pub fn get_sequence(&self) -> String {
        self.residues.join("")
    }

    /// Returns the length of the sequence.
    pub fn get_length(&self) -> usize {
        self.residues.len()
    }

    /// Returns true if the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Returns the residue at the given index (0-based).
    pub fn get_residue(&self, index: usize) -> Option<&str> {
        self.residues.get(index).map(|r| r.as_str())
    }

    /// Returns an iterator over the residues.
    pub fn iter(&self) -> std::slice::Iter<'_, String> {
        self.residues.iter()
    }

    /// Returns true if this SEQRES record is for the given chain.
    pub fn is_chain(&self, chain_id: &str) -> bool {
        self.chain_id == chain_id
    }

    /// Returns the sequence as a one-letter code string.
    /// Note: This is a simplified version and may not handle all amino acids.
    pub fn get_one_letter_code(&self) -> String {
        self.residues
            .iter()
            .map(|res| match res.as_str() {
                "ALA" => 'A',
                "ARG" => 'R',
                "ASN" => 'N',
                "ASP" => 'D',
                "CYS" => 'C',
                "GLN" => 'Q',
                "GLU" => 'E',
                "GLY" => 'G',
                "HIS" => 'H',
                "ILE" => 'I',
                "LEU" => 'L',
                "LYS" => 'K',
                "MET" => 'M',
                "PHE" => 'F',
                "PRO" => 'P',
                "SER" => 'S',
                "THR" => 'T',
                "TRP" => 'W',
                "TYR" => 'Y',
                "VAL" => 'V',
                _ => 'X',
            })
            .collect()
    }
}

impl PartialEq for SeqRes {
    fn eq(&self, other: &Self) -> bool {
        self.chain_id == other.chain_id && self.residues == other.residues
    }
}

impl Eq for SeqRes {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seqres_creation() {
        let residues = vec!["ALA".to_string(), "GLY".to_string(), "SER".to_string()];
        let seqres = SeqRes::new(1, "A".to_string(), 3, residues.clone());

        assert_eq!(seqres.serial, 1);
        assert_eq!(seqres.chain_id, "A");
        assert_eq!(seqres.num_residues, 3);
        assert_eq!(seqres.residues, residues);
    }

    #[test]
    fn test_seqres_from_pdb_line() {
        let line = "SEQRES   1 A   21  ALA GLY SER THR VAL LEU ILE PRO PHE MET TRP CYS             ";
        let seqres = SeqRes::from_pdb_line(line).unwrap();

        assert_eq!(seqres.serial, 1);
        assert_eq!(seqres.chain_id, "A");
        assert_eq!(seqres.num_residues, 21);
        assert_eq!(seqres.get_length(), 12);
        assert_eq!(seqres.get_residue(0), Some("ALA"));
        assert_eq!(seqres.get_residue(11), Some("CYS"));
    }

    #[test]
    fn test_seqres_sequence() {
        let residues = vec!["ALA".to_string(), "GLY".to_string(), "SER".to_string()];
        let seqres = SeqRes::new(1, "A".to_string(), 3, residues);

        assert_eq!(seqres.get_sequence(), "ALAGLYSER");
    }

    #[test]
    fn test_seqres_one_letter_code() {
        let residues = vec!["ALA".to_string(), "GLY".to_string(), "SER".to_string()];
        let seqres = SeqRes::new(1, "A".to_string(), 3, residues);

        assert_eq!(seqres.get_one_letter_code(), "AGS");
    }

    #[test]
    fn test_seqres_chain() {
        let seqres = SeqRes::new(1, "A".to_string(), 3, vec![]);

        assert!(seqres.is_chain("A"));
        assert!(!seqres.is_chain("B"));
    }

    #[test]
    fn test_seqres_empty() {
        let empty_seqres = SeqRes::new(1, "A".to_string(), 0, vec![]);
        let non_empty_seqres = SeqRes::new(1, "A".to_string(), 1, vec!["ALA".to_string()]);

        assert!(empty_seqres.is_empty());
        assert!(!non_empty_seqres.is_empty());
    }
}
