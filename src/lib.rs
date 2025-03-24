use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PdbError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Invalid record format: {0}")]
    InvalidRecord(String),
    #[error("Parse error: {0}")]
    ParseError(String),
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: Option<char>,
    pub residue_name: String,
    pub chain_id: String,
    pub residue_seq: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: String,
}

#[derive(Debug, Clone)]
pub struct SeqRes {
    pub serial: i32,
    pub chain_id: String,
    pub num_residues: i32,
    pub residues: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct Conect {
    pub atom_serial: i32,
    pub bonded_atoms: Vec<i32>,
}

#[derive(Debug, Clone)]
pub struct SSBond {
    pub serial: i32,
    pub residue1_name: String,
    pub chain1_id: String,
    pub residue1_seq: i32,
    pub residue2_name: String,
    pub chain2_id: String,
    pub residue2_seq: i32,
    pub distance: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct Remark {
    pub number: i32,
    pub content: String,
}

#[derive(Debug, Clone)]
pub struct Model {
    pub serial: i32,
    pub atoms: Vec<Atom>,
    pub remarks: Vec<Remark>,
}

#[derive(Debug, Default)]
pub struct PdbStructure {
    pub atoms: Vec<Atom>,
    pub header: Option<String>,
    pub title: Option<String>,
    pub seqres: Vec<SeqRes>,
    pub connects: Vec<Conect>,
    pub ssbonds: Vec<SSBond>,
    pub remarks: Vec<Remark>,
    pub models: Vec<Model>,
    pub current_model: Option<i32>,
}


impl PdbStructure {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            header: None,
            title: None,
            seqres: Vec::new(),
            connects: Vec::new(),
            ssbonds: Vec::new(),
            remarks: Vec::new(),
            models: Vec::new(),
            current_model: None,
        }
    }

    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, PdbError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut structure = PdbStructure::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                let atom = Self::parse_atom_record(&line)?;
                
                // If we're in a model, add it to the current model
                // Otherwise, add it to the general atoms list
                if let Some(model_num) = structure.current_model {
                    // Find the model and add the atom to it
                    if let Some(model) = structure.models.iter_mut()
                        .find(|m| m.serial == model_num) {
                        model.atoms.push(atom.clone());
                    }
                }
                
                // Always add to the general atoms list for backward compatibility
                structure.atoms.push(atom);
            } else if line.starts_with("HEADER") {
                structure.header = Some(line[6..].trim().to_string());
            } else if line.starts_with("TITLE ") {
                structure.title = Some(line[6..].trim().to_string());
            } else if line.starts_with("MODEL ") {
                let model_num = Self::parse_model_record(&line)?;
                structure.current_model = Some(model_num);
                structure.models.push(Model {
                    serial: model_num,
                    atoms: Vec::new(),
                    remarks: Vec::new(),
                });
            } else if line.starts_with("ENDMDL") {
                structure.current_model = None;
            } else if line.starts_with("SEQRES") {
                let seqres = Self::parse_seqres_record(&line)?;
                structure.seqres.push(seqres);
            } else if line.starts_with("CONECT") {
                let conect = Self::parse_conect_record(&line)?;
                structure.connects.push(conect);
            } else if line.starts_with("SSBOND") {
                let ssbond = Self::parse_ssbond_record(&line)?;
                structure.ssbonds.push(ssbond);
            } else if line.starts_with("REMARK") {
                let remark = Self::parse_remark_record(&line)?;
                structure.remarks.push(remark);
            }
        }

        Ok(structure)
    }

    fn parse_atom_record(line: &str) -> Result<Atom, PdbError> {
        if line.len() < 66 {
            return Err(PdbError::InvalidRecord(
                "ATOM record line too short".to_string(),
            ));
        }

        let serial = line[6..11]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid atom serial number".to_string()))?;

        let name = line[12..16].trim().to_string();
        let alt_loc = {
            let c = line.chars().nth(16).unwrap_or(' ');
            if c == ' ' {
                None
            } else {
                Some(c)
            }
        };

        let residue_name = line[17..20].trim().to_string();
        let chain_id = line[21..22].trim().to_string();
        
        let residue_seq = line[22..26]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid residue sequence number".to_string()))?;

        let x = line[30..38]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid x coordinate".to_string()))?;
        let y = line[38..46]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid y coordinate".to_string()))?;
        let z = line[46..54]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid z coordinate".to_string()))?;

        let occupancy = line[54..60]
            .trim()
            .parse()
            .unwrap_or(1.0);
        let temp_factor = line[60..66]
            .trim()
            .parse()
            .unwrap_or(0.0);

        let element = if line.len() >= 78 {
            line[76..78].trim().to_string()
        } else {
            "".to_string()
        };

        Ok(Atom {
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
        })
    }

    fn parse_model_record(line: &str) -> Result<i32, PdbError> {
        if line.len() < 14 {
            return Err(PdbError::InvalidRecord(
                "MODEL record line too short".to_string(),
            ));
        }

        line[10..14]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid model serial number".to_string()))
    }

    fn parse_seqres_record(line: &str) -> Result<SeqRes, PdbError> {
        if line.len() < 70 {
            return Err(PdbError::InvalidRecord(
                "SEQRES record line too short".to_string(),
            ));
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(PdbError::InvalidRecord(
                "SEQRES record must have at least serial, chain ID, and number of residues".to_string(),
            ));
        }

        let serial = parts[1]
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid SEQRES serial number".to_string()))?;

        let chain_id = parts[2].to_string();
        
        let num_residues = parts[3]
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid number of residues".to_string()))?;

        let residues = parts[4..].iter().map(|s| s.to_string()).collect();

        Ok(SeqRes {
            serial,
            chain_id,
            num_residues,
            residues,
        })
    }

    fn parse_conect_record(line: &str) -> Result<Conect, PdbError> {
        if line.len() < 11 {
            return Err(PdbError::InvalidRecord(
                "CONECT record line too short".to_string(),
            ));
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            return Err(PdbError::InvalidRecord(
                "CONECT record must have at least one atom serial number".to_string(),
            ));
        }

        let atom_serial = parts[1]
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid atom serial number".to_string()))?;

        let bonded_atoms: Vec<i32> = parts[2..]
            .iter()
            .filter_map(|s| s.parse().ok())
            .collect();

        Ok(Conect {
            atom_serial,
            bonded_atoms,
        })
    }

    fn parse_ssbond_record(line: &str) -> Result<SSBond, PdbError> {
        if line.len() < 38 {
            return Err(PdbError::InvalidRecord(
                "SSBOND record line too short".to_string(),
            ));
        }

        let serial = line[7..10]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid SSBOND serial number".to_string()))?;

        let residue1_name = line[11..14].trim().to_string();
        let chain1_id = line[15..16].trim().to_string();
        
        let residue1_seq = line[17..21]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid first residue sequence number".to_string()))?;

        let residue2_name = line[25..28].trim().to_string();
        let chain2_id = line[29..30].trim().to_string();
        
        let residue2_seq = line[31..35]
            .trim()
            .parse()
            .map_err(|_| PdbError::ParseError("Invalid second residue sequence number".to_string()))?;
            
        // Try to find the distance value by looking for a number at the end of the line
        let distance = line[35..]
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .next();

        Ok(SSBond {
            serial,
            residue1_name,
            chain1_id,
            residue1_seq,
            residue2_name,
            chain2_id,
            residue2_seq,
            distance,
        })
    }

    fn parse_remark_record(line: &str) -> Result<Remark, PdbError> {
        if line.len() < 11 {
            return Err(PdbError::InvalidRecord(
                "REMARK record line too short".to_string(),
            ));
        }

        // REMARK number is in columns 7-10
        let number = line[7..10]
            .trim()
            .parse()
            .unwrap_or(0); // Default to 0 if not a valid number

        // Content starts at column 11
        let content = line[11..].trim().to_string();

        Ok(Remark { number, content })
    }

    /// Get all unique chain IDs in the structure
    pub fn get_chain_ids(&self) -> Vec<String> {
        let mut chain_ids: Vec<String> = self.atoms
            .iter()
            .map(|atom| atom.chain_id.clone())
            .collect();
        chain_ids.sort();
        chain_ids.dedup();
        chain_ids
    }

    /// Get all residues for a specific chain
    pub fn get_residues_for_chain(&self, chain_id: &str) -> Vec<(i32, String)> {
        let mut residues: Vec<(i32, String)> = self.atoms
            .iter()
            .filter(|atom| atom.chain_id == chain_id)
            .map(|atom| (atom.residue_seq, atom.residue_name.clone()))
            .collect();
        residues.sort();
        residues.dedup();
        residues
    }

    /// Get sequence for a specific chain from SEQRES records
    pub fn get_sequence(&self, chain_id: &str) -> Vec<String> {
        self.seqres
            .iter()
            .filter(|sr| sr.chain_id == chain_id)
            .flat_map(|sr| sr.residues.clone())
            .collect()
    }

    /// Get all remarks with a specific number
    pub fn get_remarks_by_number(&self, number: i32) -> Vec<&Remark> {
        self.remarks
            .iter()
            .filter(|r| r.number == number)
            .collect()
    }

    /// Get all atoms connected to a specific atom
    pub fn get_connected_atoms(&self, atom_serial: i32) -> Vec<&Atom> {
        // First find all atoms that are connected to this atom
        let connected_serials: Vec<i32> = self.connects
            .iter()
            .filter(|c| c.atom_serial == atom_serial)
            .flat_map(|c| c.bonded_atoms.clone())
            .collect();

        // Then get the actual atom references
        self.atoms
            .iter()
            .filter(|a| connected_serials.contains(&a.serial))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_atom_record() {
        let atom_line = "ATOM      1  N   ASP A  30      31.904  -0.904  -0.904  1.00 20.00           N  ";
        let atom = PdbStructure::parse_atom_record(atom_line).unwrap();
        
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.residue_name, "ASP");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.residue_seq, 30);
        assert!((atom.x - 31.904).abs() < 1e-6);
        assert!((atom.y - (-0.904)).abs() < 1e-6);
        assert!((atom.z - (-0.904)).abs() < 1e-6);
        assert!((atom.occupancy - 1.0).abs() < 1e-6);
        assert!((atom.temp_factor - 20.0).abs() < 1e-6);
        assert_eq!(atom.element, "N");
    }

    #[test]
    fn test_parse_seqres_record() {
        let seqres_line = "SEQRES   1 A  20  MET ALA CYS PRO GLY VAL GLY PRO GLU CYS ARG GLU MET SER PRO";
        let seqres = PdbStructure::parse_seqres_record(seqres_line).unwrap();
        
        assert_eq!(seqres.serial, 1);
        assert_eq!(seqres.chain_id, "A");
        assert_eq!(seqres.num_residues, 20);
        let expected_residues = vec!["MET", "ALA", "CYS", "PRO", "GLY", "VAL", "GLY", "PRO", "GLU", "CYS", "ARG", "GLU", "MET", "SER", "PRO"]
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        assert_eq!(seqres.residues, expected_residues);
    }

    #[test]
    fn test_parse_conect_record() {
        let conect_line = "CONECT 1179 1174 1177 1180";
        let conect = PdbStructure::parse_conect_record(conect_line).unwrap();
        
        assert_eq!(conect.atom_serial, 1179);
        let expected_bonds = vec![1174, 1177, 1180];
        assert_eq!(conect.bonded_atoms, expected_bonds);
    }

    #[test]
    fn test_parse_ssbond_record() {
        let ssbond_line = "SSBOND   1 CYS A   85    CYS A  101                               2.03";
        let ssbond = PdbStructure::parse_ssbond_record(ssbond_line).unwrap();
        
        assert_eq!(ssbond.serial, 1);
        assert_eq!(ssbond.residue1_name, "CYS");
        assert_eq!(ssbond.chain1_id, "A");
        assert_eq!(ssbond.residue1_seq, 85);
        assert_eq!(ssbond.residue2_name, "CYS");
        assert_eq!(ssbond.chain2_id, "A");
        assert_eq!(ssbond.residue2_seq, 101);
        assert!(ssbond.distance.is_some());
        assert!((ssbond.distance.unwrap() - 2.03).abs() < 1e-6);
    }

    #[test]
    fn test_new_structure() {
        let structure = PdbStructure::new();
        assert!(structure.atoms.is_empty());
        assert!(structure.header.is_none());
        assert!(structure.title.is_none());
        assert!(structure.seqres.is_empty());
        assert!(structure.connects.is_empty());
        assert!(structure.ssbonds.is_empty());
        assert!(structure.remarks.is_empty());
        assert!(structure.models.is_empty());
        assert!(structure.current_model.is_none());
    }

    #[test]
    fn test_get_chain_ids() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "B".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
        });

        let chain_ids = structure.get_chain_ids();
        assert_eq!(chain_ids, vec!["A", "B"]);
    }

    #[test]
    fn test_get_residues_for_chain() {
        let mut structure = PdbStructure::new();
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
        });

        let residues = structure.get_residues_for_chain("A");
        assert_eq!(residues, vec![(1, "ASP".to_string())]);
    }

    #[test]
    fn test_get_sequence() {
        let mut structure = PdbStructure::new();
        structure.seqres.push(SeqRes {
            serial: 1,
            chain_id: "A".to_string(),
            num_residues: 3,
            residues: vec!["ALA".to_string(), "GLY".to_string(), "VAL".to_string()],
        });

        let sequence = structure.get_sequence("A");
        assert_eq!(sequence, vec!["ALA", "GLY", "VAL"]);
    }

    #[test]
    fn test_get_connected_atoms() {
        let mut structure = PdbStructure::new();
        
        // Add some atoms
        structure.atoms.push(Atom {
            serial: 1,
            name: "N".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "N".to_string(),
        });
        structure.atoms.push(Atom {
            serial: 2,
            name: "CA".to_string(),
            alt_loc: None,
            residue_name: "ASP".to_string(),
            chain_id: "A".to_string(),
            residue_seq: 1,
            x: 1.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: "C".to_string(),
        });

        // Add connectivity
        structure.connects.push(Conect {
            atom_serial: 1,
            bonded_atoms: vec![2],
        });

        let connected = structure.get_connected_atoms(1);
        assert_eq!(connected.len(), 1);
        assert_eq!(connected[0].serial, 2);
        assert_eq!(connected[0].name, "CA");
    }
}
