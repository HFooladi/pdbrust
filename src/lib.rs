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

#[derive(Debug, Default)]
pub struct PdbStructure {
    pub atoms: Vec<Atom>,
    pub header: Option<String>,
    pub title: Option<String>,
}

impl PdbStructure {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, PdbError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut structure = PdbStructure::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                let atom = Self::parse_atom_record(&line)?;
                structure.atoms.push(atom);
            } else if line.starts_with("HEADER") {
                structure.header = Some(line[6..].trim().to_string());
            } else if line.starts_with("TITLE ") {
                structure.title = Some(line[6..].trim().to_string());
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
}
