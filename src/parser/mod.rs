use std::path::Path;
use std::io::{self, BufRead, BufReader};
use crate::error::PdbError;
use crate::core::PdbStructure;
use crate::records::{Atom, SeqRes, Conect, SSBond};
use crate::utils::{parse_float, parse_int};

/// Functions for parsing different PDB record types
pub struct Parser;

impl Parser {
    /// Parses a PDB file and returns a PdbStructure.
    pub fn parse_pdb_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
        let file = std::fs::File::open(path)?;
        let reader = BufReader::new(file);
        let mut structure = PdbStructure::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
                structure.atoms.push(parse_atom_record(&line)?);
            } else if line.starts_with("SEQRES") {
                structure.seqres.push(parse_seqres_record(&line)?);
            } else if line.starts_with("CONECT") {
                structure.conect.push(parse_conect_record(&line)?);
            } else if line.starts_with("SSBOND") {
                structure.ssbond.push(parse_ssbond_record(&line)?);
            } else if line.starts_with("REMARK") {
                structure.remarks.push(line);
            } else if line.starts_with("HEADER") {
                structure.header = line;
            } else if line.starts_with("TITLE ") {
                structure.title = line;
            }
        }

        Ok(structure)
    }

    /// Parses an ATOM or HETATM record.
    fn parse_atom_record(line: &str) -> Result<Atom, PdbError> {
        if line.len() < 80 {
            return Err(PdbError::InvalidRecord("ATOM/HETATM record too short".to_string()));
        }

        Ok(Atom {
            serial: parse_int(&line[6..11])?,
            name: line[12..16].trim().to_string(),
            alt_loc: line[16..17].chars().next().unwrap_or(' '),
            res_name: line[17..20].trim().to_string(),
            chain_id: line[21..22].chars().next().unwrap_or(' '),
            res_seq: parse_int(&line[22..26])?,
            i_code: line[26..27].chars().next().unwrap_or(' '),
            x: parse_float(&line[30..38])?,
            y: parse_float(&line[38..46])?,
            z: parse_float(&line[46..54])?,
            occupancy: parse_float(&line[54..60])?,
            temp_factor: parse_float(&line[60..66])?,
            segment: line[72..76].trim().to_string(),
            element: line[76..78].trim().to_string(),
            charge: line[78..80].trim().to_string(),
        })
    }

    /// Parses a SEQRES record.
    fn parse_seqres_record(line: &str) -> Result<SeqRes, PdbError> {
        if line.len() < 80 {
            return Err(PdbError::InvalidRecord("SEQRES record too short".to_string()));
        }

        let mut res_names = Vec::new();
        for i in (19..70).step_by(4) {
            let res_name = line[i..i + 3].trim().to_string();
            if !res_name.is_empty() {
                res_names.push(res_name);
            }
        }

        Ok(SeqRes {
            ser_num: parse_int(&line[8..10])?,
            chain_id: line[11..12].chars().next().unwrap_or(' '),
            num_res: parse_int(&line[13..17])?,
            res_names,
        })
    }

    /// Parses a CONECT record.
    fn parse_conect_record(line: &str) -> Result<Conect, PdbError> {
        if line.len() < 31 {
            return Err(PdbError::InvalidRecord("CONECT record too short".to_string()));
        }

        Ok(Conect {
            serial1: parse_int(&line[6..11])?,
            serial2: parse_int(&line[11..16])?,
            serial3: if line[16..21].trim().is_empty() {
                None
            } else {
                Some(parse_int(&line[16..21])?)
            },
            serial4: if line[21..26].trim().is_empty() {
                None
            } else {
                Some(parse_int(&line[21..26])?)
            },
        })
    }

    /// Parses an SSBOND record.
    ///
    /// This method parses an SSBOND record from a PDB file line.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the SSBOND record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `SSBond` or a `PdbError`.
    ///
    /// # Errors
    fn parse_ssbond_record(line: &str) -> Result<SSBond, PdbError> {
        if line.len() < 80 {
            return Err(PdbError::InvalidRecord("SSBOND record too short".to_string()));
        }

        Ok(SSBond {
            ser_num: parse_int(&line[8..10])?,
            res1: line[11..14].trim().to_string(),
            chain_id1: line[15..16].chars().next().unwrap_or(' '),
            seq1: parse_int(&line[17..21])?,
            icode1: line[21..22].chars().next().unwrap_or(' '),
            res2: line[25..28].trim().to_string(),
            chain_id2: line[29..30].chars().next().unwrap_or(' '),
            seq2: parse_int(&line[31..35])?,
            icode2: line[35..36].chars().next().unwrap_or(' '),
            sym1: parse_int(&line[59..65])?,
            sym2: parse_int(&line[66..72])?,
            length: parse_float(&line[73..78])?,
        })
    }
    /// Parses a REMARK record from a PDB file line.
    ///
    /// This method extracts the remark number and content from a REMARK record.
    ///
    /// # Arguments
    ///
    /// * `line` - A string slice containing the REMARK record line
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing either the parsed `Remark` or a `PdbError`.
    ///
    /// # Errors
    ///
    /// Returns `PdbError::InvalidRecord` if the line format is invalid or
    /// `PdbError::ParseError` if the remark number cannot be parsed.
    pub fn parse_remark_record(line: &str) -> Result<Remark, PdbError> {
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

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_pdb(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file
    }

    #[test]
    fn test_parse_empty_file() {
        let file = create_test_pdb("");
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert!(structure.atoms.is_empty());
        assert!(structure.conect.is_empty());
        assert!(structure.ssbond.is_empty());
        assert!(structure.seqres.is_empty());
    }

    #[test]
    fn test_parse_single_atom() {
        let content = "ATOM      1  N   ASP A  30      27.360  44.310  48.500  1.00 13.79           N\n";
        let file = create_test_pdb(content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.atoms.len(), 1);
        let atom = &structure.atoms[0];
        assert_eq!(atom.serial, 1);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.res_name, "ASP");
        assert_eq!(atom.chain_id, 'A');
        assert_eq!(atom.res_seq, 30);
        assert_eq!(atom.alt_loc, ' ');
        assert_eq!(atom.i_code, ' ');
        assert_eq!(atom.occupancy, 1.00);
        assert_eq!(atom.temp_factor, 13.79);
        assert_eq!(atom.element, "N");
        assert_eq!(atom.segment, "");
        assert_eq!(atom.charge, "");
    }

    #[test]
    fn test_parse_seqres() {
        let content = "SEQRES   1 A   20  SER ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA\n";
        let file = create_test_pdb(content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.seqres.len(), 1);
        let seqres = &structure.seqres[0];
        assert_eq!(seqres.ser_num, 1);
        assert_eq!(seqres.chain_id, 'A');
        assert_eq!(seqres.num_res, 20);
        let expected_res_names = vec!["SER", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA"];
        assert_eq!(seqres.res_names, expected_res_names);
    }

    #[test]
    fn test_parse_conect() {
        let content = "CONECT 1179 1174 1177 1180\n";
        let file = create_test_pdb(content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.conect.len(), 1);
        let conect = &structure.conect[0];
        assert_eq!(conect.serial1, 1179);
        assert_eq!(conect.serial2, 1174);
        assert_eq!(conect.serial3, Some(1177));
        assert_eq!(conect.serial4, Some(1180));
    }

    #[test]
    fn test_parse_ssbond() {
        let content = "SSBOND   1 CYS A   85    CYS A  101    1555  1555  2.03\n";
        let file = create_test_pdb(content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.ssbond.len(), 1);
        let ssbond = &structure.ssbond[0];
        assert_eq!(ssbond.ser_num, 1);
        assert_eq!(ssbond.res1, "CYS");
        assert_eq!(ssbond.chain_id1, 'A');
        assert_eq!(ssbond.seq1, 85);
        assert_eq!(ssbond.icode1, ' ');
        assert_eq!(ssbond.res2, "CYS");
        assert_eq!(ssbond.chain_id2, 'A');
        assert_eq!(ssbond.seq2, 101);
        assert_eq!(ssbond.icode2, ' ');
        assert_eq!(ssbond.sym1, 1555);
        assert_eq!(ssbond.sym2, 1555);
        assert!((ssbond.length - 2.03).abs() < 1e-6);
    }

    #[test]
    fn test_parse_complete_structure() {
        let content = "\
HEADER    PROTEIN                                01-JAN-01   1ABC
TITLE     TEST PROTEIN
REMARK   1 THIS IS A TEST
SEQRES   1 A   20  SER ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA
ATOM      1  N   ASP A  30      27.360  44.310  48.500  1.00 13.79           N
ATOM      2  CA  ASP A  30      26.360  43.310  48.500  1.00 13.79           C
CONECT    1    2
SSBOND   1 CYS A   85    CYS A  101    1555  1555  2.03
END
";
        let file = create_test_pdb(content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.header, "PROTEIN");
        assert_eq!(structure.title, "TEST PROTEIN");
        assert_eq!(structure.remarks.len(), 1);
        assert_eq!(structure.remarks[0], "THIS IS A TEST");
        assert_eq!(structure.seqres.len(), 1);
        assert_eq!(structure.atoms.len(), 2);
        assert_eq!(structure.conect.len(), 1);
        assert_eq!(structure.ssbond.len(), 1);
    }
}