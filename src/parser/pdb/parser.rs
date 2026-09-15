use crate::core::PdbStructure;
use crate::error::PdbError;
use crate::records::{Atom, Conect, Model, Remark, SSBond, SeqRes};
use crate::utils::{decode_hybrid36, parse_float, parse_int};
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parses a PDB file and returns a PdbStructure.
pub fn parse_pdb_file<P: AsRef<Path>>(path: P) -> Result<PdbStructure, PdbError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    parse_pdb_reader(reader)
}

/// Parses PDB data from a string and returns a PdbStructure.
pub fn parse_pdb_string(content: &str) -> Result<PdbStructure, PdbError> {
    let reader = BufReader::new(content.as_bytes());
    parse_pdb_reader(reader)
}

/// Parses PDB data from any reader implementing BufRead.
///
/// This function is useful for parsing PDB data from sources other than files,
/// such as network streams, compressed files, or embedded data.
///
/// Bytes that are not valid UTF-8 (e.g. Latin-1 text in older files) are
/// replaced with U+FFFD instead of failing the whole parse.
///
/// # Arguments
/// * `reader` - Any type implementing BufRead containing PDB data
///
/// # Returns
/// * `Result<PdbStructure, PdbError>` - The parsed structure or an error
pub fn parse_pdb_reader<R: BufRead>(mut reader: R) -> Result<PdbStructure, PdbError> {
    let mut structure = PdbStructure::new();

    // Track current model if parsing multi-model file
    let mut current_model: Option<Model> = None;

    let mut buffer = Vec::new();
    loop {
        buffer.clear();
        if reader.read_until(b'\n', &mut buffer)? == 0 {
            break;
        }
        while matches!(buffer.last(), Some(b'\n' | b'\r')) {
            buffer.pop();
        }
        let text = String::from_utf8_lossy(&buffer);
        let line = Line::new(&text);

        if line.len() < 6 {
            continue; // Skip short lines
        }

        let record_type = line.col(0, 6);

        match record_type.trim() {
            "ATOM" | "HETATM" => {
                let is_hetatm = record_type.trim() == "HETATM";
                let atom = parse_atom_record(&line, is_hetatm)?;

                // Add to current model if we're in a model, otherwise add to structure
                if let Some(model) = &mut current_model {
                    model.atoms.push(atom.clone());
                }

                structure.atoms.push(atom);
            }
            "SEQRES" => {
                structure.seqres.push(parse_seqres_record(&line)?);
            }
            "CONECT" => {
                if let Some(conect) = parse_conect_record(&line)? {
                    structure.connects.push(conect);
                }
            }
            "SSBOND" => {
                structure.ssbonds.push(parse_ssbond_record(&line)?);
            }
            "REMARK" => {
                let remark = parse_remark_record(&line);

                // Add to current model if we're in a model, otherwise add to structure
                if let Some(model) = &mut current_model {
                    model.remarks.push(remark.clone());
                }

                structure.remarks.push(remark);
            }
            "HEADER" => {
                structure.header = Some(line.rest(10).trim().to_string());
            }
            "TITLE" => {
                // Continuation lines (serial number in columns 9-10) extend the title.
                let text = line.rest(10);
                let text = text.trim();
                structure.title = Some(match structure.title.take() {
                    Some(previous) if !previous.is_empty() && !text.is_empty() => {
                        format!("{previous} {text}")
                    }
                    Some(previous) if !previous.is_empty() => previous,
                    _ => text.to_string(),
                });
            }
            "MODEL" => {
                // Finish previous model if there was one
                if let Some(model) = current_model.take() {
                    structure.models.push(model);
                }

                // The model number belongs in columns 11-14; some programs write
                // it right after the record name ("MODEL 1"). Without a number,
                // models are numbered sequentially.
                let serial = parse_int(&line.col(10, 14))
                    .or_else(|_| parse_int(&line.rest(6)))
                    .unwrap_or(structure.models.len() as i32 + 1);

                current_model = Some(Model {
                    serial,
                    atoms: Vec::new(),
                    remarks: Vec::new(),
                });
                structure.current_model = Some(serial);
            }
            "ENDMDL" => {
                // Finish current model if there is one
                if let Some(model) = current_model.take() {
                    structure.models.push(model);
                }
                structure.current_model = None;
            }
            "END" => {
                // End of file, finish current model if there is one
                if let Some(model) = current_model.take() {
                    structure.models.push(model);
                }
                break;
            }
            _ => {
                // Ignore other record types
            }
        }
    }

    // If we have a model that wasn't closed with ENDMDL
    if let Some(model) = current_model {
        structure.models.push(model);
    }

    Ok(structure)
}

/// A line of a fixed-column PDB file, addressed by 0-based column positions.
///
/// Columns count characters. Lines containing non-ASCII characters are handled
/// by character position, so they neither panic nor shift the fields; ASCII
/// lines (almost all real PDB lines) are sliced directly.
struct Line<'a> {
    text: &'a str,
    ascii: bool,
    len: usize,
}

impl<'a> Line<'a> {
    fn new(text: &'a str) -> Self {
        let ascii = text.is_ascii();
        let len = if ascii {
            text.len()
        } else {
            text.chars().count()
        };
        Self { text, ascii, len }
    }

    /// Number of columns (characters) in the line.
    fn len(&self) -> usize {
        self.len
    }

    /// Columns `start..end`, shortened (possibly to "") when the line is shorter.
    fn col(&self, start: usize, end: usize) -> Cow<'a, str> {
        if self.ascii {
            let end = end.min(self.text.len());
            Cow::Borrowed(self.text.get(start..end).unwrap_or(""))
        } else {
            Cow::Owned(
                self.text
                    .chars()
                    .skip(start)
                    .take(end.saturating_sub(start))
                    .collect(),
            )
        }
    }

    /// Everything from column `start` to the end of the line.
    fn rest(&self, start: usize) -> Cow<'a, str> {
        self.col(start, usize::MAX)
    }

    /// The character in column `index`, or `None` if blank or beyond the line.
    fn non_blank_char(&self, index: usize) -> Option<char> {
        let c = if self.ascii {
            self.text.as_bytes().get(index).map(|&b| b as char)
        } else {
            self.text.chars().nth(index)
        };
        c.filter(|&c| c != ' ')
    }
}

/// Parses an optional floating-point field; a blank field yields `default`.
fn parse_float_or(field: &str, default: f64) -> Result<f64, PdbError> {
    if field.trim().is_empty() {
        Ok(default)
    } else {
        Ok(parse_float(field)?)
    }
}

/// Parses an optional integer field; a blank field yields `default`.
fn parse_int_or(field: &str, default: i32) -> Result<i32, PdbError> {
    if field.trim().is_empty() {
        Ok(default)
    } else {
        Ok(parse_int(field)?)
    }
}

/// Parses an ATOM or HETATM record.
fn parse_atom_record(line: &Line, is_hetatm: bool) -> Result<Atom, PdbError> {
    if line.len() < 54 {
        return Err(PdbError::InvalidRecord(
            "ATOM/HETATM record too short".to_string(),
        ));
    }

    let serial = decode_hybrid36(&line.col(6, 11), 5)?;
    let name = line.col(12, 16).trim().to_string();
    let alt_loc = line.non_blank_char(16);
    let residue_name = line.col(17, 20).trim().to_string();
    let chain_id = line.col(21, 22).into_owned();
    let residue_seq = decode_hybrid36(&line.col(22, 26), 4)?;
    let ins_code = line.non_blank_char(26);

    let x = parse_float(&line.col(30, 38))?;
    let y = parse_float(&line.col(38, 46))?;
    let z = parse_float(&line.col(46, 54))?;

    // Occupancy and B-factor are optional; blank or missing fields get defaults.
    let occupancy = if line.len() >= 60 {
        parse_float_or(&line.col(54, 60), 1.0)?
    } else {
        1.0
    };
    let temp_factor = if line.len() >= 66 {
        parse_float_or(&line.col(60, 66), 0.0)?
    } else {
        0.0
    };

    let element = if line.len() >= 78 {
        line.col(76, 78).trim().to_string()
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
        ins_code,
        is_hetatm,
    })
}

/// Parses a SEQRES record.
fn parse_seqres_record(line: &Line) -> Result<SeqRes, PdbError> {
    if line.len() < 19 {
        return Err(PdbError::InvalidRecord(
            "SEQRES record too short".to_string(),
        ));
    }

    let serial = parse_int(&line.col(7, 10))?;
    let chain_id = line.col(11, 12).into_owned();
    let num_residues = parse_int(&line.col(13, 17))?;

    // Residue names are in 4-column fields starting at column 20.
    let section = line.rest(19);
    let section = Line::new(&section);
    let residues = (0..section.len())
        .step_by(4)
        .filter(|&i| i + 3 <= section.len())
        .map(|i| section.col(i, i + 3).trim().to_string())
        .filter(|residue| !residue.is_empty())
        .collect();

    Ok(SeqRes {
        serial,
        chain_id,
        num_residues,
        residues,
    })
}

/// Parses a CONECT record. Returns `None` for a record that lists no bonded atoms.
fn parse_conect_record(line: &Line) -> Result<Option<Conect>, PdbError> {
    if line.len() < 11 {
        return Err(PdbError::InvalidRecord(
            "CONECT record too short".to_string(),
        ));
    }

    let atom1 = decode_hybrid36(&line.col(6, 11), 5)?;

    let atom2_field = line.col(11, 16);
    if line.len() < 16 || atom2_field.trim().is_empty() {
        return Ok(None);
    }
    let atom2 = decode_hybrid36(&atom2_field, 5)?;

    // Further bonded atoms are read only from complete 5-column fields.
    let optional_atom = |start: usize| -> Result<Option<i32>, PdbError> {
        let field = line.col(start, start + 5);
        if line.len() < start + 5 || field.trim().is_empty() {
            Ok(None)
        } else {
            decode_hybrid36(&field, 5).map(Some)
        }
    };
    let atom3 = optional_atom(16)?;
    let atom4 = optional_atom(21)?;

    Ok(Some(Conect {
        atom1,
        atom2,
        atom3,
        atom4,
    }))
}

/// Parses an SSBOND record.
fn parse_ssbond_record(line: &Line) -> Result<SSBond, PdbError> {
    if line.len() < 35 {
        // Changed from 70 to 35 to handle shorter SSBOND records
        return Err(PdbError::InvalidRecord(
            "SSBOND record too short".to_string(),
        ));
    }

    let serial = parse_int(&line.col(7, 10))?;
    let residue1_name = line.col(11, 14).trim().to_string();
    let chain1_id = line.col(15, 16).into_owned();
    let residue1_seq = parse_int(&line.col(17, 21))?;
    let icode1 = line.non_blank_char(21);

    let residue2_name = line.col(25, 28).trim().to_string();
    let chain2_id = line.col(29, 30).into_owned();
    let residue2_seq = parse_int(&line.col(31, 35))?;
    let icode2 = line.non_blank_char(35);

    // Symmetry operators (columns 60-65, 67-72) and bond length (74-78) are
    // optional; missing or blank fields get the identity operator and a
    // typical disulfide length. Shorter lines (some programs put other text
    // in these columns) keep the defaults.
    let (sym1, sym2, length) = if line.len() >= 72 {
        (
            parse_int_or(&line.col(59, 65), 1555)?,
            parse_int_or(&line.col(66, 72), 1555)?,
            if line.len() >= 78 {
                parse_float_or(&line.col(73, 78), 2.04)?
            } else {
                2.04
            },
        )
    } else {
        (1555, 1555, 2.04)
    };

    Ok(SSBond {
        serial,
        residue1_name,
        chain1_id,
        residue1_seq,
        icode1,
        residue2_name,
        chain2_id,
        residue2_seq,
        icode2,
        sym1,
        sym2,
        length,
    })
}

/// Parses a REMARK record.
///
/// Standard remarks carry a number in columns 8-10 and text from column 12.
/// Unnumbered free-text remarks (e.g. from GROMACS) are kept with number 0.
fn parse_remark_record(line: &Line) -> Remark {
    match parse_int(&line.col(6, 10)) {
        Ok(number) => Remark {
            number,
            content: line.rest(11).trim().to_string(),
        },
        Err(_) => Remark {
            number: 0,
            content: line.rest(6).trim().to_string(),
        },
    }
}
