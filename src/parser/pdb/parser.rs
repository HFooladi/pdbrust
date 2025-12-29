use crate::core::PdbStructure;
use crate::error::PdbError;
use crate::records::{Atom, Conect, Model, Remark, SSBond, SeqRes};
use crate::utils::{parse_float, parse_int};
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

/// Internal function to parse PDB from any reader.
fn parse_pdb_reader<R: BufRead>(reader: R) -> Result<PdbStructure, PdbError> {
    let mut structure = PdbStructure::new();

    // Track current model if parsing multi-model file
    let mut current_model: Option<Model> = None;

    for line in reader.lines() {
        let line = line?;

        if line.len() < 6 {
            continue; // Skip short lines
        }

        let record_type = line[0..6].trim();

        match record_type {
            "ATOM" | "HETATM" => {
                let atom = parse_atom_record(&line)?;

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
                structure.connects.push(parse_conect_record(&line)?);
            }
            "SSBOND" => {
                structure.ssbonds.push(parse_ssbond_record(&line)?);
            }
            "REMARK" => {
                let remark = parse_remark_record(&line)?;

                // Add to current model if we're in a model, otherwise add to structure
                if let Some(model) = &mut current_model {
                    model.remarks.push(remark.clone());
                }

                structure.remarks.push(remark);
            }
            "HEADER" => {
                structure.header = Some(line[10..].trim().to_string());
            }
            "TITLE" => {
                structure.title = Some(line[10..].trim().to_string());
            }
            "MODEL" => {
                // Finish previous model if there was one
                if let Some(model) = current_model.take() {
                    structure.models.push(model);
                }

                // Start a new model
                let serial = if line.len() >= 14 {
                    let serial_str = line[10..14].trim();
                    if serial_str.is_empty() {
                        1 // Default model number if empty
                    } else {
                        parse_int(serial_str)?
                    }
                } else {
                    1 // Default model number if line is too short
                };

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

/// Parses an ATOM or HETATM record.
fn parse_atom_record(line: &str) -> Result<Atom, PdbError> {
    if line.len() < 54 {
        return Err(PdbError::InvalidRecord(
            "ATOM/HETATM record too short".to_string(),
        ));
    }

    let serial = parse_int(&line[6..11])?;
    let name = line[12..16].trim().to_string();

    let alt_loc = if line.len() > 16 {
        let c = line.chars().nth(16).unwrap_or(' ');
        if c == ' ' {
            None
        } else {
            Some(c)
        }
    } else {
        None
    };

    let residue_name = line[17..20].trim().to_string();
    let chain_id = line[21..22].to_string();
    let residue_seq = parse_int(&line[22..26])?;

    let ins_code = if line.len() > 26 {
        let c = line.chars().nth(26).unwrap_or(' ');
        if c == ' ' {
            None
        } else {
            Some(c)
        }
    } else {
        None
    };

    let x = parse_float(&line[30..38])?;
    let y = parse_float(&line[38..46])?;
    let z = parse_float(&line[46..54])?;

    let occupancy = if line.len() >= 60 {
        parse_float(&line[54..60])?
    } else {
        1.0
    };
    let temp_factor = if line.len() >= 66 {
        parse_float(&line[60..66])?
    } else {
        0.0
    };

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
        ins_code,
    })
}

/// Parses a SEQRES record.
fn parse_seqres_record(line: &str) -> Result<SeqRes, PdbError> {
    if line.len() < 19 {
        return Err(PdbError::InvalidRecord(
            "SEQRES record too short".to_string(),
        ));
    }

    let serial = parse_int(&line[8..10])?;
    let chain_id = line[11..12].to_string();
    let num_residues = parse_int(&line[13..17])?;

    let mut residues = Vec::new();
    let residue_section = if line.len() > 19 { &line[19..] } else { "" };

    // Parse residues (each takes 4 characters, with the residue name being 3 chars)
    for i in (0..residue_section.len()).step_by(4) {
        if i + 3 <= residue_section.len() {
            let residue = residue_section[i..i + 3].trim();
            if !residue.is_empty() {
                residues.push(residue.to_string());
            }
        }
    }

    Ok(SeqRes {
        serial,
        chain_id,
        num_residues,
        residues,
    })
}

/// Parses a CONECT record.
fn parse_conect_record(line: &str) -> Result<Conect, PdbError> {
    if line.len() < 11 {
        return Err(PdbError::InvalidRecord(
            "CONECT record too short".to_string(),
        ));
    }

    let atom1 = parse_int(&line[6..11])?;

    // Parse atom2 (required)
    let atom2 = if line.len() >= 16 {
        parse_int(&line[11..16])?
    } else {
        return Err(PdbError::InvalidRecord(
            "CONECT record missing second atom".to_string(),
        ));
    };

    // Parse atom3 (optional)
    let atom3 = if line.len() >= 21 && !line[16..21].trim().is_empty() {
        Some(parse_int(&line[16..21])?)
    } else {
        None
    };

    // Parse atom4 (optional)
    let atom4 = if line.len() >= 26 && !line[21..26].trim().is_empty() {
        Some(parse_int(&line[21..26])?)
    } else {
        None
    };

    Ok(Conect {
        atom1,
        atom2,
        atom3,
        atom4,
    })
}

/// Parses an SSBOND record.
fn parse_ssbond_record(line: &str) -> Result<SSBond, PdbError> {
    if line.len() < 35 {
        // Changed from 70 to 35 to handle shorter SSBOND records
        return Err(PdbError::InvalidRecord(
            "SSBOND record too short".to_string(),
        ));
    }

    let serial = parse_int(&line[7..10])?;
    let residue1_name = line[11..14].trim().to_string();
    let chain1_id = line[15..16].to_string();
    let residue1_seq = parse_int(&line[17..21])?;

    let icode1 = if line.len() > 21 {
        let c = line.chars().nth(21).unwrap_or(' ');
        if c == ' ' {
            None
        } else {
            Some(c)
        }
    } else {
        None
    };

    let residue2_name = line[25..28].trim().to_string();
    let chain2_id = line[29..30].to_string();
    let residue2_seq = parse_int(&line[31..35])?;

    let icode2 = if line.len() > 35 {
        let c = line.chars().nth(35).unwrap_or(' ');
        if c == ' ' {
            None
        } else {
            Some(c)
        }
    } else {
        None
    };

    // Parse symmetry operators and length if available
    let (sym1, sym2, length) = if line.len() >= 70 {
        (
            parse_int(&line[59..63])?,
            parse_int(&line[66..70])?,
            if line.len() >= 78 {
                parse_float(&line[73..78])?
            } else {
                2.04
            }, // Default to 2.04 if length not specified
        )
    } else {
        (1555, 1555, 2.04) // Default values if not specified
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
fn parse_remark_record(line: &str) -> Result<Remark, PdbError> {
    if line.len() < 10 {
        return Err(PdbError::InvalidRecord(
            "REMARK record too short".to_string(),
        ));
    }

    let number = parse_int(&line[6..10])?;
    let content = line[11..].trim().to_string();

    Ok(Remark { number, content })
}
