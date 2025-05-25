//! Converts parsed mmCIF data into PdbStructure format

use crate::core::mmcif::{MmcifParser, Category};
use crate::core::PdbStructure;
use crate::records::{Atom, SeqRes, SSBond, Remark};
use crate::error::PdbError;
use std::collections::HashMap;

/// Converts mmCIF data to PdbStructure
pub fn mmcif_to_pdb_structure(parser: &MmcifParser) -> Result<PdbStructure, PdbError> {
    let mut structure = PdbStructure::new();

    // Parse header and title information
    parse_header_info(parser, &mut structure)?;
    
    // Parse atoms from _atom_site category
    parse_atoms(parser, &mut structure)?;
    
    // Parse sequence information from _entity_poly_seq
    parse_sequences(parser, &mut structure)?;
    
    // Parse disulfide bonds from _struct_disulfid
    parse_disulfide_bonds(parser, &mut structure)?;
    
    // Parse connectivity from _struct_conn
    parse_connectivity(parser, &mut structure)?;

    // Parse remarks/comments
    parse_remarks(parser, &mut structure)?;

    Ok(structure)
}

/// Parse header and title information
fn parse_header_info(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    // Parse entry information
    if let Some(entry_category) = parser.get_category("entry") {
        if let Some(row) = entry_category.get_row(0) {
            if let Some(id) = row.get("id") {
                structure.header = Some(format!("mmCIF entry: {}", id));
            }
        }
    }

    // Parse structure title
    if let Some(struct_category) = parser.get_category("struct") {
        if let Some(row) = struct_category.get_row(0) {
            if let Some(title) = row.get("title") {
                structure.title = Some(title.to_string());
            }
        }
    }

    // Parse database references for additional header info
    if let Some(database_category) = parser.get_category("database_PDB_rev") {
        if let Some(row) = database_category.get_row(0) {
            let mut header_parts = Vec::new();
            
            if let Some(date) = row.get("date_original") {
                header_parts.push(format!("Original date: {}", date));
            }
            if let Some(id) = row.get("pdb_id") {
                header_parts.push(format!("PDB ID: {}", id));
            }
            
            if !header_parts.is_empty() {
                let existing_header = structure.header.as_deref().unwrap_or("");
                structure.header = Some(format!("{} | {}", existing_header, header_parts.join(", ")));
            }
        }
    }

    Ok(())
}

/// Parse atoms from _atom_site category
fn parse_atoms(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    let atom_site = match parser.get_category("atom_site") {
        Some(category) => category,
        None => return Ok(()), // No atoms found, not an error
    };

    // Get column indices for efficient access
    let col_indices = get_atom_column_indices(atom_site)?;

    for (row_idx, row) in atom_site.rows.iter().enumerate() {
        let atom = parse_atom_row(row, &col_indices)?;
        structure.atoms.push(atom);
    }

    Ok(())
}

/// Helper struct to store column indices for atom parsing
struct AtomColumnIndices {
    id: usize,
    type_symbol: usize,
    label_atom_id: usize,
    label_alt_id: Option<usize>,
    label_comp_id: usize,
    label_asym_id: usize,
    label_seq_id: usize,
    pdbx_pdb_ins_code: Option<usize>,
    cartn_x: usize,
    cartn_y: usize,
    cartn_z: usize,
    occupancy: Option<usize>,
    b_iso_or_equiv: Option<usize>,
    group_pdb: usize,
}

fn get_atom_column_indices(category: &Category) -> Result<AtomColumnIndices, PdbError> {
    let headers = &category.headers;
    
    let find_required = |name: &str| -> Result<usize, PdbError> {
        headers.iter().position(|h| h == name)
            .ok_or_else(|| PdbError::InvalidRecord(format!("Required column '{}' not found in _atom_site", name)))
    };
    
    let find_optional = |name: &str| -> Option<usize> {
        headers.iter().position(|h| h == name)
    };

    Ok(AtomColumnIndices {
        id: find_required("id")?,
        type_symbol: find_required("type_symbol")?,
        label_atom_id: find_required("label_atom_id")?,
        label_alt_id: find_optional("label_alt_id"),
        label_comp_id: find_required("label_comp_id")?,
        label_asym_id: find_required("label_asym_id")?,
        label_seq_id: find_required("label_seq_id")?,
        pdbx_pdb_ins_code: find_optional("pdbx_PDB_ins_code"),
        cartn_x: find_required("Cartn_x")?,
        cartn_y: find_required("Cartn_y")?,
        cartn_z: find_required("Cartn_z")?,
        occupancy: find_optional("occupancy"),
        b_iso_or_equiv: find_optional("B_iso_or_equiv"),
        group_pdb: find_required("group_PDB")?,
    })
}

fn parse_atom_row(row: &[String], indices: &AtomColumnIndices) -> Result<Atom, PdbError> {
    let serial: i32 = row[indices.id].parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid atom serial: {}", row[indices.id])))?;
    
    let name = row[indices.label_atom_id].clone();
    
    let alt_loc = if let Some(idx) = indices.label_alt_id {
        let alt_str = &row[idx];
        if alt_str == "." || alt_str == "?" || alt_str.is_empty() {
            None
        } else {
            alt_str.chars().next()
        }
    } else {
        None
    };
    
    let residue_name = row[indices.label_comp_id].clone();
    let chain_id = row[indices.label_asym_id].clone();
    
    let residue_seq: i32 = row[indices.label_seq_id].parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid residue sequence: {}", row[indices.label_seq_id])))?;
    
    let ins_code = if let Some(idx) = indices.pdbx_pdb_ins_code {
        let ins_str = &row[idx];
        if ins_str == "." || ins_str == "?" || ins_str.is_empty() {
            None
        } else {
            ins_str.chars().next()
        }
    } else {
        None
    };
    
    let x: f64 = row[indices.cartn_x].parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid x coordinate: {}", row[indices.cartn_x])))?;
    
    let y: f64 = row[indices.cartn_y].parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid y coordinate: {}", row[indices.cartn_y])))?;
    
    let z: f64 = row[indices.cartn_z].parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid z coordinate: {}", row[indices.cartn_z])))?;
    
    let occupancy = if let Some(idx) = indices.occupancy {
        row[idx].parse().unwrap_or(1.0)
    } else {
        1.0
    };
    
    let temp_factor = if let Some(idx) = indices.b_iso_or_equiv {
        row[idx].parse().unwrap_or(0.0)
    } else {
        0.0
    };
    
    let element = row[indices.type_symbol].clone();

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

/// Parse sequence information from _entity_poly_seq
fn parse_sequences(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    let entity_poly_seq = match parser.get_category("entity_poly_seq") {
        Some(category) => category,
        None => return Ok(()), // No sequence information, not an error
    };

    // Group sequences by entity_id and build SEQRES records
    let mut sequences: HashMap<String, Vec<String>> = HashMap::new();
    
    if let (Some(entity_col), Some(mon_id_col), Some(num_col)) = (
        entity_poly_seq.get_column("entity_id"),
        entity_poly_seq.get_column("mon_id"),
        entity_poly_seq.get_column("num")
    ) {
        for ((entity_id, mon_id), num) in entity_col.iter().zip(mon_id_col.iter()).zip(num_col.iter()) {
            let seq_entry = sequences.entry(entity_id.to_string()).or_insert_with(Vec::new);
            
            // Parse the sequence number to ensure proper ordering
            let seq_num: usize = num.parse().unwrap_or(0);
            
            // Extend vector if necessary
            if seq_entry.len() <= seq_num {
                seq_entry.resize(seq_num + 1, String::new());
            }
            
            if seq_num > 0 {
                seq_entry[seq_num - 1] = mon_id.to_string(); // Convert to 0-based indexing
            }
        }
    }

    // Convert to SEQRES records
    // Map entity IDs to chain IDs using _struct_asym if available
    let chain_mapping = get_entity_to_chain_mapping(parser);
    
    let mut serial = 1;
    for (entity_id, residues) in sequences {
        let residues: Vec<String> = residues.into_iter().filter(|r| !r.is_empty()).collect();
        
        if residues.is_empty() {
            continue;
        }

        // Try to find the corresponding chain ID
        let chain_id = chain_mapping.get(&entity_id)
            .cloned()
            .unwrap_or_else(|| entity_id.clone());

        let seqres = SeqRes {
            serial,
            chain_id,
            num_residues: residues.len() as i32,
            residues,
        };
        
        structure.seqres.push(seqres);
        serial += 1;
    }

    Ok(())
}

/// Get mapping from entity ID to chain ID using _struct_asym
fn get_entity_to_chain_mapping(parser: &MmcifParser) -> HashMap<String, String> {
    let mut mapping = HashMap::new();
    
    if let Some(struct_asym) = parser.get_category("struct_asym") {
        if let (Some(entity_col), Some(id_col)) = (
            struct_asym.get_column("entity_id"),
            struct_asym.get_column("id")
        ) {
            for (entity_id, chain_id) in entity_col.iter().zip(id_col.iter()) {
                mapping.insert(entity_id.to_string(), chain_id.to_string());
            }
        }
    }
    
    mapping
}

/// Parse disulfide bonds from _struct_disulfid
fn parse_disulfide_bonds(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    let struct_disulfid = match parser.get_category("struct_disulfid") {
        Some(category) => category,
        None => return Ok(()), // No disulfide bonds, not an error
    };

    for (serial, row) in struct_disulfid.rows.iter().enumerate() {
        if let Some(ssbond) = parse_disulfide_row(row, &struct_disulfid.headers, serial as i32 + 1)? {
            structure.ssbonds.push(ssbond);
        }
    }

    Ok(())
}

fn parse_disulfide_row(row: &[String], headers: &[String], serial: i32) -> Result<Option<SSBond>, PdbError> {
    let find_col = |name: &str| headers.iter().position(|h| h == name);
    
    let ptnr1_auth_asym_id = find_col("ptnr1_auth_asym_id").and_then(|i| row.get(i));
    let ptnr1_auth_comp_id = find_col("ptnr1_auth_comp_id").and_then(|i| row.get(i));
    let ptnr1_auth_seq_id = find_col("ptnr1_auth_seq_id").and_then(|i| row.get(i));
    
    let ptnr2_auth_asym_id = find_col("ptnr2_auth_asym_id").and_then(|i| row.get(i));
    let ptnr2_auth_comp_id = find_col("ptnr2_auth_comp_id").and_then(|i| row.get(i));
    let ptnr2_auth_seq_id = find_col("ptnr2_auth_seq_id").and_then(|i| row.get(i));

    if let (Some(chain1), Some(res1_name), Some(res1_seq), 
            Some(chain2), Some(res2_name), Some(res2_seq)) = 
        (ptnr1_auth_asym_id, ptnr1_auth_comp_id, ptnr1_auth_seq_id,
         ptnr2_auth_asym_id, ptnr2_auth_comp_id, ptnr2_auth_seq_id) {
        
        let residue1_seq: i32 = res1_seq.parse().map_err(|_| 
            PdbError::ParseError(format!("Invalid residue sequence: {}", res1_seq)))?;
        let residue2_seq: i32 = res2_seq.parse().map_err(|_| 
            PdbError::ParseError(format!("Invalid residue sequence: {}", res2_seq)))?;

        // Try to get distance if available
        let length = find_col("ptnr1_symmetry")
            .and_then(|i| row.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(2.04); // Default S-S bond length

        Ok(Some(SSBond {
            serial,
            residue1_name: res1_name.clone(),
            chain1_id: chain1.clone(),
            residue1_seq,
            icode1: None, // mmCIF doesn't typically use insertion codes the same way
            residue2_name: res2_name.clone(),
            chain2_id: chain2.clone(),
            residue2_seq,
            icode2: None,
            sym1: 1555, // Default symmetry
            sym2: 1555,
            length,
        }))
    } else {
        Ok(None)
    }
}

/// Parse connectivity information from _struct_conn
fn parse_connectivity(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    let struct_conn = match parser.get_category("struct_conn") {
        Some(category) => category,
        None => return Ok(()), // No connectivity info, not an error
    };

    // This is more complex in mmCIF as connectivity can be of different types
    // For now, we'll focus on covalent bonds that can be converted to CONECT records
    // In a full implementation, you might want to filter by conn_type_id

    // Note: mmCIF connectivity is more complex than PDB CONECT records
    // This is a simplified conversion
    for row in &struct_conn.rows {
        // This would need more sophisticated parsing based on the specific
        // mmCIF connectivity format requirements
        // Left as a placeholder for now as it requires more detailed mmCIF specification knowledge
    }

    Ok(())
}

/// Parse remarks and comments
fn parse_remarks(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    // Parse from various comment categories
    if let Some(audit_conform) = parser.get_category("audit_conform") {
        for (i, row) in audit_conform.rows.iter().enumerate() {
            if let Some(row_map) = audit_conform.get_row(i) {
                if let Some(dict_name) = row_map.get("dict_name") {
                    let remark = Remark {
                        number: 1,
                        content: format!("Conforms to dictionary: {}", dict_name),
                    };
                    structure.remarks.push(remark);
                }
            }
        }
    }

    // Parse resolution information
    if let Some(refine) = parser.get_category("refine") {
        if let Some(row) = refine.get_row(0) {
            if let Some(resolution) = row.get("ls_d_res_high") {
                let remark = Remark {
                    number: 2,
                    content: format!("RESOLUTION. {} ANGSTROMS.", resolution),
                };
                structure.remarks.push(remark);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::mmcif::MmcifParser;
    use std::io::Cursor;

    #[test]
    fn test_mmcif_atom_parsing() {
        let mmcif_data = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N . MET A 1 20.154 10.000 5.000 1.00 25.00
ATOM 2 C CA . MET A 1 21.500 10.500 5.500 1.00 24.50
"#;

        let mut parser = MmcifParser::new();
        parser.parse_reader(Cursor::new(mmcif_data)).unwrap();
        
        let structure = mmcif_to_pdb_structure(&parser).unwrap();
        
        assert_eq!(structure.atoms.len(), 2);
        
        let atom1 = &structure.atoms[0];
        assert_eq!(atom1.serial, 1);
        assert_eq!(atom1.name, "N");
        assert_eq!(atom1.residue_name, "MET");
        assert_eq!(atom1.chain_id, "A");
        assert_eq!(atom1.residue_seq, 1);
        assert_eq!(atom1.x, 20.154);
        assert_eq!(atom1.y, 10.000);
        assert_eq!(atom1.z, 5.000);
    }

    #[test]
    fn test_mmcif_sequence_parsing() {
        let mmcif_data = r#"
data_test
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 MET
1 2 ALA
1 3 GLY
"#;

        let mut parser = MmcifParser::new();
        parser.parse_reader(Cursor::new(mmcif_data)).unwrap();
        
        let structure = mmcif_to_pdb_structure(&parser).unwrap();
        
        assert_eq!(structure.seqres.len(), 1);
        
        let seqres = &structure.seqres[0];
        assert_eq!(seqres.residues, vec!["MET", "ALA", "GLY"]);
        assert_eq!(seqres.num_residues, 3);
    }
}