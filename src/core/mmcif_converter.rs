//! Converts parsed mmCIF data into PdbStructure format

use crate::core::PdbStructure;
use crate::core::mmcif::{Category, MmcifParser};
use crate::error::PdbError;
use crate::records::{Atom, Remark, SSBond, SeqRes};
use std::collections::{HashMap, HashSet};

/// Converts mmCIF data to PdbStructure
pub fn mmcif_to_pdb_structure(parser: &MmcifParser) -> Result<PdbStructure, PdbError> {
    let mut structure = PdbStructure::new();

    // Parse header and title information
    parse_header_info(parser, &mut structure)?;

    // Parse atoms from _atom_site category
    parse_atoms(parser, &mut structure)?;

    // Parse sequence information from _entity_poly_seq
    parse_sequences(parser, &mut structure)?;

    // Parse disulfide bonds from _struct_conn and _struct_disulfid
    parse_disulfide_bonds(parser, &mut structure)?;

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
                // Remove surrounding quotes (CIF values may use ' or ")
                let is_quoted = title.len() > 1
                    && ((title.starts_with('"') && title.ends_with('"'))
                        || (title.starts_with('\'') && title.ends_with('\'')));
                let clean_title = if is_quoted {
                    &title[1..title.len() - 1]
                } else {
                    title
                };
                // "?" (unknown) and "." (inapplicable) are not text.
                if !clean_title.is_empty() && clean_title != "?" && clean_title != "." {
                    structure.title = Some(clean_title.to_string());
                }
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
                structure.header =
                    Some(format!("{} | {}", existing_header, header_parts.join(", ")));
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

    for row in atom_site.rows.iter() {
        let atom = parse_atom_row(row, &col_indices)?;
        structure.atoms.push(atom);
    }

    Ok(())
}

/// Helper struct to store column indices for atom parsing
#[allow(dead_code)]
struct AtomColumnIndices {
    id: usize,
    type_symbol: usize,
    label_atom_id: usize,
    label_alt_id: Option<usize>,
    label_comp_id: usize,
    label_asym_id: usize,
    label_seq_id: usize,
    auth_seq_id: Option<usize>,
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
        headers.iter().position(|h| h == name).ok_or_else(|| {
            PdbError::InvalidRecord(format!(
                "Required column '{}' not found in _atom_site",
                name
            ))
        })
    };

    let find_optional = |name: &str| -> Option<usize> { headers.iter().position(|h| h == name) };

    Ok(AtomColumnIndices {
        id: find_required("id")?,
        type_symbol: find_required("type_symbol")?,
        label_atom_id: find_required("label_atom_id")?,
        label_alt_id: find_optional("label_alt_id"),
        label_comp_id: find_required("label_comp_id")?,
        label_asym_id: find_required("label_asym_id")?,
        label_seq_id: find_required("label_seq_id")?,
        auth_seq_id: find_optional("auth_seq_id"),
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
    // Check row length against required indices
    let max_required = [
        indices.id,
        indices.type_symbol,
        indices.label_atom_id,
        indices.label_comp_id,
        indices.label_asym_id,
        indices.label_seq_id,
        indices.cartn_x,
        indices.cartn_y,
        indices.cartn_z,
        indices.group_pdb,
    ]
    .into_iter()
    .max()
    .unwrap();

    if row.len() <= max_required {
        return Err(PdbError::ParseError(format!(
            "Row too short: has {} columns, needs at least {}",
            row.len(),
            max_required + 1
        )));
    }

    let serial: i32 = row[indices.id]
        .parse()
        .map_err(|_| PdbError::ParseError(format!("Invalid atom serial: {}", row[indices.id])))?;

    let name = row[indices.label_atom_id].clone();

    let alt_loc = if let Some(idx) = indices.label_alt_id {
        if idx >= row.len() {
            None
        } else {
            let alt_str = &row[idx];
            if alt_str == "." || alt_str == "?" || alt_str.is_empty() {
                None
            } else {
                alt_str.chars().next()
            }
        }
    } else {
        None
    };

    let residue_name = row[indices.label_comp_id].clone();
    let chain_id = row[indices.label_asym_id].clone();

    let residue_seq: i32 = {
        let label_seq = &row[indices.label_seq_id];
        if label_seq == "." || label_seq == "?" {
            // For non-polymer atoms (HETATM), label_seq_id is "." - use auth_seq_id
            if let Some(auth_idx) = indices.auth_seq_id {
                row[auth_idx].parse().map_err(|_| {
                    PdbError::ParseError(format!("Invalid auth_seq_id: {}", row[auth_idx]))
                })?
            } else {
                return Err(PdbError::ParseError(
                    "label_seq_id is '.' and auth_seq_id not available".to_string(),
                ));
            }
        } else {
            label_seq.parse().map_err(|_| {
                PdbError::ParseError(format!("Invalid residue sequence: {}", label_seq))
            })?
        }
    };

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

    let x: f64 = row[indices.cartn_x].parse().map_err(|_| {
        PdbError::ParseError(format!("Invalid x coordinate: {}", row[indices.cartn_x]))
    })?;

    let y: f64 = row[indices.cartn_y].parse().map_err(|_| {
        PdbError::ParseError(format!("Invalid y coordinate: {}", row[indices.cartn_y]))
    })?;

    let z: f64 = row[indices.cartn_z].parse().map_err(|_| {
        PdbError::ParseError(format!("Invalid z coordinate: {}", row[indices.cartn_z]))
    })?;

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

    // Determine if this is a HETATM record from group_PDB column
    let is_hetatm = row[indices.group_pdb] == "HETATM";

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

/// Parse sequence information from _entity_poly_seq
fn parse_sequences(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    let entity_poly_seq = match parser.get_category("entity_poly_seq") {
        Some(category) => category,
        None => return Ok(()), // No sequence information, not an error
    };

    // Group sequences by entity_id and build SEQRES records, keeping the
    // entities in the order they appear in the file.
    let mut sequences: HashMap<String, Vec<String>> = HashMap::new();
    let mut entity_order: Vec<String> = Vec::new();

    if let (Some(entity_col), Some(mon_id_col), Some(num_col)) = (
        entity_poly_seq.get_column("entity_id"),
        entity_poly_seq.get_column("mon_id"),
        entity_poly_seq.get_column("num"),
    ) {
        for ((entity_id, mon_id), num) in
            entity_col.iter().zip(mon_id_col.iter()).zip(num_col.iter())
        {
            if !sequences.contains_key(*entity_id) {
                entity_order.push(entity_id.to_string());
            }
            let seq_entry = sequences.entry(entity_id.to_string()).or_default();

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
    for entity_id in entity_order {
        let residues = sequences.remove(&entity_id).unwrap_or_default();
        let residues: Vec<String> = residues.into_iter().filter(|r| !r.is_empty()).collect();

        if residues.is_empty() {
            continue;
        }

        // Try to find the corresponding chain ID
        let chain_id = chain_mapping
            .get(&entity_id)
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
            struct_asym.get_column("id"),
        ) {
            for (entity_id, chain_id) in entity_col.iter().zip(id_col.iter()) {
                mapping.insert(entity_id.to_string(), chain_id.to_string());
            }
        }
    }

    mapping
}

/// Parse disulfide bonds.
///
/// Two categories can hold them:
/// - `_struct_conn` rows with `conn_type_id` `disulf`, which is where RCSB (and
///   our own mmCIF writer) records them;
/// - the older `_struct_disulfid` category.
///
/// A bond listed in both is kept once, and serials are numbered from 1.
fn parse_disulfide_bonds(
    parser: &MmcifParser,
    structure: &mut PdbStructure,
) -> Result<(), PdbError> {
    let mut bonds = Vec::new();

    if let Some(struct_disulfid) = parser.get_category("struct_disulfid") {
        for row in &struct_disulfid.rows {
            if let Some(bond) = parse_disulfide_row(row, &struct_disulfid.headers)? {
                bonds.push(bond);
            }
        }
    }

    if let Some(struct_conn) = parser.get_category("struct_conn") {
        let conn_type = struct_conn.headers.iter().position(|h| h == "conn_type_id");
        for row in &struct_conn.rows {
            let is_disulfide = conn_type
                .and_then(|i| row.get(i))
                .is_some_and(|t| t.eq_ignore_ascii_case("disulf"));
            if !is_disulfide {
                continue;
            }
            if let Some(bond) = parse_disulfide_row(row, &struct_conn.headers)? {
                bonds.push(bond);
            }
        }
    }

    let mut seen = HashSet::new();
    for mut bond in bonds {
        let end1 = (bond.chain1_id.clone(), bond.residue1_seq, bond.icode1);
        let end2 = (bond.chain2_id.clone(), bond.residue2_seq, bond.icode2);
        let key = if end1 <= end2 {
            (end1, end2)
        } else {
            (end2, end1)
        };
        if seen.insert(key) {
            bond.serial = structure.ssbonds.len() as i32 + 1;
            structure.ssbonds.push(bond);
        }
    }

    Ok(())
}

/// One end of a disulfide bond as read from an mmCIF row.
struct DisulfidePartner {
    chain_id: String,
    residue_name: String,
    residue_seq: i32,
    icode: Option<char>,
    symmetry: i32,
}

/// Builds an [`SSBond`] from one `_struct_disulfid` or `_struct_conn` row.
///
/// Author chain/residue identifiers are preferred and label identifiers used
/// when the author ones are absent. Returns `Ok(None)` when a partner cannot be
/// identified (missing columns, or `?`/`.` values), and an error for a residue
/// number that is present but not an integer. The serial is set by the caller.
fn parse_disulfide_row(row: &[String], headers: &[String]) -> Result<Option<SSBond>, PdbError> {
    let value = |name: &str| -> Option<&str> {
        let i = headers.iter().position(|h| h == name)?;
        let v = row.get(i)?.as_str();
        (v != "?" && v != ".").then_some(v)
    };
    let id = |partner: u8, item: &str| {
        value(&format!("ptnr{partner}_auth_{item}"))
            .or_else(|| value(&format!("ptnr{partner}_label_{item}")))
    };
    let partner = |p: u8| -> Result<Option<DisulfidePartner>, PdbError> {
        let (Some(chain), Some(name), Some(seq)) =
            (id(p, "asym_id"), id(p, "comp_id"), id(p, "seq_id"))
        else {
            return Ok(None);
        };
        let seq: i32 = seq
            .parse()
            .map_err(|_| PdbError::ParseError(format!("Invalid residue sequence: {}", seq)))?;
        let icode = value(&format!("pdbx_ptnr{p}_PDB_ins_code")).and_then(|s| s.chars().next());
        // "1_555" -> 1555, the PDB SSBOND convention.
        let sym = value(&format!("ptnr{p}_symmetry"))
            .and_then(|s| s.replace('_', "").parse().ok())
            .unwrap_or(1555);
        Ok(Some(DisulfidePartner {
            chain_id: chain.to_string(),
            residue_name: name.to_string(),
            residue_seq: seq,
            icode,
            symmetry: sym,
        }))
    };

    let (Some(p1), Some(p2)) = (partner(1)?, partner(2)?) else {
        return Ok(None);
    };
    let length = value("pdbx_dist_value")
        .or_else(|| value("distance"))
        .and_then(|s| s.parse().ok())
        .unwrap_or(2.04); // Typical S-S bond length

    Ok(Some(SSBond {
        serial: 0,
        residue1_name: p1.residue_name,
        chain1_id: p1.chain_id,
        residue1_seq: p1.residue_seq,
        icode1: p1.icode,
        residue2_name: p2.residue_name,
        chain2_id: p2.chain_id,
        residue2_seq: p2.residue_seq,
        icode2: p2.icode,
        sym1: p1.symmetry,
        sym2: p2.symmetry,
        length,
    }))
}

/// The first row's value of `category.field` if it is a number ("?" and "."
/// mean unknown or inapplicable).
fn numeric_item<'a>(parser: &'a MmcifParser, category: &str, field: &str) -> Option<&'a str> {
    let value = parser
        .get_category(category)?
        .get_row(0)?
        .get(field)
        .copied()?;
    value.parse::<f64>().is_ok().then_some(value)
}

/// Parse remarks and comments
fn parse_remarks(parser: &MmcifParser, structure: &mut PdbStructure) -> Result<(), PdbError> {
    // Parse from various comment categories
    if let Some(audit_conform) = parser.get_category("audit_conform") {
        for (i, _row) in audit_conform.rows.iter().enumerate() {
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

    // Parse resolution information: from refinement (X-ray, neutron) or,
    // for cryo-EM entries, from the 3D reconstruction.
    let resolution = numeric_item(parser, "refine", "ls_d_res_high")
        .or_else(|| numeric_item(parser, "em_3d_reconstruction", "resolution"));
    if let Some(resolution) = resolution {
        structure.remarks.push(Remark {
            number: 2,
            content: format!("RESOLUTION. {} ANGSTROMS.", resolution),
        });
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
