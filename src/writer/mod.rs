use crate::PdbStructure;
use crate::error::PdbError;
use std::borrow::Cow;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[cfg(feature = "gzip")]
use flate2::Compression;
#[cfg(feature = "gzip")]
use flate2::write::GzEncoder;

/// Writes a PDB structure to a file.
pub fn write_pdb_file<P: AsRef<Path>>(structure: &PdbStructure, path: P) -> Result<(), PdbError> {
    let file = File::create(path)?;
    write_pdb(structure, file)
}

/// Writes a PDB structure to a writer.
pub fn write_pdb<W: Write>(structure: &PdbStructure, mut writer: W) -> Result<(), PdbError> {
    // Write header
    if let Some(header) = &structure.header {
        writeln!(writer, "HEADER    {}", header)?;
    }

    // Write title
    if let Some(title) = &structure.title {
        writeln!(writer, "TITLE     {}", title)?;
    }

    // Write remarks
    for remark in &structure.remarks {
        writeln!(writer, "REMARK {:3} {}", remark.number, remark.content)?;
    }

    // Write SEQRES records
    for seqres in &structure.seqres {
        writeln!(
            writer,
            "SEQRES {:3} {} {:4}  {}",
            seqres.serial,
            seqres.chain_id,
            seqres.num_residues,
            seqres.residues.join(" ")
        )?;
    }

    // Write SSBOND records (header section, before the coordinates)
    for ssbond in &structure.ssbonds {
        let icode1 = ssbond.icode1.unwrap_or(' ');
        let icode2 = ssbond.icode2.unwrap_or(' ');

        writeln!(
            writer,
            "SSBOND {:>3} {:>3} {} {:>4}{}   {:>3} {} {:>4}{}{:23}{:>6} {:>6} {:>5.2}",
            ssbond.serial,
            ssbond.residue1_name,
            ssbond.chain1_id,
            ssbond.residue1_seq,
            icode1,
            ssbond.residue2_name,
            ssbond.chain2_id,
            ssbond.residue2_seq,
            icode2,
            "",
            ssbond.sym1,
            ssbond.sym2,
            ssbond.length
        )?;
    }

    // Write MODEL/ATOM/ENDMDL records
    if !structure.models.is_empty() {
        // Write models if present
        for model in &structure.models {
            writeln!(writer, "MODEL     {:4}", model.serial)?;

            // Write atoms for this model
            for atom in &model.atoms {
                write_atom_record(&mut writer, atom)?;
            }

            writeln!(writer, "ENDMDL")?;
        }
    } else {
        // Write atoms directly if no models
        for atom in &structure.atoms {
            write_atom_record(&mut writer, atom)?;
        }
    }

    // Write CONECT records
    for conect in &structure.connects {
        let atom3_str = conect
            .atom3
            .map_or("     ".to_string(), |a| format!("{:5}", a));
        let atom4_str = conect
            .atom4
            .map_or("     ".to_string(), |a| format!("{:5}", a));

        writeln!(
            writer,
            "CONECT{:5}{:5}{}{}",
            conect.atom1, conect.atom2, atom3_str, atom4_str
        )?;
    }

    // Write END record
    writeln!(writer, "END")?;

    Ok(())
}

/// Helper function to write an ATOM or HETATM record
fn write_atom_record<W: Write>(writer: &mut W, atom: &crate::records::Atom) -> io::Result<()> {
    let alt_loc = atom.alt_loc.unwrap_or(' ');
    let ins_code = atom.ins_code.unwrap_or(' ');
    let record_type = if atom.is_hetatm { "HETATM" } else { "ATOM  " };
    let chain_id = if atom.chain_id.is_empty() {
        " "
    } else {
        atom.chain_id.as_str()
    };

    writeln!(
        writer,
        "{}{:5} {}{}{:>3} {}{:4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}  ",
        record_type,
        atom.serial,
        pdb_atom_name(&atom.name, &atom.element),
        alt_loc,
        atom.residue_name,
        chain_id,
        atom.residue_seq,
        ins_code,
        atom.x,
        atom.y,
        atom.z,
        atom.occupancy,
        atom.temp_factor,
        atom.element
    )
}

/// Pads an atom name to the 4-character PDB atom-name field (columns 13-16).
///
/// By convention, names of one-letter elements start in column 14 (`" CA "` is
/// an alpha carbon), while names of two-letter elements and 4-character names
/// start in column 13 (`"CA  "` is calcium, `"HD21"`).
fn pdb_atom_name(name: &str, element: &str) -> String {
    let starts_with_letter = name.chars().next().is_some_and(|c| c.is_ascii_alphabetic());
    if name.len() < 4 && element.len() < 2 && starts_with_letter {
        format!(" {name:<3}")
    } else {
        format!("{name:<4}")
    }
}

// ============================================================================
// mmCIF Writing Functions
// ============================================================================

/// Writes a PDB structure to an mmCIF file.
///
/// # Arguments
///
/// * `structure` - The PDB structure to write
/// * `path` - The file path to write to
///
/// # Examples
///
/// ```ignore
/// use pdbrust::{parse_pdb_file, write_mmcif_file};
///
/// let structure = parse_pdb_file("input.pdb")?;
/// write_mmcif_file(&structure, "output.cif")?;
/// ```
pub fn write_mmcif_file<P: AsRef<Path>>(structure: &PdbStructure, path: P) -> Result<(), PdbError> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    write_mmcif(structure, writer)
}

/// Writes a PDB structure to a gzip-compressed mmCIF file.
///
/// # Arguments
///
/// * `structure` - The PDB structure to write
/// * `path` - The file path to write to (typically ending in `.cif.gz`)
///
/// # Examples
///
/// ```ignore
/// use pdbrust::{parse_pdb_file, write_gzip_mmcif_file};
///
/// let structure = parse_pdb_file("input.pdb")?;
/// write_gzip_mmcif_file(&structure, "output.cif.gz")?;
/// ```
#[cfg(feature = "gzip")]
pub fn write_gzip_mmcif_file<P: AsRef<Path>>(
    structure: &PdbStructure,
    path: P,
) -> Result<(), PdbError> {
    let file = File::create(path)?;
    let encoder = GzEncoder::new(file, Compression::default());
    write_mmcif(structure, encoder)
}

/// Writes a PDB structure to a writer in mmCIF format.
///
/// # Arguments
///
/// * `structure` - The PDB structure to write
/// * `writer` - Any type implementing Write
///
/// # Examples
///
/// ```ignore
/// use pdbrust::{parse_pdb_file, write_mmcif};
/// use std::io::BufWriter;
/// use std::fs::File;
///
/// let structure = parse_pdb_file("input.pdb")?;
/// let file = File::create("output.cif")?;
/// write_mmcif(&structure, BufWriter::new(file))?;
/// ```
pub fn write_mmcif<W: Write>(structure: &PdbStructure, mut writer: W) -> Result<(), PdbError> {
    // Extract structure ID from header or use placeholder
    let structure_id = structure
        .header
        .as_ref()
        .and_then(|h| h.split_whitespace().last())
        .unwrap_or("XXXX");

    // Write data block header
    writeln!(writer, "data_{}", structure_id)?;
    writeln!(writer, "#")?;

    // Write _entry category
    write_entry_info(&mut writer, structure_id)?;

    // Write _struct category (title)
    if let Some(title) = &structure.title {
        write_struct_info(&mut writer, title)?;
    }

    // Write _atom_site loop (main content)
    write_atom_site_loop(&mut writer, structure)?;

    // Write _entity_poly_seq if SEQRES data exists
    if !structure.seqres.is_empty() {
        write_entity_poly_seq(&mut writer, structure)?;
    }

    // Write _struct_conn_type and _struct_disulfid if disulfide bonds exist
    if !structure.ssbonds.is_empty() {
        write_struct_disulfid(&mut writer, structure)?;
    }

    // Final comment
    writeln!(writer, "#")?;

    Ok(())
}

/// Writes a PDB structure to a String in mmCIF format.
///
/// # Arguments
///
/// * `structure` - The PDB structure to write
///
/// # Returns
///
/// A String containing the mmCIF formatted structure
///
/// # Examples
///
/// ```ignore
/// use pdbrust::{parse_pdb_file, write_mmcif_string};
///
/// let structure = parse_pdb_file("input.pdb")?;
/// let mmcif_content = write_mmcif_string(&structure)?;
/// println!("{}", mmcif_content);
/// ```
pub fn write_mmcif_string(structure: &PdbStructure) -> Result<String, PdbError> {
    let mut buffer = Vec::new();
    write_mmcif(structure, &mut buffer)?;
    Ok(String::from_utf8_lossy(&buffer).into_owned())
}

/// Helper function to write the _entry category.
fn write_entry_info<W: Write>(writer: &mut W, structure_id: &str) -> io::Result<()> {
    writeln!(writer, "_entry.id   {}", cif_value(structure_id))?;
    writeln!(writer, "#")?;
    Ok(())
}

/// Helper function to write the _struct category.
fn write_struct_info<W: Write>(writer: &mut W, title: &str) -> io::Result<()> {
    writeln!(writer, "_struct.title   {}", cif_value(title))?;
    writeln!(writer, "#")?;
    Ok(())
}

/// Helper function to write the _atom_site loop.
fn write_atom_site_loop<W: Write>(writer: &mut W, structure: &PdbStructure) -> io::Result<()> {
    writeln!(writer, "loop_")?;
    writeln!(writer, "_atom_site.group_PDB")?;
    writeln!(writer, "_atom_site.id")?;
    writeln!(writer, "_atom_site.type_symbol")?;
    writeln!(writer, "_atom_site.label_atom_id")?;
    writeln!(writer, "_atom_site.label_alt_id")?;
    writeln!(writer, "_atom_site.label_comp_id")?;
    writeln!(writer, "_atom_site.label_asym_id")?;
    writeln!(writer, "_atom_site.label_seq_id")?;
    writeln!(writer, "_atom_site.auth_seq_id")?;
    writeln!(writer, "_atom_site.pdbx_PDB_ins_code")?;
    writeln!(writer, "_atom_site.Cartn_x")?;
    writeln!(writer, "_atom_site.Cartn_y")?;
    writeln!(writer, "_atom_site.Cartn_z")?;
    writeln!(writer, "_atom_site.occupancy")?;
    writeln!(writer, "_atom_site.B_iso_or_equiv")?;
    writeln!(writer, "_atom_site.pdbx_PDB_model_num")?;

    // Determine which atoms to write
    if !structure.models.is_empty() {
        // Write atoms from all models
        for model in &structure.models {
            for atom in &model.atoms {
                write_mmcif_atom_record(writer, atom, model.serial)?;
            }
        }
    } else {
        // Write atoms directly (single model, model_num = 1)
        for atom in &structure.atoms {
            write_mmcif_atom_record(writer, atom, 1)?;
        }
    }

    writeln!(writer, "#")?;
    Ok(())
}

/// Helper function to write a single atom record in mmCIF format.
fn write_mmcif_atom_record<W: Write>(
    writer: &mut W,
    atom: &crate::records::Atom,
    model_num: i32,
) -> io::Result<()> {
    // Use is_hetatm field to determine record type
    // Per mmCIF convention: HETATM records have "." for label_seq_id
    let (group_pdb, label_seq_id) = if atom.is_hetatm {
        ("HETATM", ".".to_string())
    } else {
        ("ATOM", atom.residue_seq.to_string())
    };

    // auth_seq_id always has the numeric residue sequence
    let auth_seq_id = atom.residue_seq;

    // Handle optional fields
    let alt_loc = atom
        .alt_loc
        .map_or(".".to_string(), |c| cif_value(&c.to_string()).into_owned());
    let ins_code = atom
        .ins_code
        .map_or("?".to_string(), |c| cif_value(&c.to_string()).into_owned());

    // Element symbol (use first character of atom name if element is empty)
    let element = if atom.element.is_empty() {
        atom.name.chars().next().unwrap_or('X').to_string()
    } else {
        atom.element.clone()
    };

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {} {} {:.3} {:.3} {:.3} {:.2} {:.2} {}",
        group_pdb,
        atom.serial,
        cif_value(&element),
        cif_value(&atom.name),
        alt_loc,
        cif_value(&atom.residue_name),
        cif_value(&atom.chain_id),
        label_seq_id,
        auth_seq_id,
        ins_code,
        atom.x,
        atom.y,
        atom.z,
        atom.occupancy,
        atom.temp_factor,
        model_num
    )
}

/// Helper function to write the _entity_poly_seq loop from SEQRES records.
fn write_entity_poly_seq<W: Write>(writer: &mut W, structure: &PdbStructure) -> io::Result<()> {
    writeln!(writer, "loop_")?;
    writeln!(writer, "_entity_poly_seq.entity_id")?;
    writeln!(writer, "_entity_poly_seq.num")?;
    writeln!(writer, "_entity_poly_seq.mon_id")?;

    // Group SEQRES by chain and assign entity IDs
    let mut chain_to_entity: std::collections::HashMap<&str, i32> =
        std::collections::HashMap::new();
    let mut next_entity_id = 1;

    for seqres in &structure.seqres {
        let entity_id = *chain_to_entity.entry(&seqres.chain_id).or_insert_with(|| {
            let id = next_entity_id;
            next_entity_id += 1;
            id
        });

        // Write residues from this SEQRES record
        // SEQRES records are numbered starting at 1 for each chain
        let base_num = (seqres.serial - 1) * 13; // Each SEQRES line has up to 13 residues
        for (i, residue) in seqres.residues.iter().enumerate() {
            writeln!(
                writer,
                "{} {} {}",
                entity_id,
                base_num + i as i32 + 1,
                cif_value(residue)
            )?;
        }
    }

    writeln!(writer, "#")?;
    Ok(())
}

/// Formats a string as a CIF value, quoting it when a CIF reader would
/// otherwise misparse it (empty, containing whitespace, starting with a
/// reserved character, `.`/`?`, or a reserved word).
fn cif_value(value: &str) -> Cow<'_, str> {
    let is_reserved_word = value.get(..5).is_some_and(|prefix| {
        prefix.eq_ignore_ascii_case("data_") || prefix.eq_ignore_ascii_case("save_")
    }) || ["loop_", "stop_", "global_"]
        .iter()
        .any(|word| value.eq_ignore_ascii_case(word));
    let needs_quotes = value.is_empty()
        || value.chars().any(char::is_whitespace)
        || value.starts_with(['_', '#', '$', '\'', '"', '[', ']', ';'])
        || value == "."
        || value == "?"
        || is_reserved_word;

    if !needs_quotes {
        Cow::Borrowed(value)
    } else if !value.contains('"') {
        Cow::Owned(format!("\"{value}\""))
    } else if !value.contains('\'') {
        Cow::Owned(format!("'{value}'"))
    } else {
        Cow::Owned(format!("\"{}\"", value.replace('"', "'")))
    }
}

/// Helper function to write the _struct_conn_type and _struct_disulfid loops.
fn write_struct_disulfid<W: Write>(writer: &mut W, structure: &PdbStructure) -> io::Result<()> {
    // Write connection type
    writeln!(writer, "loop_")?;
    writeln!(writer, "_struct_conn_type.id")?;
    writeln!(writer, "_struct_conn_type.criteria")?;
    writeln!(writer, "disulf ?")?;
    writeln!(writer, "#")?;

    // Write disulfide bonds
    writeln!(writer, "loop_")?;
    writeln!(writer, "_struct_conn.id")?;
    writeln!(writer, "_struct_conn.conn_type_id")?;
    writeln!(writer, "_struct_conn.ptnr1_label_asym_id")?;
    writeln!(writer, "_struct_conn.ptnr1_label_comp_id")?;
    writeln!(writer, "_struct_conn.ptnr1_label_seq_id")?;
    writeln!(writer, "_struct_conn.ptnr2_label_asym_id")?;
    writeln!(writer, "_struct_conn.ptnr2_label_comp_id")?;
    writeln!(writer, "_struct_conn.ptnr2_label_seq_id")?;
    writeln!(writer, "_struct_conn.pdbx_dist_value")?;

    for ssbond in &structure.ssbonds {
        writeln!(
            writer,
            "disulf{} disulf {} {} {} {} {} {} {:.3}",
            ssbond.serial,
            cif_value(&ssbond.chain1_id),
            cif_value(&ssbond.residue1_name),
            ssbond.residue1_seq,
            cif_value(&ssbond.chain2_id),
            cif_value(&ssbond.residue2_name),
            ssbond.residue2_seq,
            ssbond.length
        )?;
    }

    writeln!(writer, "#")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::cif_value;

    #[test]
    fn cif_value_leaves_plain_values_bare() {
        for value in ["A", "CA", "HOH", "O5'", "1.5", "C1'"] {
            assert_eq!(cif_value(value), value);
        }
    }

    #[test]
    fn cif_value_quotes_values_a_cif_reader_would_misparse() {
        let cases = [
            ("", "\"\""),
            (" ", "\" \""),
            ("a b", "\"a b\""),
            ("a\tb", "\"a\tb\""),
            ("_x", "\"_x\""),
            ("#1", "\"#1\""),
            ("$x", "\"$x\""),
            ("'x", "\"'x\""),
            ("\"x", "'\"x'"),
            ("[x", "\"[x\""),
            ("]x", "\"]x\""),
            (";x", "\";x\""),
            (".", "\".\""),
            ("?", "\"?\""),
            ("data_x", "\"data_x\""),
            ("DATA_x", "\"DATA_x\""),
            ("save_x", "\"save_x\""),
            ("loop_", "\"loop_\""),
            ("stop_", "\"stop_\""),
            ("global_", "\"global_\""),
        ];
        for (value, expected) in cases {
            assert_eq!(cif_value(value), expected, "value {value:?}");
        }
    }

    #[test]
    fn cif_value_picks_a_quote_character_the_value_does_not_contain() {
        assert_eq!(cif_value("say \"hi\""), "'say \"hi\"'");
        assert_eq!(cif_value("it's here"), "\"it's here\"");
        // A value containing both quote characters cannot be quoted safely on
        // one line; double quotes are replaced so the output stays parseable.
        assert_eq!(cif_value("it's \"x\""), "\"it's 'x'\"");
    }
}
