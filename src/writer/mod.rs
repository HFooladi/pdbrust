use crate::PdbStructure;
use crate::error::PdbError;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

#[cfg(feature = "gzip")]
use flate2::write::GzEncoder;
#[cfg(feature = "gzip")]
use flate2::Compression;

/// Standard amino acid residue names (3-letter codes) for ATOM/HETATM distinction.
const STANDARD_AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Standard nucleotide residue names for ATOM/HETATM distinction.
const STANDARD_NUCLEOTIDES: &[&str] = &[
    "A", "C", "G", "U", // RNA
    "DA", "DC", "DG", "DT", // DNA
];

/// Check if a residue name is a standard biological residue (amino acid or nucleotide).
#[inline]
fn is_standard_residue(residue_name: &str) -> bool {
    let name = residue_name.trim();
    STANDARD_AMINO_ACIDS.contains(&name) || STANDARD_NUCLEOTIDES.contains(&name)
}

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

    // Write SSBOND records
    for ssbond in &structure.ssbonds {
        let icode1 = ssbond.icode1.unwrap_or(' ');
        let icode2 = ssbond.icode2.unwrap_or(' ');

        writeln!(
            writer,
            "SSBOND {:3} {:3} {}{:4}{} {:3} {}{:4}{} {:5} {:5} {:6.2}",
            ssbond.serial,
            ssbond.residue1_name,
            ssbond.chain1_id,
            ssbond.residue1_seq,
            icode1,
            ssbond.residue2_name,
            ssbond.chain2_id,
            ssbond.residue2_seq,
            icode2,
            ssbond.sym1,
            ssbond.sym2,
            ssbond.length
        )?;
    }

    // Write END record
    writeln!(writer, "END")?;

    Ok(())
}

/// Helper function to write an ATOM record
fn write_atom_record<W: Write>(writer: &mut W, atom: &crate::records::Atom) -> io::Result<()> {
    let alt_loc = atom.alt_loc.unwrap_or(' ');
    let ins_code = atom.ins_code.unwrap_or(' ');

    writeln!(
        writer,
        "ATOM  {:5} {:4}{}{:3} {}{:4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}      {:2}  ",
        atom.serial,
        atom.name,
        alt_loc,
        atom.residue_name,
        atom.chain_id,
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
    writeln!(writer, "_entry.id   {}", structure_id)?;
    writeln!(writer, "#")?;
    Ok(())
}

/// Helper function to write the _struct category.
fn write_struct_info<W: Write>(writer: &mut W, title: &str) -> io::Result<()> {
    // Quote the title if it contains spaces or special characters
    let quoted_title = if title.contains(' ') || title.contains('\'') {
        format!("\"{}\"", title.replace('"', "'"))
    } else {
        title.to_string()
    };
    writeln!(writer, "_struct.title   {}", quoted_title)?;
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
    // Determine ATOM or HETATM
    let group_pdb = if is_standard_residue(&atom.residue_name) {
        "ATOM"
    } else {
        "HETATM"
    };

    // Handle optional fields
    let alt_loc = atom.alt_loc.map_or(".".to_string(), |c| c.to_string());
    let ins_code = atom.ins_code.map_or("?".to_string(), |c| c.to_string());

    // Element symbol (use first character of atom name if element is empty)
    let element = if atom.element.is_empty() {
        atom.name.chars().next().unwrap_or('X').to_string()
    } else {
        atom.element.clone()
    };

    writeln!(
        writer,
        "{} {} {} {} {} {} {} {} {} {:.3} {:.3} {:.3} {:.2} {:.2} {}",
        group_pdb,
        atom.serial,
        element,
        atom.name,
        alt_loc,
        atom.residue_name,
        atom.chain_id,
        atom.residue_seq,
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
    let mut chain_to_entity: std::collections::HashMap<&str, i32> = std::collections::HashMap::new();
    let mut next_entity_id = 1;

    for seqres in &structure.seqres {
        let entity_id = *chain_to_entity
            .entry(&seqres.chain_id)
            .or_insert_with(|| {
                let id = next_entity_id;
                next_entity_id += 1;
                id
            });

        // Write residues from this SEQRES record
        // SEQRES records are numbered starting at 1 for each chain
        let base_num = (seqres.serial - 1) * 13; // Each SEQRES line has up to 13 residues
        for (i, residue) in seqres.residues.iter().enumerate() {
            writeln!(writer, "{} {} {}", entity_id, base_num + i as i32 + 1, residue)?;
        }
    }

    writeln!(writer, "#")?;
    Ok(())
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
            ssbond.chain1_id,
            ssbond.residue1_name,
            ssbond.residue1_seq,
            ssbond.chain2_id,
            ssbond.residue2_name,
            ssbond.residue2_seq,
            ssbond.length
        )?;
    }

    writeln!(writer, "#")?;
    Ok(())
}
