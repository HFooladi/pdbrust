use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use crate::PdbStructure;
use crate::error::PdbError;

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
        let atom3_str = conect.atom3.map_or("     ".to_string(), |a| format!("{:5}", a));
        let atom4_str = conect.atom4.map_or("     ".to_string(), |a| format!("{:5}", a));
        
        writeln!(
            writer,
            "CONECT{:5}{:5}{}{}",
            conect.atom1,
            conect.atom2,
            atom3_str,
            atom4_str
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