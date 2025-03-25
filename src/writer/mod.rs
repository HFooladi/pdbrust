use crate::models::{Atom, Conect, SeqRes, SSBond};
use crate::PdbStructure;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

/// Writes a PDB structure to a file.
pub fn write_pdb_file<P: AsRef<Path>>(structure: &PdbStructure, path: P) -> io::Result<()> {
    let mut file = File::create(path)?;
    write_pdb(structure, &mut file)
}

/// Writes a PDB structure to a writer.
pub fn write_pdb<W: Write>(structure: &PdbStructure, writer: &mut W) -> io::Result<()> {
    // Write header
    if !structure.header.is_empty() {
        writeln!(writer, "HEADER    {}", structure.header)?;
    }

    // Write title
    if !structure.title.is_empty() {
        writeln!(writer, "TITLE     {}", structure.title)?;
    }

    // Write remarks
    for remark in &structure.remarks {
        writeln!(writer, "REMARK   1 {}", remark)?;
    }

    // Write SEQRES records
    for seqres in &structure.seqres {
        writeln!(
            writer,
            "SEQRES {:>4} {} {:>4}  {}",
            seqres.serial,
            seqres.chain_id,
            seqres.num_residues,
            seqres.residues.join(" ")
        )?;
    }

    // Write ATOM records
    for atom in &structure.atoms {
        writeln!(
            writer,
            "ATOM  {:>5} {:>4}{}{:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}      {:>4}{:>2}{:>2}",
            atom.serial,
            atom.name,
            atom.alt_loc.unwrap_or(' '),
            atom.residue_name,
            atom.chain_id,
            atom.residue_number,
            atom.ins_code.unwrap_or(' '),
            atom.x,
            atom.y,
            atom.z,
            atom.occupancy,
            atom.b_factor,
            "",
            atom.element,
            ""
        )?;
    }

    // Write CONECT records
    for conect in &structure.conect {
        writeln!(
            writer,
            "CONECT{:>5}{:>5}{:>5}{:>5}",
            conect.atom_serial,
            conect.bonded_atoms[0],
            conect.bonded_atoms.get(1).map_or("", |&n| &format!("{:>5}", n)),
            conect.bonded_atoms.get(2).map_or("", |&n| &format!("{:>5}", n))
        )?;
    }

    // Write SSBOND records
    for ssbond in &structure.ssbond {
        writeln!(
            writer,
            "SSBOND{:>4} {:>3} {} {:>4}{} {:>3} {} {:>4}{} {:>5} {:>5} {:>6.2}",
            ssbond.serial,
            ssbond.residue1_name,
            ssbond.chain1_id,
            ssbond.residue1_seq,
            ssbond.icode1.unwrap_or(' '),
            ssbond.residue2_name,
            ssbond.chain2_id,
            ssbond.residue2_seq,
            ssbond.icode2.unwrap_or(' '),
            ssbond.sym1,
            ssbond.sym2,
            ssbond.length
        )?;
    }

    // Write END record
    writeln!(writer, "END")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_structure() -> PdbStructure {
        PdbStructure {
            header: "PROTEIN".to_string(),
            title: "TEST PROTEIN".to_string(),
            remarks: vec!["THIS IS A TEST".to_string()],
            seqres: vec![SeqRes {
                serial: 1,
                chain_id: "A".to_string(),
                num_residues: 20,
                residues: vec!["SER".to_string(); 20],
            }],
            atoms: vec![
                Atom {
                    serial: 1,
                    name: "N".to_string(),
                    alt_loc: Some(' '),
                    residue_name: "ASP".to_string(),
                    chain_id: "A".to_string(),
                    residue_number: 30,
                    ins_code: Some(' '),
                    x: 27.360,
                    y: 44.310,
                    z: 48.500,
                    occupancy: 1.00,
                    b_factor: 13.79,
                    element: "N".to_string(),
                },
                Atom {
                    serial: 2,
                    name: "CA".to_string(),
                    alt_loc: Some(' '),
                    residue_name: "ASP".to_string(),
                    chain_id: "A".to_string(),
                    residue_number: 30,
                    ins_code: Some(' '),
                    x: 26.360,
                    y: 43.310,
                    z: 48.500,
                    occupancy: 1.00,
                    b_factor: 13.79,
                    element: "C".to_string(),
                },
            ],
            conect: vec![Conect {
                atom_serial: 1,
                bonded_atoms: vec![2],
            }],
            ssbond: vec![SSBond {
                serial: 1,
                residue1_name: "CYS".to_string(),
                chain1_id: "A".to_string(),
                residue1_seq: 85,
                icode1: Some(' '),
                residue2_name: "CYS".to_string(),
                chain2_id: "A".to_string(),
                residue2_seq: 101,
                icode2: Some(' '),
                sym1: 1555,
                sym2: 1555,
                length: 2.03,
            }],
        }
    }

    #[test]
    fn test_write_pdb() {
        let structure = create_test_structure();
        let mut file = NamedTempFile::new().unwrap();
        write_pdb(&structure, &mut file).unwrap();
        file.flush().unwrap();

        let mut content = String::new();
        file.seek(std::io::SeekFrom::Start(0)).unwrap();
        std::io::Read::read_to_string(&mut file, &mut content).unwrap();

        let expected = "\
HEADER    PROTEIN
TITLE     TEST PROTEIN
REMARK   1 THIS IS A TEST
SEQRES    1 A   20  SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER SER
ATOM      1  N   ASP A  30      27.360  44.310  48.500  1.00 13.79           N
ATOM      2  CA  ASP A  30      26.360  43.310  48.500  1.00 13.79           C
CONECT    1    2
SSBOND   1 CYS A   85  CYS A  101  1555  1555   2.03
END
";
        assert_eq!(content, expected);
    }
}
