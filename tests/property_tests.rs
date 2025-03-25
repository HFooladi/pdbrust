use pdbrust::{parse_pdb_file, Atom, Conect};
use proptest::prelude::*;
use std::fs::File;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_pdb(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    file.write_all(content.as_bytes()).unwrap();
    file
}

fn generate_atom_record(serial: i32, x: f64, y: f64, z: f64) -> String {
    format!(
        "ATOM  {:>5} {:>4}{}{:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}      {:>4}{:>2}{:>2}",
        serial,
        "N",
        " ",
        "ALA",
        "A",
        serial,
        " ",
        x,
        y,
        z,
        1.00,
        13.79,
        "",
        "N",
        ""
    )
}

proptest! {
    #[test]
    fn test_valid_atom_records(
        serial in 1..1000i32,
        x in -1000.0..1000.0f64,
        y in -1000.0..1000.0f64,
        z in -1000.0..1000.0f64,
    ) {
        let content = generate_atom_record(serial, x, y, z);
        let file = create_test_pdb(&content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.atoms.len(), 1);
        let atom = &structure.atoms[0];
        assert_eq!(atom.serial, serial);
        assert!((atom.x - x).abs() < 1e-6);
        assert!((atom.y - y).abs() < 1e-6);
        assert!((atom.z - z).abs() < 1e-6);
        assert_eq!(atom.name, "N");
        assert_eq!(atom.residue_name, "ALA");
        assert_eq!(atom.chain_id, "A");
        assert_eq!(atom.residue_seq, serial);
        assert_eq!(atom.alt_loc, None);
        assert_eq!(atom.ins_code, None);
        assert_eq!(atom.occupancy, 1.00);
        assert_eq!(atom.temp_factor, 13.79);
        assert_eq!(atom.element, "N");
    }

    #[test]
    fn test_valid_conect_records(
        atom1 in 1..1000i32,
        atom2 in 1..1000i32,
        atom3 in 1..1000i32,
        atom4 in 1..1000i32,
    ) {
        let content = format!(
            "CONECT{:>5}{:>5}{:>5}{:>5}",
            atom1, atom2, atom3, atom4
        );
        let file = create_test_pdb(&content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.connects.len(), 1);
        let conect = &structure.connects[0];
        assert_eq!(conect.atom1, atom1);
        assert_eq!(conect.atom2, atom2);
        assert_eq!(conect.atom3, Some(atom3));
        assert_eq!(conect.atom4, Some(atom4));
    }

    #[test]
    fn test_multiple_atoms(
        atoms in prop::collection::vec(
            (1..1000i32, -1000.0..1000.0f64, -1000.0..1000.0f64, -1000.0..1000.0f64),
            1..10
        ),
    ) {
        let mut content = String::new();
        for (serial, x, y, z) in &atoms {
            content.push_str(&generate_atom_record(*serial, *x, *y, *z));
            content.push('\n');
        }
        let file = create_test_pdb(&content);
        let result = parse_pdb_file(file.path());
        assert!(result.is_ok());
        let structure = result.unwrap();
        assert_eq!(structure.atoms.len(), atoms.len());
        for (i, (serial, x, y, z)) in atoms.iter().enumerate() {
            let atom = &structure.atoms[i];
            assert_eq!(atom.serial, *serial);
            assert!((atom.x - x).abs() < 1e-6);
            assert!((atom.y - y).abs() < 1e-6);
            assert!((atom.z - z).abs() < 1e-6);
            assert_eq!(atom.name, "N");
            assert_eq!(atom.residue_name, "ALA");
            assert_eq!(atom.chain_id, "A");
            assert_eq!(atom.residue_seq, *serial);
            assert_eq!(atom.alt_loc, None);
            assert_eq!(atom.ins_code, None);
            assert_eq!(atom.occupancy, 1.00);
            assert_eq!(atom.temp_factor, 13.79);
            assert_eq!(atom.element, "N");
        }
    }
}
