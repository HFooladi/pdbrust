//! Tests for the PDB and mmCIF writers: record layout, value quoting, and
//! round-trips through the parsers.

use pdbrust::{
    PdbStructure, parse_mmcif_string, parse_pdb_file, parse_pdb_string, write_mmcif_string,
    write_pdb,
};

/// Writes `structure` as PDB-format text.
fn pdb_text(structure: &PdbStructure) -> String {
    let mut buffer = Vec::new();
    write_pdb(structure, &mut buffer).unwrap();
    String::from_utf8(buffer).unwrap()
}

/// Returns the lines of `text` that start with `record`.
fn records<'a>(text: &'a str, record: &str) -> Vec<&'a str> {
    text.lines()
        .filter(|line| line.starts_with(record))
        .collect()
}

/// Asserts that two structures contain the same atoms, field by field.
fn assert_same_atoms(expected: &PdbStructure, actual: &PdbStructure, context: &str) {
    assert_eq!(
        expected.atoms.len(),
        actual.atoms.len(),
        "{context}: atom count"
    );
    for (e, a) in expected.atoms.iter().zip(&actual.atoms) {
        assert_eq!(
            (
                e.serial,
                &e.name,
                e.alt_loc,
                &e.residue_name,
                &e.chain_id,
                e.residue_seq,
                e.ins_code,
                &e.element,
                e.is_hetatm
            ),
            (
                a.serial,
                &a.name,
                a.alt_loc,
                &a.residue_name,
                &a.chain_id,
                a.residue_seq,
                a.ins_code,
                &a.element,
                a.is_hetatm
            ),
            "{context}: atom {}",
            e.serial
        );
        for (field, x, y) in [
            ("x", e.x, a.x),
            ("y", e.y, a.y),
            ("z", e.z, a.z),
            ("occupancy", e.occupancy, a.occupancy),
            ("temp_factor", e.temp_factor, a.temp_factor),
        ] {
            assert!(
                (x - y).abs() < 1e-6,
                "{context}: atom {} {field}: {x} != {y}",
                e.serial
            );
        }
    }
}

const EXAMPLE_PDB_FILES: &[&str] = &[
    "1UBQ.pdb",
    "1HSG.pdb",
    "1L2Y.pdb",
    "8HM2.pdb",
    "AF-P62987-F1.pdb",
    "multi_model.pdb",
    "test.pdb",
];

// ============================================================================
// PDB writer: ATOM/HETATM layout
// ============================================================================

/// Coordinate records laid out as in wwPDB files (PDB format v3.3). Writing a
/// structure parsed from one of these lines must reproduce it exactly.
const CANONICAL_ATOM_RECORDS: &[&str] = &[
    // Names of one-letter elements start in column 14.
    "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N  ",
    "ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C  ",
    "ATOM      2  O5'  DA B   1      13.000  23.000  33.000  1.00 18.00           O  ",
    // Four-character names, and names starting with a digit (older hydrogen
    // naming), start in column 13.
    "ATOM     15 HD21 ASN A   1     -11.572   3.791  -4.444  1.00  0.00           H  ",
    "ATOM      5 1HB  ALA A   1      10.000  20.000  30.000  1.00  0.00           H  ",
    // Names of two-letter elements start in column 13, so a calcium ion is
    // distinguishable from an alpha carbon.
    "HETATM  605 ZN    ZN A 101      10.000  20.000  30.000  1.00 15.00          ZN  ",
    "HETATM  606 CA    CA A 102      11.000  21.000  31.000  1.00 16.00          CA  ",
    // Residue names are right-justified in columns 18-20.
    "ATOM      1  P     A B   1      12.000  22.000  32.000  1.00 17.00           P  ",
    "HETATM  604  O   HOH A  77      45.747  30.081  19.708  1.00 12.43           O  ",
];

#[test]
fn writing_canonical_atom_records_reproduces_them_exactly() {
    for &line in CANONICAL_ATOM_RECORDS {
        let structure = parse_pdb_string(line).unwrap();
        let text = pdb_text(&structure);
        let record = &line[..6];
        assert_eq!(records(&text, record), vec![line]);
    }
}

#[test]
fn empty_chain_id_is_written_as_a_blank_column() {
    let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    structure.atoms[0].chain_id = String::new();

    let text = pdb_text(&structure);

    assert_eq!(
        records(&text, "ATOM"),
        vec!["ATOM      2  CA  MET     1      26.266  25.413   2.842  1.00 10.38           C  "]
    );
}

#[test]
fn pdb_round_trip_preserves_atoms_of_example_files() {
    for name in EXAMPLE_PDB_FILES {
        let original = parse_pdb_file(format!("examples/pdb_files/{name}")).unwrap();
        let reparsed = parse_pdb_string(&pdb_text(&original)).unwrap();
        assert_same_atoms(&original, &reparsed, name);
    }
}

// ============================================================================
// PDB reader and writer: SSBOND
// ============================================================================

const CANONICAL_SSBOND_RECORDS: &[&str] = &[
    "SSBOND   1 CYS A    6    CYS A   11                          1555   1555  2.05  ",
    "SSBOND   3 CYS L   23    CYS L   88A                         1555   1555  2.04  ",
];

#[test]
fn ssbond_symmetry_operators_are_read_from_columns_60_to_72() {
    let structure = parse_pdb_string(CANONICAL_SSBOND_RECORDS[0]).unwrap();

    let bond = &structure.ssbonds[0];
    assert_eq!((bond.sym1, bond.sym2), (1555, 1555));
    assert!((bond.length - 2.05).abs() < 1e-9);
}

#[test]
fn writing_canonical_ssbond_records_reproduces_them() {
    for &line in CANONICAL_SSBOND_RECORDS {
        let structure = parse_pdb_string(line).unwrap();
        let text = pdb_text(&structure);
        let written: Vec<&str> = records(&text, "SSBOND")
            .into_iter()
            .map(str::trim_end)
            .collect();
        assert_eq!(written, vec![line.trim_end()]);
    }
}

#[test]
fn ssbond_records_are_written_before_coordinate_records() {
    let input = format!(
        "{}\n{}\n",
        CANONICAL_SSBOND_RECORDS[0], CANONICAL_ATOM_RECORDS[1]
    );
    let structure = parse_pdb_string(&input).unwrap();

    let text = pdb_text(&structure);

    let position = |record: &str| text.lines().position(|line| line.starts_with(record));
    assert!(position("SSBOND") < position("ATOM"), "{text}");
}

// ============================================================================
// mmCIF writer
// ============================================================================

#[test]
fn mmcif_writer_quotes_blank_and_empty_chain_ids() {
    for chain in [" ", ""] {
        let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
        structure.atoms[0].chain_id = chain.to_string();

        let reparsed = parse_mmcif_string(&write_mmcif_string(&structure).unwrap()).unwrap();

        let atom = &reparsed.atoms[0];
        assert_eq!(atom.chain_id, chain);
        assert_eq!((atom.x, atom.y, atom.z), (26.266, 25.413, 2.842));
        assert_eq!((atom.occupancy, atom.temp_factor), (1.0, 10.38));
    }
}

#[test]
fn mmcif_round_trip_preserves_titles() {
    for title in [
        "PLAIN",
        "WITH SPACES",
        "O'NEIL'S PROTEIN",
        "THE \"X\" COMPLEX",
    ] {
        let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
        structure.title = Some(title.to_string());

        let reparsed = parse_mmcif_string(&write_mmcif_string(&structure).unwrap()).unwrap();

        assert_eq!(reparsed.title.as_deref(), Some(title));
    }
}

#[test]
fn mmcif_round_trip_preserves_atoms_of_example_files() {
    for name in EXAMPLE_PDB_FILES {
        let original = parse_pdb_file(format!("examples/pdb_files/{name}")).unwrap();
        let reparsed = parse_mmcif_string(&write_mmcif_string(&original).unwrap()).unwrap();
        assert_same_atoms(&original, &reparsed, name);
    }
}
