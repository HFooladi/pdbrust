//! Tests for the PDB and mmCIF writers: record layout, value quoting, and
//! round-trips through the parsers.

use pdbrust::{
    PdbError, PdbStructure, parse_mmcif_string, parse_pdb_file, parse_pdb_string,
    write_mmcif_string, write_pdb, write_pdb_file,
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
    // Serial numbers above 99,999 and residue numbers above 9,999 use hybrid-36
    // (the first line as written by gemmi).
    "HETATMA0000  O   HOH ABXG0      99.000  99.000   9.000  1.00 20.00           O  ",
    "ATOM  A0000  CA  GLY AA000       1.000   2.000   3.000  1.00  0.00           C  ",
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

#[test]
fn conect_records_with_hybrid36_serials_are_reproduced() {
    let structure = parse_pdb_string("CONECTA0000A0001\n").unwrap();

    let text = pdb_text(&structure);

    let written: Vec<&str> = records(&text, "CONECT")
        .into_iter()
        .map(str::trim_end)
        .collect();
    assert_eq!(written, vec!["CONECTA0000A0001"]);
}

/// Changes a structure so that it no longer fits the PDB format.
type MakeUnrepresentable = fn(&mut PdbStructure);

#[test]
fn values_that_do_not_fit_the_pdb_format_are_rejected_before_writing() {
    let base = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    let cases: [(&str, MakeUnrepresentable); 6] = [
        ("2-character chain ID", |s| {
            s.atoms[0].chain_id = "AB".into()
        }),
        ("5-character residue name", |s| {
            s.atoms[0].residue_name = "A1AAA".into()
        }),
        ("5-character atom name", |s| {
            s.atoms[0].name = "C1234".into()
        }),
        ("3-character element", |s| s.atoms[0].element = "XYZ".into()),
        ("serial beyond hybrid-36", |s| {
            s.atoms[0].serial = 87_440_032
        }),
        ("residue number below -999", |s| {
            s.atoms[0].residue_seq = -1_000
        }),
    ];

    for (case, modify) in cases {
        let mut structure = base.clone();
        modify(&mut structure);
        let mut buffer = Vec::new();

        let result = write_pdb(&structure, &mut buffer);

        assert!(
            matches!(result, Err(PdbError::InvalidRecord(_))),
            "{case}: {result:?}"
        );
        assert!(buffer.is_empty(), "{case}: output was written");
    }
}

#[test]
fn write_pdb_file_creates_no_file_for_structures_that_do_not_fit() {
    let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    structure.atoms[0].chain_id = "AB".into();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("out.pdb");

    assert!(write_pdb_file(&structure, &path).is_err());
    assert!(!path.exists());
}

// ============================================================================
// PDB writer: TITLE and SEQRES
// ============================================================================

#[test]
fn long_titles_are_written_as_continuation_records() {
    let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    structure.title = Some(
        "CRYSTAL STRUCTURE OF THE COMPLEX OF CYCLOPHILIN A WITH A HEXAPEPTIDE INHIBITOR"
            .to_string(),
    );

    let text = pdb_text(&structure);

    assert_eq!(
        records(&text, "TITLE"),
        vec![
            "TITLE     CRYSTAL STRUCTURE OF THE COMPLEX OF CYCLOPHILIN A WITH A HEXAPEPTIDE",
            "TITLE    2 INHIBITOR",
        ]
    );
}

#[test]
fn long_seqres_records_are_split_into_lines_of_13_residues() {
    let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    structure.seqres = vec![pdbrust::records::SeqRes {
        serial: 1,
        chain_id: "A".to_string(),
        num_residues: 15,
        residues: "MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR ILE THR LEU"
            .split(' ')
            .map(String::from)
            .collect(),
    }];

    let text = pdb_text(&structure);

    assert_eq!(
        records(&text, "SEQRES"),
        vec![
            "SEQRES   1 A   15  MET GLN ILE PHE VAL LYS THR LEU THR GLY LYS THR ILE",
            "SEQRES   2 A   15  THR LEU",
        ]
    );
}

#[test]
fn seqres_residue_names_are_right_justified() {
    let mut structure = parse_pdb_string(CANONICAL_ATOM_RECORDS[1]).unwrap();
    structure.seqres = vec![pdbrust::records::SeqRes {
        serial: 1,
        chain_id: "B".to_string(),
        num_residues: 3,
        residues: vec!["DA".to_string(), "DG".to_string(), "C".to_string()],
    }];

    let text = pdb_text(&structure);

    assert_eq!(
        records(&text, "SEQRES"),
        vec!["SEQRES   1 B    3   DA  DG   C"]
    );
}

#[test]
fn seqres_records_of_example_files_are_reproduced() {
    // test.pdb is left out: its SEQRES line is not in the standard columns.
    for name in [
        "1UBQ.pdb",
        "1HSG.pdb",
        "1L2Y.pdb",
        "8HM2.pdb",
        "AF-P62987-F1.pdb",
    ] {
        assert_records_reproduced(name, "SEQRES");
    }
}

#[test]
fn titles_of_wwpdb_files_are_wrapped_like_the_originals() {
    // wwPDB wraps TITLE text at word boundaries within columns 11-80.
    for name in ["1UBQ.pdb", "1HSG.pdb", "1L2Y.pdb", "8HM2.pdb"] {
        assert_records_reproduced(name, "TITLE");
    }
}

#[test]
fn titles_of_example_files_fit_in_80_columns_and_round_trip() {
    for name in EXAMPLE_PDB_FILES {
        let original = parse_pdb_file(format!("examples/pdb_files/{name}")).unwrap();

        let text = pdb_text(&original);

        assert!(
            records(&text, "TITLE").iter().all(|line| line.len() <= 80),
            "{name}"
        );
        let reparsed = parse_pdb_string(&text).unwrap();
        assert_eq!(reparsed.title, original.title, "{name}");
    }
}

/// Asserts that writing a parsed example file reproduces its `record` lines.
fn assert_records_reproduced(name: &str, record: &str) {
    let path = format!("examples/pdb_files/{name}");
    let original = std::fs::read_to_string(&path).unwrap();
    let text = pdb_text(&parse_pdb_file(&path).unwrap());

    let expected: Vec<&str> = records(&original, record)
        .into_iter()
        .map(str::trim_end)
        .collect();
    let written: Vec<&str> = records(&text, record)
        .into_iter()
        .map(str::trim_end)
        .collect();
    assert_eq!(written, expected, "{name}: {record}");
}

#[test]
fn seqres_from_mmcif_fits_in_80_columns_and_round_trips() {
    let original = pdbrust::parse_mmcif_file("examples/pdb_files/1CRN.cif").unwrap();

    let text = pdb_text(&original);

    assert!(records(&text, "SEQRES").iter().all(|line| line.len() <= 80));
    let reparsed = parse_pdb_string(&text).unwrap();
    let residues = |s: &PdbStructure| -> Vec<String> {
        s.seqres.iter().flat_map(|r| r.residues.clone()).collect()
    };
    assert_eq!(residues(&reparsed), residues(&original));
    assert_eq!(residues(&original).len(), 46);
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
