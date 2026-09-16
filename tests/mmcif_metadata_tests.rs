//! Single-value mmCIF items (title, resolution, ...) must be read from real
//! files, where they appear after the first `loop_` and often on the line
//! following the tag.

use pdbrust::{parse_mmcif_file, parse_mmcif_string};

/// A minimal mmCIF document: one atom, then a single-value category.
fn mmcif_with_trailing_items(items: &str) -> String {
    format!(
        "data_test\n\
         loop_\n\
         _atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
         _atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n\
         _atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n\
         _atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
         _atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n\
         _atom_site.pdbx_PDB_model_num\n\
         ATOM 1 C CA . GLY A 1 1 ? 1.000 2.000 3.000 1.00 10.00 1\n\
         #\n{items}"
    )
}

#[test]
fn titles_of_real_files_are_read() {
    let cases = [
        (
            "1CRN.cif",
            "WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION. \
             PENTAGON RINGS OF WATER MOLECULES IN CRYSTALS OF CRAMBIN",
        ),
        (
            "2GB1.cif",
            "A NOVEL, HIGHLY STABLE FOLD OF THE IMMUNOGLOBULIN BINDING DOMAIN OF \
             STREPTOCOCCAL PROTEIN G",
        ),
        (
            "5HH6.cif",
            "Crystal structure of B3 metallo-beta-lactamase L1 in complex with a \
             phosphonate-based inhibitor",
        ),
    ];

    for (name, expected) in cases {
        let structure = parse_mmcif_file(format!("examples/pdb_files/{name}")).unwrap();
        assert_eq!(structure.title.as_deref(), Some(expected), "{name}");
    }
}

#[test]
fn items_following_a_loop_are_read() {
    let structure = parse_mmcif_string(&mmcif_with_trailing_items(
        "_struct.title   'A SHORT TITLE'\n",
    ))
    .unwrap();

    assert_eq!(structure.title.as_deref(), Some("A SHORT TITLE"));
    assert_eq!(structure.atoms.len(), 1);
}

#[test]
fn values_on_the_line_after_the_tag_are_read() {
    let structure = parse_mmcif_string(&mmcif_with_trailing_items(
        "_struct.title\n'A TITLE ON THE NEXT LINE'\n",
    ))
    .unwrap();

    assert_eq!(structure.title.as_deref(), Some("A TITLE ON THE NEXT LINE"));
}

#[test]
fn semicolon_delimited_text_fields_are_read() {
    let structure = parse_mmcif_string(&mmcif_with_trailing_items(
        "_struct.title\n;A TITLE IN A TEXT FIELD\n;\n",
    ))
    .unwrap();

    assert_eq!(structure.title.as_deref(), Some("A TITLE IN A TEXT FIELD"));
}

#[test]
fn unknown_and_inapplicable_values_are_not_read_as_text() {
    for value in ["?", "."] {
        let structure = parse_mmcif_string(&mmcif_with_trailing_items(&format!(
            "_struct.title  {value}\n"
        )))
        .unwrap();

        assert_eq!(structure.title, None, "value {value}");
    }
}

#[cfg(feature = "quality")]
#[test]
fn resolution_of_real_files_is_read() {
    // 1CRN and 5HH6 state it in _refine.ls_d_res_high;
    // 2GB1 is an NMR structure without a resolution.
    for (name, expected) in [
        ("1CRN.cif", Some(1.5)),
        ("5HH6.cif", Some(1.8)),
        ("2GB1.cif", None),
    ] {
        let structure = parse_mmcif_file(format!("examples/pdb_files/{name}")).unwrap();
        assert_eq!(structure.get_resolution(), expected, "{name}");
    }
}

#[test]
fn sequences_follow_the_entity_order_of_the_file() {
    let mut text = String::from(
        "data_test\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n",
    );
    for entity in 1..=6 {
        text.push_str(&format!("{entity} 1 GLY\n{entity} 2 ALA\n"));
    }
    text.push_str("#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\n");
    for (entity, chain) in ["A", "B", "C", "D", "E", "F"].iter().enumerate() {
        text.push_str(&format!("{chain} {}\n", entity + 1));
    }

    let structure = parse_mmcif_string(&text).unwrap();

    let chains: Vec<&str> = structure
        .seqres
        .iter()
        .map(|record| record.chain_id.as_str())
        .collect();
    assert_eq!(chains, vec!["A", "B", "C", "D", "E", "F"]);
}

#[cfg(feature = "quality")]
#[test]
fn cryo_em_resolution_is_read() {
    // Cryo-EM entries state the resolution in _em_3d_reconstruction, not _refine.
    let structure = parse_mmcif_string(&mmcif_with_trailing_items(
        "_exptl.method   'ELECTRON MICROSCOPY'\n#\n_em_3d_reconstruction.resolution   3.2\n",
    ))
    .unwrap();

    assert_eq!(structure.get_resolution(), Some(3.2));
}
