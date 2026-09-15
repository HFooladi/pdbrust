"""Files written by pdbrust must be read correctly by other tools.

gemmi is used as an independent reader: for every example structure, the
first model gemmi reads from pdbrust's output must match what it reads from
the original file.
"""

from pathlib import Path

import gemmi
import pytest

import pdbrust

EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "pdb_files"
PDB_FILES = [
    "1UBQ.pdb",
    "1HSG.pdb",
    "1L2Y.pdb",
    "8HM2.pdb",
    "AF-P62987-F1.pdb",
    "multi_model.pdb",
    "test.pdb",
]


def first_model_atoms(path):
    """Atoms of the first model as gemmi reads them."""
    structure = gemmi.read_structure(str(path))
    return [
        (
            chain.name,
            residue.seqid.num,
            residue.seqid.icode,
            residue.name,
            atom.name,
            atom.altloc,
            atom.element.name,
            round(atom.pos.x, 3),
            round(atom.pos.y, 3),
            round(atom.pos.z, 3),
            round(atom.occ, 2),
            round(atom.b_iso, 2),
        )
        for chain in structure[0]
        for residue in chain
        for atom in residue
    ]


@pytest.mark.parametrize("name", PDB_FILES)
def test_gemmi_reads_pdb_written_by_pdbrust_like_the_original(name, tmp_path):
    source = EXAMPLES / name
    output = tmp_path / "output.pdb"

    pdbrust.write_pdb_file(pdbrust.parse_pdb_file(str(source)), str(output))

    assert first_model_atoms(output) == first_model_atoms(source)


def write_large_water_box(path, n_atoms=100_010):
    """Writes, with gemmi, a structure with more than 99,999 atoms and 9,999
    residues, so its PDB file uses hybrid-36 serial and residue numbers."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for i in range(n_atoms):
        residue = gemmi.Residue()
        residue.name = "HOH"
        residue.seqid = gemmi.SeqId(i + 1, " ")
        residue.het_flag = "H"
        atom = gemmi.Atom()
        atom.name = "O"
        atom.element = gemmi.Element("O")
        atom.pos = gemmi.Position(i % 100, (i // 100) % 100, i // 10_000)
        atom.occ = 1.0
        atom.b_iso = 20.0
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.write_pdb(str(path))


def test_gemmi_reads_large_structures_written_by_pdbrust(tmp_path):
    source = tmp_path / "large.pdb"
    write_large_water_box(source)
    output = tmp_path / "output.pdb"

    pdbrust.write_pdb_file(pdbrust.parse_pdb_file(str(source)), str(output))

    atoms = first_model_atoms(output)
    assert len(atoms) == 100_010
    assert atoms == first_model_atoms(source)


def test_structures_that_do_not_fit_the_pdb_format_raise_value_error(tmp_path):
    # A two-character chain ID cannot be written to PDB columns.
    structure = pdbrust.parse_mmcif_string(
        "data_test\n"
        "loop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 C CA . GLY AB 1 1 ? 1.000 2.000 3.000 1.00 10.00 1\n"
    )
    output = tmp_path / "output.pdb"

    with pytest.raises(ValueError, match="chain ID"):
        pdbrust.write_pdb_file(structure, str(output))
    assert not output.exists()


@pytest.mark.parametrize("name", PDB_FILES)
def test_gemmi_reads_mmcif_written_by_pdbrust_like_the_original(name, tmp_path):
    source = EXAMPLES / name
    output = tmp_path / "output.cif"

    pdbrust.write_mmcif_file(pdbrust.parse_pdb_file(str(source)), str(output))

    assert first_model_atoms(output) == first_model_atoms(source)
