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


@pytest.mark.parametrize("name", PDB_FILES)
def test_gemmi_reads_mmcif_written_by_pdbrust_like_the_original(name, tmp_path):
    source = EXAMPLES / name
    output = tmp_path / "output.cif"

    pdbrust.write_mmcif_file(pdbrust.parse_pdb_file(str(source)), str(output))

    assert first_model_atoms(output) == first_model_atoms(source)
