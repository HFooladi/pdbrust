"""File functions accept path objects and report missing files clearly."""

import gzip
import shutil
from pathlib import Path

import pytest

import pdbrust

EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "pdb_files"


@pytest.fixture
def gzipped_examples(tmp_path):
    """Gzip-compressed copies of a PDB and an mmCIF example file."""
    copies = {}
    for name in ["1UBQ.pdb", "1CRN.cif"]:
        target = tmp_path / f"{name}.gz"
        with open(EXAMPLES / name, "rb") as source, gzip.open(target, "wb") as compressed:
            shutil.copyfileobj(source, compressed)
        copies[name] = target
    return copies


def test_parse_functions_accept_pathlib_paths(gzipped_examples):
    cases = [
        (pdbrust.parse_pdb_file, EXAMPLES / "1UBQ.pdb"),
        (pdbrust.parse_mmcif_file, EXAMPLES / "1CRN.cif"),
        (pdbrust.parse_structure_file, EXAMPLES / "1CRN.cif"),
        (pdbrust.parse_gzip_pdb_file, gzipped_examples["1UBQ.pdb"]),
        (pdbrust.parse_gzip_mmcif_file, gzipped_examples["1CRN.cif"]),
        (pdbrust.parse_gzip_structure_file, gzipped_examples["1UBQ.pdb"]),
    ]

    for parse, path in cases:
        assert parse(path).num_atoms > 0, parse.__name__


def test_write_functions_accept_pathlib_paths(tmp_path):
    structure = pdbrust.parse_pdb_file(str(EXAMPLES / "1UBQ.pdb"))

    pdbrust.write_pdb_file(structure, tmp_path / "written.pdb")
    pdbrust.write_mmcif_file(structure, tmp_path / "written.cif")
    pdbrust.write_gzip_mmcif_file(structure, tmp_path / "written.cif.gz")
    structure.to_file(tmp_path / "method.pdb")

    assert pdbrust.parse_pdb_file(tmp_path / "written.pdb").num_atoms == structure.num_atoms
    assert pdbrust.parse_mmcif_file(tmp_path / "written.cif").num_atoms == structure.num_atoms
    assert (
        pdbrust.parse_gzip_mmcif_file(tmp_path / "written.cif.gz").num_atoms
        == structure.num_atoms
    )
    assert pdbrust.parse_pdb_file(tmp_path / "method.pdb").num_atoms == structure.num_atoms


@pytest.mark.parametrize(
    "parse",
    [
        pdbrust.parse_pdb_file,
        pdbrust.parse_mmcif_file,
        pdbrust.parse_structure_file,
        pdbrust.parse_gzip_pdb_file,
        pdbrust.parse_gzip_mmcif_file,
        pdbrust.parse_gzip_structure_file,
    ],
    ids=lambda parse: parse.__name__,
)
def test_missing_files_raise_file_not_found_error(parse, tmp_path):
    with pytest.raises(FileNotFoundError):
        parse(str(tmp_path / "does-not-exist.pdb"))


def test_file_not_found_error_is_still_caught_as_ioerror(tmp_path):
    # Existing code that catches IOError (an alias of OSError) keeps working.
    with pytest.raises(IOError):
        pdbrust.parse_pdb_file(str(tmp_path / "does-not-exist.pdb"))
