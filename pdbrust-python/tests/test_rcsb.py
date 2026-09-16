"""RCSB search and download are available in every wheel.

Tests that contact RCSB only run with PDBRUST_NETWORK_TESTS=1.
"""

import os

import pytest

import pdbrust

RCSB_API = [
    "FileFormat",
    "ExperimentalMethod",
    "PolymerType",
    "SearchQuery",
    "SearchResult",
    "rcsb_search",
    "download_structure",
    "download_pdb_string",
    "download_to_file",
    "AsyncDownloadOptions",
    "DownloadResult",
    "download_multiple",
]

needs_network = pytest.mark.skipif(
    os.environ.get("PDBRUST_NETWORK_TESTS") != "1",
    reason="set PDBRUST_NETWORK_TESTS=1 to run tests that contact RCSB",
)


def test_rcsb_api_is_available():
    missing = [name for name in RCSB_API if not hasattr(pdbrust, name)]
    assert missing == []


@needs_network
def test_download_structure_over_https():
    structure = pdbrust.download_structure("1UBQ", pdbrust.FileFormat.pdb())
    assert structure.num_atoms > 500


@needs_network
def test_download_to_file_over_https(tmp_path):
    path = tmp_path / "1UBQ.cif"
    pdbrust.download_to_file("1UBQ", path, pdbrust.FileFormat.cif())
    assert pdbrust.parse_mmcif_file(path).num_atoms > 500


@needs_network
def test_rcsb_search_over_https():
    results = pdbrust.rcsb_search(pdbrust.SearchQuery().with_text("ubiquitin"), 5)
    assert len(results.pdb_ids) > 0


@needs_network
def test_download_multiple_over_https():
    results = pdbrust.download_multiple(["1UBQ", "1CRN"], pdbrust.FileFormat.pdb())
    assert [result.success for result in results] == [True, True]
