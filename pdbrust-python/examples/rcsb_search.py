#!/usr/bin/env python3
"""
RCSB PDB search and download example for PDBRust Python bindings.

This example demonstrates:
- Downloading structures from RCSB PDB
- Building search queries with various filters
- Searching for structures by different criteria
- Working with search results

NOTE: This example requires network access to RCSB PDB.
"""

import pdbrust
from pdbrust import (
    SearchQuery,
    rcsb_search,
    download_structure,
    FileFormat,
    ExperimentalMethod,
    PolymerType,
)
import tempfile
import os


def main():
    print("=" * 60)
    print("PDBRust RCSB Search & Download Example")
    print("=" * 60)

    # --- Download a Structure ---
    print("\n1. DOWNLOADING STRUCTURES")
    print("-" * 40)

    # Download ubiquitin (1UBQ) in PDB format
    print("Downloading 1UBQ...")
    try:
        structure = download_structure("1UBQ", FileFormat.pdb())
        print(f"Downloaded: {structure.num_atoms} atoms, {structure.num_chains} chains")
    except Exception as e:
        print(f"Download failed (network issue?): {e}")
        print("Continuing with other examples...")
        structure = None

    # Download in mmCIF format
    print("\nDownloading 4HHB (hemoglobin) in mmCIF format...")
    try:
        structure_cif = download_structure("4HHB", FileFormat.cif())
        print(f"Downloaded: {structure_cif.num_atoms} atoms, {structure_cif.num_chains} chains")
    except Exception as e:
        print(f"Download failed: {e}")

    # --- Download to File ---
    print("\n2. DOWNLOAD TO FILE")
    print("-" * 40)

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "1ubq.pdb")
        print(f"Downloading to {output_path}...")
        try:
            pdbrust.download_to_file("1UBQ", output_path, FileFormat.pdb())
            print(f"Saved! File size: {os.path.getsize(output_path):,} bytes")
        except Exception as e:
            print(f"Download failed: {e}")

    # --- Download as String ---
    print("\n3. DOWNLOAD AS STRING")
    print("-" * 40)

    try:
        pdb_string = pdbrust.download_pdb_string("1UBQ", FileFormat.pdb())
        lines = pdb_string.split('\n')
        print(f"Downloaded {len(lines)} lines")
        print("First 3 lines:")
        for line in lines[:3]:
            print(f"  {line}")
    except Exception as e:
        print(f"Download failed: {e}")

    # --- Basic Search ---
    print("\n4. BASIC SEARCH")
    print("-" * 40)

    # Search for kinases
    print("Searching for 'kinase'...")
    try:
        query = SearchQuery().with_text("kinase")
        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} structures")
        print("First 5 PDB IDs:")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- Advanced Search with Multiple Filters ---
    print("\n5. ADVANCED SEARCH")
    print("-" * 40)

    # Human proteins with high resolution
    print("Searching for human proteins with resolution < 2.0 Å...")
    try:
        query = (SearchQuery()
            .with_organism("Homo sapiens")
            .with_resolution_max(2.0)
            .with_experimental_method(ExperimentalMethod.xray()))

        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} structures")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- Search by Sequence Length ---
    print("\n6. SEARCH BY SEQUENCE LENGTH")
    print("-" * 40)

    # Small proteins (50-100 residues)
    print("Searching for small proteins (50-100 residues)...")
    try:
        query = (SearchQuery()
            .with_sequence_length_min(50)
            .with_sequence_length_max(100)
            .with_resolution_max(2.5))

        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} structures")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- Search by Experimental Method ---
    print("\n7. SEARCH BY EXPERIMENTAL METHOD")
    print("-" * 40)

    # NMR structures
    print("Searching for NMR structures...")
    try:
        query = (SearchQuery()
            .with_experimental_method(ExperimentalMethod.nmr())
            .with_text("ubiquitin"))

        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} NMR ubiquitin structures")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # Cryo-EM structures
    print("\nSearching for Cryo-EM structures...")
    try:
        query = (SearchQuery()
            .with_experimental_method(ExperimentalMethod.em())
            .with_text("ribosome"))

        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} Cryo-EM ribosome structures")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- Search by Polymer Type ---
    print("\n8. SEARCH BY POLYMER TYPE")
    print("-" * 40)

    # DNA structures
    print("Searching for DNA structures...")
    try:
        query = (SearchQuery()
            .with_polymer_type(PolymerType.dna())
            .with_resolution_max(3.0))

        results = rcsb_search(query, 5)
        print(f"Found {results.total_count} DNA structures")
        for pdb_id in results.pdb_ids:
            print(f"  {pdb_id}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- Combined Search and Download ---
    print("\n9. SEARCH AND DOWNLOAD WORKFLOW")
    print("-" * 40)

    print("Finding high-resolution insulin structures...")
    try:
        query = (SearchQuery()
            .with_text("insulin")
            .with_resolution_max(1.5)
            .with_experimental_method(ExperimentalMethod.xray()))

        results = rcsb_search(query, 3)
        print(f"Found {results.total_count} structures, downloading first 3...")

        for pdb_id in results.pdb_ids:
            try:
                struct = download_structure(pdb_id, FileFormat.pdb())
                rg = struct.radius_of_gyration()
                print(f"  {pdb_id}: {struct.num_atoms} atoms, Rg = {rg:.2f} Å")
            except Exception as e:
                print(f"  {pdb_id}: download failed - {e}")
    except Exception as e:
        print(f"Search failed: {e}")

    # --- File Format Options ---
    print("\n10. FILE FORMAT OPTIONS")
    print("-" * 40)

    print("Available file formats:")
    print(f"  PDB format: {FileFormat.pdb()}")
    print(f"  mmCIF format: {FileFormat.cif()}")

    print("\nFile extensions:")
    print(f"  PDB: {FileFormat.extension(FileFormat.pdb())}")
    print(f"  mmCIF: {FileFormat.extension(FileFormat.cif())}")
    print(f"  PDB compressed: {FileFormat.compressed_extension(FileFormat.pdb())}")
    print(f"  mmCIF compressed: {FileFormat.compressed_extension(FileFormat.cif())}")

    print("\n" + "=" * 60)
    print("Example completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
