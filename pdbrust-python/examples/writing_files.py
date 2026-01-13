#!/usr/bin/env python3
"""
Writing files example for PDBRust Python bindings.

This example demonstrates:
- Writing structures to PDB format
- Writing structures to mmCIF format
- Writing compressed files
- Round-trip parsing and writing
- Getting structure data as strings
"""

import pdbrust
import tempfile
import os

# Sample PDB files from parent examples directory
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"
CIF_FILE = "../../examples/pdb_files/5HH6.cif"


def main():
    print("=" * 60)
    print("PDBRust Writing Files Example")
    print("=" * 60)

    # Create a temporary directory for output files
    with tempfile.TemporaryDirectory() as tmpdir:
        # --- Write PDB Format ---
        print("\n1. WRITING PDB FORMAT")
        print("-" * 40)

        # Parse a structure
        structure = pdbrust.parse_pdb_file(PDB_FILE)
        print(f"Loaded structure: {structure.num_atoms} atoms")

        # Write to PDB file
        output_pdb = os.path.join(tmpdir, "output.pdb")
        pdbrust.write_pdb_file(structure, output_pdb)
        print(f"Written to: {output_pdb}")

        # Verify the output
        reloaded = pdbrust.parse_pdb_file(output_pdb)
        print(f"Reloaded: {reloaded.num_atoms} atoms (matches: {reloaded.num_atoms == structure.num_atoms})")

        # --- Write mmCIF Format ---
        print("\n2. WRITING mmCIF FORMAT")
        print("-" * 40)

        # Write to mmCIF file
        output_cif = os.path.join(tmpdir, "output.cif")
        pdbrust.write_mmcif_file(structure, output_cif)
        print(f"Written to: {output_cif}")

        # Verify the output
        reloaded_cif = pdbrust.parse_mmcif_file(output_cif)
        print(f"Reloaded: {reloaded_cif.num_atoms} atoms")

        # --- Write Compressed mmCIF ---
        print("\n3. WRITING COMPRESSED mmCIF")
        print("-" * 40)

        output_gz = os.path.join(tmpdir, "output.cif.gz")
        pdbrust.write_gzip_mmcif_file(structure, output_gz)
        print(f"Written to: {output_gz}")

        # Check file size
        original_size = os.path.getsize(output_cif)
        compressed_size = os.path.getsize(output_gz)
        print(f"Original size: {original_size:,} bytes")
        print(f"Compressed size: {compressed_size:,} bytes")
        print(f"Compression ratio: {original_size / compressed_size:.1f}x")

        # Verify compressed file can be read back
        reloaded_gz = pdbrust.parse_gzip_mmcif_file(output_gz)
        print(f"Reloaded from gzip: {reloaded_gz.num_atoms} atoms")

        # --- Get as String ---
        print("\n4. GETTING STRUCTURE AS STRING")
        print("-" * 40)

        # Get mmCIF format string
        mmcif_string = pdbrust.write_mmcif_string(structure)
        lines = mmcif_string.split('\n')
        print(f"mmCIF string: {len(lines)} lines")
        print("First 5 lines:")
        for line in lines[:5]:
            print(f"  {line}")

        # --- Round-Trip with Modifications ---
        print("\n5. ROUND-TRIP WITH MODIFICATIONS")
        print("-" * 40)

        # Load, modify, save
        structure = pdbrust.parse_pdb_file(PDB_FILE)
        print(f"Original: {structure.num_atoms} atoms")

        # Filter to CA only
        ca_structure = structure.keep_only_ca()
        print(f"After keep_only_ca(): {ca_structure.num_atoms} atoms")

        # Center the structure
        ca_structure.center_structure()
        centroid = ca_structure.get_centroid()
        print(f"After centering, centroid: ({centroid[0]:.6f}, {centroid[1]:.6f}, {centroid[2]:.6f})")

        # Save modified structure
        output_modified = os.path.join(tmpdir, "modified.pdb")
        pdbrust.write_pdb_file(ca_structure, output_modified)
        print(f"Saved modified structure to: {output_modified}")

        # --- mmCIF Round-Trip (with HETATM) ---
        print("\n6. mmCIF ROUND-TRIP (WITH HETATM)")
        print("-" * 40)

        # Load mmCIF with ligands/ions
        if os.path.exists(CIF_FILE):
            structure_cif = pdbrust.parse_mmcif_file(CIF_FILE)
            print(f"Loaded mmCIF: {structure_cif.num_atoms} atoms")

            # Count HETATM
            hetatm_residues = set()
            for atom in structure_cif.atoms:
                if atom.residue_name not in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                              'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                              'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                              'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    hetatm_residues.add(atom.residue_name)
            print(f"Non-standard residues: {hetatm_residues}")

            # Write and reload
            output_round = os.path.join(tmpdir, "roundtrip.cif")
            pdbrust.write_mmcif_file(structure_cif, output_round)
            reloaded_round = pdbrust.parse_mmcif_file(output_round)
            print(f"After round-trip: {reloaded_round.num_atoms} atoms")
            print(f"Atom count preserved: {reloaded_round.num_atoms == structure_cif.num_atoms}")
        else:
            print(f"Skipping mmCIF round-trip test (file not found: {CIF_FILE})")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
