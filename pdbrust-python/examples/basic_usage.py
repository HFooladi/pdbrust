#!/usr/bin/env python3
"""
Basic usage example for PDBRust Python bindings.

This example demonstrates:
- Parsing PDB and mmCIF files
- Accessing structure properties
- Working with atoms, chains, and residues
- Basic filtering operations
"""

import pdbrust

# For this example, we'll use a sample file from the parent examples directory
# You can replace this with any PDB file path
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Basic Usage Example")
    print("=" * 60)

    # --- Parsing ---
    print("\n1. PARSING FILES")
    print("-" * 40)

    # Parse a PDB file
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded: {PDB_FILE}")
    print(f"  Atoms: {structure.num_atoms}")
    print(f"  Chains: {structure.num_chains}")
    print(f"  Residues: {structure.num_residues}")

    # You can also parse from strings
    # structure = pdbrust.parse_pdb_string(pdb_content)

    # Or parse mmCIF files
    # structure = pdbrust.parse_mmcif_file("protein.cif")

    # Auto-detect format
    # structure = pdbrust.parse_structure_file("protein.ent")

    # --- Structure Properties ---
    print("\n2. STRUCTURE PROPERTIES")
    print("-" * 40)

    # Get chain IDs
    chains = structure.get_chain_ids()
    print(f"Chain IDs: {chains}")

    # Get header and title (if present)
    if structure.header:
        print(f"Header: {structure.header}")
    if structure.title:
        print(f"Title: {structure.title}")

    # --- Working with Atoms ---
    print("\n3. WORKING WITH ATOMS")
    print("-" * 40)

    # Access all atoms
    atoms = structure.atoms
    print(f"Total atoms: {len(atoms)}")

    # Print first 5 atoms
    print("\nFirst 5 atoms:")
    for atom in atoms[:5]:
        print(f"  {atom.serial:4d} {atom.name:4s} {atom.residue_name:3s} "
              f"{atom.chain_id}{atom.residue_seq:4d} "
              f"({atom.x:7.3f}, {atom.y:7.3f}, {atom.z:7.3f})")

    # Atom properties
    atom = atoms[0]
    print(f"\nAtom properties for atom {atom.serial}:")
    print(f"  Name: {atom.name}")
    print(f"  Element: {atom.element}")
    print(f"  Residue: {atom.residue_name} {atom.residue_seq}")
    print(f"  Chain: {atom.chain_id}")
    print(f"  Coordinates: ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f})")
    print(f"  Occupancy: {atom.occupancy}")
    print(f"  B-factor: {atom.temp_factor}")
    print(f"  Is backbone: {atom.is_backbone()}")
    print(f"  Is hydrogen: {atom.is_hydrogen()}")

    # Calculate distance between atoms
    if len(atoms) > 1:
        dist = atoms[0].distance_to(atoms[1])
        print(f"\nDistance between atoms 1 and 2: {dist:.3f} Å")

    # --- Working with Residues ---
    print("\n4. WORKING WITH RESIDUES")
    print("-" * 40)

    # Get all residues
    residues = structure.get_residues()
    print(f"Total residues: {len(residues)}")

    # Get residues for a specific chain
    if chains:
        chain_residues = structure.get_residues_for_chain(chains[0])
        print(f"Residues in chain {chains[0]}: {len(chain_residues)}")
        print(f"First 5: {chain_residues[:5]}")

    # Get sequence from SEQRES records (if present)
    if chains:
        sequence = structure.get_sequence(chains[0])
        if sequence:
            print(f"\nSequence for chain {chains[0]}:")
            print(f"  {sequence[:50]}..." if len(sequence) > 50 else f"  {sequence}")

    # --- Filtering ---
    print("\n5. FILTERING OPERATIONS")
    print("-" * 40)

    # Method chaining for filtering
    print(f"Original structure: {structure.num_atoms} atoms")

    # Keep only CA atoms
    ca_only = structure.keep_only_ca()
    print(f"After keep_only_ca(): {ca_only.num_atoms} atoms")

    # Reload and try other filters
    structure = pdbrust.parse_pdb_file(PDB_FILE)

    # Keep only backbone atoms
    backbone = structure.keep_only_backbone()
    print(f"After keep_only_backbone(): {backbone.num_atoms} atoms")

    # Reload and chain multiple operations
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    if len(chains) > 0:
        filtered = structure.keep_only_chain(chains[0]).remove_hydrogens()
        print(f"After keep_only_chain('{chains[0]}').remove_hydrogens(): {filtered.num_atoms} atoms")

    # --- CA Coordinates ---
    print("\n6. CA COORDINATES")
    print("-" * 40)

    structure = pdbrust.parse_pdb_file(PDB_FILE)
    ca_coords = structure.get_ca_coords()
    print(f"Number of CA atoms: {len(ca_coords)}")
    if ca_coords:
        print(f"First CA coordinate: {ca_coords[0]}")

    # --- Center of Mass ---
    print("\n7. CENTER OF MASS")
    print("-" * 40)

    centroid = structure.get_centroid()
    print(f"Structure centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")

    ca_centroid = structure.get_ca_centroid()
    print(f"CA centroid: ({ca_centroid[0]:.3f}, {ca_centroid[1]:.3f}, {ca_centroid[2]:.3f})")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
