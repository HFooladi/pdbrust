#!/usr/bin/env python3
"""
Advanced filtering and manipulation example for PDBRust Python bindings.

This example demonstrates:
- Method chaining for filtering operations
- Chain normalization and residue reindexing
- Atom renumbering
- Structure centering and translation
- Residue and sequence extraction
- Practical preprocessing workflows
"""

import pdbrust

# Sample PDB files
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"
MULTI_CHAIN_FILE = "../../examples/pdb_files/4INS.pdb"


def main():
    print("=" * 60)
    print("PDBRust Advanced Filtering Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_chains} chains")

    # --- Method Chaining ---
    print("\n1. METHOD CHAINING")
    print("-" * 40)

    # Single operation
    ca_only = structure.keep_only_ca()
    print(f"Original: {structure.num_atoms} atoms")
    print(f"After keep_only_ca(): {ca_only.num_atoms} atoms")

    # Chain multiple operations
    print("\nChaining multiple operations:")
    print("  structure.remove_ligands().keep_only_chain('A').remove_hydrogens()")

    # Note: We reload to start fresh
    s = pdbrust.parse_pdb_file(PDB_FILE)
    cleaned = s.remove_ligands().keep_only_chain("A").remove_hydrogens()
    print(f"  Result: {cleaned.num_atoms} atoms")

    # More complex chain
    s = pdbrust.parse_pdb_file(PDB_FILE)
    backbone_clean = s.remove_ligands().remove_hydrogens().keep_only_backbone()
    print(f"\n  remove_ligands().remove_hydrogens().keep_only_backbone()")
    print(f"  Result: {backbone_clean.num_atoms} atoms")

    # --- Available Filter Methods ---
    print("\n2. AVAILABLE FILTER METHODS")
    print("-" * 40)

    print("""
Filtering methods (return new structure):
  remove_ligands()     - Remove HETATM records (ligands, waters)
  remove_hydrogens()   - Remove hydrogen atoms
  keep_only_chain(id)  - Keep only specified chain
  keep_only_ca()       - Keep only CA (alpha carbon) atoms
  keep_only_backbone() - Keep only backbone atoms (N, CA, C, O)
  select(query)        - PyMOL-style selection language

Comparison of approaches:
""")

    s = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"  Original: {s.num_atoms} atoms")
    print(f"  keep_only_ca(): {s.keep_only_ca().num_atoms} atoms")
    print(f"  keep_only_backbone(): {s.keep_only_backbone().num_atoms} atoms")
    print(f"  remove_ligands(): {s.remove_ligands().num_atoms} atoms")
    print(f"  remove_hydrogens(): {s.remove_hydrogens().num_atoms} atoms")

    # --- In-Place Modifications ---
    print("\n3. IN-PLACE MODIFICATIONS")
    print("-" * 40)

    # Load multi-chain structure
    try:
        multi = pdbrust.parse_pdb_file(MULTI_CHAIN_FILE)
        print(f"Loaded multi-chain structure: {multi.num_atoms} atoms")
        print(f"Original chains: {multi.get_chain_ids()}")
    except Exception as e:
        print(f"Could not load multi-chain file: {e}")
        multi = None

    # --- Chain Normalization ---
    print("\n4. CHAIN NORMALIZATION")
    print("-" * 40)

    if multi:
        # Normalize chain IDs to A, B, C, ...
        print(f"Original chains: {multi.get_chain_ids()}")
        multi.normalize_chain_ids()
        print(f"After normalize_chain_ids(): {multi.get_chain_ids()}")
    else:
        # Use single-chain structure
        s = pdbrust.parse_pdb_file(PDB_FILE)
        print(f"Original chains: {s.get_chain_ids()}")
        s.normalize_chain_ids()
        print(f"After normalize_chain_ids(): {s.get_chain_ids()}")

    # --- Residue Reindexing ---
    print("\n5. RESIDUE REINDEXING")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Get residues before reindexing
    residues_before = s.get_residues()[:5]
    print(f"Before reindex_residues():")
    for seq, name in residues_before:
        print(f"  {name} {seq}")

    # Reindex residues starting from 1
    s.reindex_residues()

    residues_after = s.get_residues()[:5]
    print(f"\nAfter reindex_residues():")
    for seq, name in residues_after:
        print(f"  {name} {seq}")

    # --- Atom Renumbering ---
    print("\n6. ATOM RENUMBERING")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Get first few atom serials before
    atoms_before = s.atoms[:5]
    print(f"Before renumber_atoms():")
    for a in atoms_before:
        print(f"  Atom {a.serial}: {a.name} {a.residue_name}")

    # Renumber atoms sequentially from 1
    s.renumber_atoms()

    atoms_after = s.atoms[:5]
    print(f"\nAfter renumber_atoms():")
    for a in atoms_after:
        print(f"  Atom {a.serial}: {a.name} {a.residue_name}")

    # --- Structure Centering ---
    print("\n7. STRUCTURE CENTERING")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Get centroid before centering
    centroid_before = s.get_centroid()
    print(f"Centroid before: ({centroid_before[0]:.3f}, {centroid_before[1]:.3f}, {centroid_before[2]:.3f})")

    # Center structure at origin
    s.center_structure()

    centroid_after = s.get_centroid()
    print(f"Centroid after center_structure(): ({centroid_after[0]:.6f}, {centroid_after[1]:.6f}, {centroid_after[2]:.6f})")

    # Also get CA centroid
    ca_centroid = s.get_ca_centroid()
    print(f"CA centroid: ({ca_centroid[0]:.6f}, {ca_centroid[1]:.6f}, {ca_centroid[2]:.6f})")

    # --- Translation ---
    print("\n8. STRUCTURE TRANSLATION")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Get centroid before
    before = s.get_centroid()
    print(f"Before translation: ({before[0]:.3f}, {before[1]:.3f}, {before[2]:.3f})")

    # Translate by (10, 20, 30)
    s.translate(10.0, 20.0, 30.0)

    after = s.get_centroid()
    print(f"After translate(10, 20, 30): ({after[0]:.3f}, {after[1]:.3f}, {after[2]:.3f})")

    # Verify the difference
    diff = (after[0] - before[0], after[1] - before[1], after[2] - before[2])
    print(f"Difference: ({diff[0]:.3f}, {diff[1]:.3f}, {diff[2]:.3f})")

    # --- Residue and Sequence Access ---
    print("\n9. RESIDUE AND SEQUENCE ACCESS")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Get all residues
    all_residues = s.get_residues()
    print(f"Total residues: {len(all_residues)}")
    print(f"First 5 residues: {all_residues[:5]}")

    # Get residues for specific chain
    chains = s.get_chain_ids()
    if chains:
        chain_residues = s.get_residues_for_chain(chains[0])
        print(f"\nResidues in chain {chains[0]}: {len(chain_residues)}")
        print(f"First 5: {chain_residues[:5]}")

    # Get sequence from SEQRES records
    if chains:
        sequence = s.get_sequence(chains[0])
        if sequence:
            print(f"\nSequence for chain {chains[0]} (from SEQRES):")
            seq_str = " ".join(sequence[:10])
            print(f"  {seq_str}... ({len(sequence)} residues)")
        else:
            print(f"\nNo SEQRES records for chain {chains[0]}")

    # --- Practical Preprocessing Workflow ---
    print("\n10. PRACTICAL PREPROCESSING WORKFLOW")
    print("-" * 40)

    print("Typical preprocessing for structure analysis:")

    # Step 1: Load
    s = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"\n1. Load: {s.num_atoms} atoms")

    # Step 2: Clean
    s = s.remove_ligands().remove_hydrogens()
    print(f"2. Remove ligands & hydrogens: {s.num_atoms} atoms")

    # Step 3: Select chain (if multi-chain)
    chains = s.get_chain_ids()
    if len(chains) > 1:
        s = s.keep_only_chain(chains[0])
        print(f"3. Keep chain {chains[0]}: {s.num_atoms} atoms")
    else:
        print(f"3. Single chain, no filtering needed")

    # Step 4: Normalize numbering
    s.normalize_chain_ids()
    s.reindex_residues()
    s.renumber_atoms()
    print("4. Normalize chain IDs, reindex residues, renumber atoms")

    # Step 5: Center
    s.center_structure()
    print(f"5. Center at origin")

    # Final state
    print(f"\nFinal structure:")
    print(f"  Atoms: {s.num_atoms}")
    print(f"  Chains: {s.get_chain_ids()}")
    print(f"  Residues: {len(s.get_residues())}")
    print(f"  Centroid: {s.get_centroid()}")

    # --- Workflow for ML Preprocessing ---
    print("\n11. ML PREPROCESSING WORKFLOW")
    print("-" * 40)

    print("""
Common preprocessing for machine learning:

1. Parse structure
   structure = pdbrust.parse_pdb_file(path)

2. Clean and filter
   cleaned = structure.remove_ligands().remove_hydrogens()

3. Select protein atoms only
   protein = cleaned.select("protein")

4. Center at origin
   protein.center_structure()

5. Extract coordinates
   coords = protein.get_coords_array()  # numpy array (N, 3)
   ca_coords = protein.get_ca_coords_array()  # CA only

6. Normalize (optional)
   protein.normalize_chain_ids()
   protein.reindex_residues()

7. Write processed structure
   protein.to_file("processed.pdb")
""")

    # --- Selection vs Method Chaining ---
    print("\n12. SELECTION LANGUAGE VS METHOD CHAINING")
    print("-" * 40)

    s = pdbrust.parse_pdb_file(PDB_FILE)

    # Method chaining
    result1 = s.keep_only_chain("A").keep_only_backbone()
    print(f"Method chaining: keep_only_chain('A').keep_only_backbone()")
    print(f"  Result: {result1.num_atoms} atoms")

    # Selection language (more flexible)
    result2 = s.select("chain A and backbone")
    print(f"\nSelection: select('chain A and backbone')")
    print(f"  Result: {result2.num_atoms} atoms")

    # Selection language advantages
    print("\nSelection language advantages:")
    print("  - More expressive: (chain A or chain B) and backbone")
    print("  - Numeric filters: bfactor < 30.0")
    print("  - Complex logic: protein and not hydrogen and bfactor < 40.0")

    # Method chaining advantages
    print("\nMethod chaining advantages:")
    print("  - Simpler syntax for common operations")
    print("  - In-place modifications available")
    print("  - No query string parsing overhead")

    # --- Summary ---
    print("\n13. SUMMARY")
    print("-" * 40)

    print("""
Key filtering methods:
  remove_ligands()      - Remove non-protein atoms
  remove_hydrogens()    - Remove H atoms
  keep_only_chain(id)   - Single chain
  keep_only_ca()        - CA atoms only
  keep_only_backbone()  - Backbone atoms only
  select(query)         - PyMOL-style selections

Key manipulation methods:
  normalize_chain_ids() - Rename chains to A, B, C...
  reindex_residues()    - Renumber from 1
  renumber_atoms()      - Sequential atom numbers
  center_structure()    - Center at origin
  translate(dx, dy, dz) - Move structure

Access methods:
  get_chain_ids()           - List chain IDs
  get_residues()            - All residues
  get_residues_for_chain()  - Chain-specific
  get_sequence()            - From SEQRES
  get_centroid()            - Center of mass
  get_ca_centroid()         - CA center of mass
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
