#!/usr/bin/env python3
"""
Geometry and RMSD example for PDBRust Python bindings.

This example demonstrates:
- RMSD calculation between structures
- Structure alignment using Kabsch algorithm
- Different atom selections (CA, backbone, all atoms)
- Per-residue RMSD for flexibility analysis
"""

import pdbrust
from pdbrust import AtomSelection

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Geometry & RMSD Example")
    print("=" * 60)

    # Load structure
    structure1 = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure1.num_atoms} atoms")

    # Create a modified copy for comparison
    # (In real use, you'd load two different structures)
    structure2 = pdbrust.parse_pdb_file(PDB_FILE)

    # Translate structure2 to create a difference
    structure2.translate(1.0, 0.5, 0.2)
    print("Created second structure by translating 1.0, 0.5, 0.2 Å")

    # --- Basic RMSD ---
    print("\n1. BASIC RMSD CALCULATION")
    print("-" * 40)

    # RMSD without alignment (measures structural difference)
    rmsd = structure1.rmsd_to(structure2)
    print(f"RMSD (no alignment): {rmsd:.4f} Å")

    # --- Atom Selections ---
    print("\n2. RMSD WITH DIFFERENT ATOM SELECTIONS")
    print("-" * 40)

    # CA atoms only (default)
    rmsd_ca = structure1.rmsd_to(structure2, AtomSelection.ca_only())
    print(f"RMSD (CA only): {rmsd_ca:.4f} Å")

    # Backbone atoms (N, CA, C, O)
    rmsd_bb = structure1.rmsd_to(structure2, AtomSelection.backbone())
    print(f"RMSD (backbone): {rmsd_bb:.4f} Å")

    # All atoms
    rmsd_all = structure1.rmsd_to(structure2, AtomSelection.all_atoms())
    print(f"RMSD (all atoms): {rmsd_all:.4f} Å")

    # Custom atom selection
    rmsd_custom = structure1.rmsd_to(structure2, AtomSelection.custom(["CA", "CB"]))
    print(f"RMSD (CA + CB): {rmsd_custom:.4f} Å")

    # --- Structure Alignment ---
    print("\n3. STRUCTURE ALIGNMENT (KABSCH ALGORITHM)")
    print("-" * 40)

    # Align structure1 to structure2
    # This returns:
    # - aligned: the transformed structure
    # - result: AlignmentResult with RMSD, rotation, translation
    aligned, result = structure1.align_to(structure2)

    print(f"Alignment RMSD: {result.rmsd:.4f} Å")
    print(f"Number of atoms used: {result.num_atoms}")

    # The aligned structure is now superimposed on structure2
    # Verify by calculating RMSD again (should be ~0 after alignment)
    rmsd_after = aligned.rmsd_to(structure2)
    print(f"RMSD after alignment: {rmsd_after:.6f} Å")

    # Alignment with specific atom selection
    aligned_bb, result_bb = structure1.align_to(structure2, AtomSelection.backbone())
    print(f"\nBackbone alignment RMSD: {result_bb.rmsd:.4f} Å ({result_bb.num_atoms} atoms)")

    # --- Per-Residue RMSD ---
    print("\n4. PER-RESIDUE RMSD (FLEXIBILITY ANALYSIS)")
    print("-" * 40)

    # Create a structure with local deviations for demonstration
    structure3 = pdbrust.parse_pdb_file(PDB_FILE)

    # In real analysis, you might compare:
    # - Two conformations of the same protein
    # - NMR models
    # - MD simulation snapshots
    # - Homologous structures

    # For this demo, we'll add varying translations to different regions
    # (This is artificial - real differences would come from different structures)

    # Calculate per-residue RMSD
    per_res = structure1.per_residue_rmsd_to(structure2)

    print(f"Number of residues analyzed: {len(per_res)}")
    print("\nPer-residue RMSD (first 10):")
    print(f"{'Chain':<6} {'ResSeq':<8} {'Name':<5} {'RMSD (Å)':<10} {'Atoms':<6}")
    print("-" * 40)

    for r in per_res[:10]:
        print(f"{r.chain_id:<6} {r.residue_seq:<8} {r.residue_name:<5} {r.rmsd:<10.4f} {r.num_atoms:<6}")

    # Find flexible regions (high RMSD)
    threshold = 0.5  # Å
    flexible = [r for r in per_res if r.rmsd > threshold]
    if flexible:
        print(f"\nFlexible regions (RMSD > {threshold} Å): {len(flexible)}")
        for r in flexible[:5]:
            print(f"  {r.chain_id}{r.residue_seq} {r.residue_name}: {r.rmsd:.4f} Å")
    else:
        print(f"\nNo residues with RMSD > {threshold} Å (structures are very similar)")

    # --- Practical Example: Self-Comparison ---
    print("\n5. SELF-COMPARISON (SHOULD BE ZERO)")
    print("-" * 40)

    # RMSD of structure with itself should be exactly 0
    self_rmsd = structure1.rmsd_to(structure1)
    print(f"Self-RMSD: {self_rmsd:.10f} Å")

    # Self-alignment should also give 0 RMSD
    self_aligned, self_result = structure1.align_to(structure1)
    print(f"Self-alignment RMSD: {self_result.rmsd:.10f} Å")

    # --- Summary Statistics ---
    print("\n6. SUMMARY")
    print("-" * 40)

    print(f"""
Structure Comparison Summary:
  - Structure 1: {structure1.num_atoms} atoms
  - Structure 2: {structure2.num_atoms} atoms (translated)
  - RMSD before alignment: {rmsd:.4f} Å
  - RMSD after alignment: {rmsd_after:.6f} Å
  - Atoms used for alignment: {result.num_atoms}
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
