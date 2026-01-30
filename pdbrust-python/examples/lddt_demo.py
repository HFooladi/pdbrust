"""
LDDT (Local Distance Difference Test) demonstration.

This example demonstrates how to calculate LDDT scores between
protein structures. LDDT is a superposition-free metric widely
used in AlphaFold (pLDDT) and CASP evaluations.

Run with:
    cd pdbrust-python
    maturin develop --release
    python examples/lddt_demo.py
"""

import pdbrust
from pdbrust import parse_pdb_file, AtomSelection, LddtOptions
import math
from pathlib import Path


def get_test_file(name: str) -> Path:
    """Get path to test PDB file."""
    return Path(__file__).parent.parent.parent / "examples" / "pdb_files" / name


def main():
    print("=== LDDT (Local Distance Difference Test) Demo ===\n")

    # Load the structure
    path = get_test_file("1UBQ.pdb")
    reference = parse_pdb_file(str(path))
    print(f"Loaded reference: 1UBQ.pdb ({reference.num_atoms} atoms)")

    # =========================================================================
    # Example 1: Self-comparison (should be perfect)
    # =========================================================================
    print("\n--- Example 1: Self-LDDT ---")

    result = reference.lddt_to(reference)
    print(f"LDDT Score: {result.score:.4f}")
    print(f"Distance pairs evaluated: {result.num_pairs}")
    print(f"Residues evaluated: {result.num_residues}")
    print(f"Per-threshold scores: {[f'{s:.4f}' for s in result.per_threshold_scores]}")

    # =========================================================================
    # Example 2: Translated structure (LDDT should still be 1.0)
    # =========================================================================
    print("\n--- Example 2: Translation Invariance ---")

    # Create a translated copy
    translated = parse_pdb_file(str(path))
    translated.translate(100.0, 50.0, 25.0)

    result = translated.lddt_to(reference)
    print("After translation of (100, 50, 25) Angstroms:")
    print(f"  LDDT Score: {result.score:.4f} (should be 1.0)")

    # Compare with RMSD (which requires alignment)
    rmsd = translated.rmsd_to(reference)
    print(f"  Direct RMSD: {rmsd:.2f} Angstroms (large without alignment)")

    # =========================================================================
    # Example 3: Per-residue LDDT for quality analysis
    # =========================================================================
    print("\n--- Example 3: Per-Residue LDDT ---")

    per_res = reference.per_residue_lddt_to(reference)

    print(f"Total residues analyzed: {len(per_res)}")
    print("\nFirst 5 residues:")
    for r in per_res[:5]:
        print(f"  {r.chain_id}{r.residue_seq} {r.residue_name}: "
              f"LDDT = {r.score:.3f} ({r.num_pairs} pairs)")

    # =========================================================================
    # Example 4: Custom options
    # =========================================================================
    print("\n--- Example 4: Custom Options ---")

    # Default options
    default_opts = LddtOptions()
    print(f"Default options: radius={default_opts.inclusion_radius}, "
          f"thresholds={default_opts.thresholds}")

    # Stricter thresholds
    strict_opts = LddtOptions(thresholds=[0.25, 0.5, 1.0])
    print(f"Strict options: radius={strict_opts.inclusion_radius}, "
          f"thresholds={strict_opts.thresholds}")

    # Smaller radius
    small_radius = LddtOptions(inclusion_radius=8.0)
    result_small = reference.lddt_to(reference, options=small_radius)
    print(f"\nWith smaller inclusion radius (8.0 A):")
    print(f"  LDDT Score: {result_small.score:.4f} ({result_small.num_pairs} pairs)")

    # =========================================================================
    # Example 5: Different atom selections
    # =========================================================================
    print("\n--- Example 5: Different Atom Selections ---")

    result_ca = reference.lddt_to(reference, AtomSelection.ca_only())
    result_bb = reference.lddt_to(reference, AtomSelection.backbone())

    print("Self-LDDT with different selections:")
    print(f"  CA atoms only: {result_ca.score:.4f} ({result_ca.num_pairs} pairs)")
    print(f"  Backbone atoms: {result_bb.score:.4f} ({result_bb.num_pairs} pairs)")

    # =========================================================================
    # Example 6: LDDT for structure quality assessment
    # =========================================================================
    print("\n--- Example 6: Quality Assessment Workflow ---")

    # In practice, you would compare a model to a reference
    # Here we simulate a "perturbed" model by creating a second copy
    model = parse_pdb_file(str(path))

    # Calculate LDDT
    result = model.lddt_to(reference)
    print(f"Model vs Reference:")
    print(f"  Global LDDT: {result.score:.4f}")

    # Find problematic regions
    per_res = model.per_residue_lddt_to(reference)
    low_quality = [r for r in per_res if r.score < 0.7]

    print(f"\nResidues with LDDT < 0.7: {len(low_quality)}")
    if low_quality:
        print("  Low quality regions:")
        for r in low_quality[:5]:
            print(f"    {r.chain_id}{r.residue_seq} {r.residue_name}: {r.score:.3f}")

    # =========================================================================
    # Example 7: LDDT vs RMSD comparison
    # =========================================================================
    print("\n--- Example 7: LDDT vs RMSD ---")

    print("Key differences:")
    print("  - LDDT is superposition-free (invariant to rotation/translation)")
    print("  - RMSD requires alignment for meaningful comparison")
    print("  - LDDT focuses on local distance preservation")
    print("  - LDDT is used in AlphaFold (pLDDT) and CASP")

    # Demonstrate with translated structure
    translated = parse_pdb_file(str(path))
    translated.translate(50.0, 50.0, 50.0)

    lddt = translated.lddt_to(reference)
    rmsd_direct = translated.rmsd_to(reference)
    aligned, align_result = translated.align_to(reference)

    print(f"\nFor translated structure (+50A in all directions):")
    print(f"  LDDT: {lddt.score:.4f} (superposition-free)")
    print(f"  Direct RMSD: {rmsd_direct:.2f} A (without alignment)")
    print(f"  Aligned RMSD: {align_result.rmsd:.4f} A (after Kabsch alignment)")

    print("\n=== Demo Complete ===")


if __name__ == "__main__":
    main()
