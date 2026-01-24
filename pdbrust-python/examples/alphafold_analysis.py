#!/usr/bin/env python3
"""
AlphaFold/pLDDT confidence analysis example for PDBRust Python bindings.

This example demonstrates:
- Detecting predicted structures (AlphaFold/ESMFold)
- Interpreting B-factors as pLDDT confidence scores
- Per-residue pLDDT analysis with confidence categories
- Identifying low/high confidence regions
- pLDDT distribution analysis
"""

import pdbrust

# Sample PDB file - using 1UBQ as example (not actually AlphaFold)
# For real AlphaFold analysis, use an AlphaFold model (e.g., from AlphaFold DB)
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust AlphaFold/pLDDT Analysis Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_residues} residues")

    # --- Predicted Structure Detection ---
    print("\n1. PREDICTED STRUCTURE DETECTION")
    print("-" * 40)

    is_predicted = structure.is_predicted()
    print(f"Is predicted structure: {is_predicted}")

    if is_predicted:
        print("  Structure appears to be from AlphaFold/ESMFold")
        print("  B-factors interpreted as pLDDT confidence scores")
    else:
        print("  Structure appears to be experimental (X-ray/NMR)")
        print("  B-factors represent crystallographic temperature factors")
        print("\n  Note: For this demo, we'll still analyze as if it were pLDDT")

    # --- Mean pLDDT ---
    print("\n2. MEAN pLDDT SCORE")
    print("-" * 40)

    mean_plddt = structure.plddt_mean()
    print(f"Mean pLDDT: {mean_plddt:.1f}")

    # Interpretation
    if mean_plddt > 90:
        print("  Interpretation: Very high confidence (excellent prediction)")
    elif mean_plddt > 70:
        print("  Interpretation: High confidence (good prediction)")
    elif mean_plddt > 50:
        print("  Interpretation: Low confidence (treat with caution)")
    else:
        print("  Interpretation: Very low confidence (unreliable)")

    # --- Per-Residue pLDDT ---
    print("\n3. PER-RESIDUE pLDDT ANALYSIS")
    print("-" * 40)

    residue_plddt = structure.per_residue_plddt()
    print(f"Total residues: {len(residue_plddt)}")

    print(f"\n{'Chain':<6} {'Res#':<6} {'Name':<5} {'pLDDT':<8} {'Category':<15} {'Status'}")
    print("-" * 60)

    # Print first 10 residues
    for res in residue_plddt[:10]:
        status = "Confident" if res.is_confident() else ("Disordered" if res.is_disordered() else "Caution")
        print(f"{res.chain_id:<6} {res.residue_seq:<6} {res.residue_name:<5} "
              f"{res.plddt:<8.1f} {str(res.confidence_category):<15} {status}")

    if len(residue_plddt) > 10:
        print(f"... ({len(residue_plddt) - 10} more residues)")

    # --- Confidence Categories ---
    print("\n4. CONFIDENCE CATEGORY REFERENCE")
    print("-" * 40)

    print("""
pLDDT Confidence Categories:
  VeryHigh (>90):  High accuracy, well-modeled backbone and sidechains
  Confident (70-90): Generally good backbone prediction
  Low (50-70):      Should be treated with caution
  VeryLow (<50):    Likely disordered, should not be interpreted

Methods:
  is_reliable(): True for VeryHigh and Confident categories
  needs_caution(): True for Low and VeryLow categories
""")

    # --- pLDDT Distribution ---
    print("\n5. pLDDT DISTRIBUTION")
    print("-" * 40)

    very_high, confident, low, very_low = structure.plddt_distribution()

    print("Distribution of residues by confidence category:")
    print(f"  VeryHigh (>90):   {very_high * 100:5.1f}%")
    print(f"  Confident (70-90): {confident * 100:5.1f}%")
    print(f"  Low (50-70):       {low * 100:5.1f}%")
    print(f"  VeryLow (<50):     {very_low * 100:5.1f}%")

    reliable = very_high + confident
    caution = low + very_low
    print(f"\n  Total reliable: {reliable * 100:.1f}%")
    print(f"  Total needs caution: {caution * 100:.1f}%")

    # --- Low Confidence Regions ---
    print("\n6. LOW CONFIDENCE REGIONS")
    print("-" * 40)

    # Find residues below pLDDT 70 (standard threshold)
    low_conf = structure.low_confidence_regions(70.0)
    print(f"Residues with pLDDT < 70: {len(low_conf)}")

    if low_conf:
        print("\nLow confidence residues:")
        for res in low_conf[:10]:
            print(f"  {res.chain_id}{res.residue_seq} {res.residue_name}: pLDDT = {res.plddt:.1f}")
        if len(low_conf) > 10:
            print(f"  ... ({len(low_conf) - 10} more)")

    # Find very disordered regions (pLDDT < 50)
    very_low_conf = structure.low_confidence_regions(50.0)
    print(f"\nResidues with pLDDT < 50: {len(very_low_conf)}")

    # --- High Confidence Regions ---
    print("\n7. HIGH CONFIDENCE REGIONS")
    print("-" * 40)

    # Find residues above pLDDT 90 (very high confidence)
    high_conf = structure.high_confidence_regions(90.0)
    print(f"Residues with pLDDT >= 90: {len(high_conf)}")

    # Find well-predicted regions (pLDDT >= 70)
    reliable_regions = structure.high_confidence_regions(70.0)
    print(f"Residues with pLDDT >= 70: {len(reliable_regions)}")

    # --- Identify Disordered Segments ---
    print("\n8. IDENTIFYING DISORDERED SEGMENTS")
    print("-" * 40)

    # Find contiguous low-confidence regions
    disordered_segments = []
    current_segment = []

    for res in residue_plddt:
        if res.is_disordered():  # pLDDT < 50
            current_segment.append(res)
        else:
            if len(current_segment) >= 3:  # At least 3 residues
                disordered_segments.append(current_segment)
            current_segment = []

    if len(current_segment) >= 3:
        disordered_segments.append(current_segment)

    if disordered_segments:
        print(f"Found {len(disordered_segments)} disordered segments (pLDDT < 50, length >= 3):")
        for i, segment in enumerate(disordered_segments, 1):
            start = segment[0]
            end = segment[-1]
            avg_plddt = sum(r.plddt for r in segment) / len(segment)
            print(f"  Segment {i}: {start.chain_id}{start.residue_seq}-{end.residue_seq} "
                  f"({len(segment)} residues, avg pLDDT = {avg_plddt:.1f})")
    else:
        print("No significant disordered segments found")

    # --- Summary ---
    print("\n9. SUMMARY")
    print("-" * 40)

    print(f"""
Structure Analysis Summary:
  File: {PDB_FILE}
  Predicted structure: {is_predicted}
  Total residues: {len(residue_plddt)}
  Mean pLDDT: {mean_plddt:.1f}

Confidence Distribution:
  Very High (>90):   {very_high * 100:5.1f}%
  Confident (70-90): {confident * 100:5.1f}%
  Low (50-70):       {low * 100:5.1f}%
  Very Low (<50):    {very_low * 100:5.1f}%

Quality Assessment:
  Reliable residues (pLDDT >= 70): {len(reliable_regions)}
  Low confidence residues: {len(low_conf)}
  Disordered segments: {len(disordered_segments)}
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
