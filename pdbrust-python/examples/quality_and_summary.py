#!/usr/bin/env python3
"""
Quality assessment and structure summary example for PDBRust Python bindings.

This example demonstrates:
- Quality reports for structure assessment
- Analysis readiness checks
- Unified structure summaries (quality + descriptors)
- Structure descriptors computation
- Dictionary serialization
- CSV export functionality
"""

import pdbrust

# Sample PDB files
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"
MMCIF_FILE = "../../examples/pdb_files/1CRN.cif"


def main():
    print("=" * 60)
    print("PDBRust Quality & Summary Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms")

    # --- Quality Report ---
    print("\n1. QUALITY REPORT")
    print("-" * 40)

    # Get comprehensive quality assessment
    quality = structure.quality_report()

    # Structure characteristics
    print("Structure characteristics:")
    print(f"  Atoms:    {quality.num_atoms}")
    print(f"  Chains:   {quality.num_chains}")
    print(f"  Models:   {quality.num_models}")
    print(f"  Residues: {quality.num_residues}")

    # Quality flags
    print("\nQuality flags:")
    print(f"  CA-only structure: {quality.has_ca_only}")
    print(f"  Multiple models:   {quality.has_multiple_models}")
    print(f"  Alternate locs:    {quality.has_altlocs}")
    print(f"  HETATM records:    {quality.has_hetatm}")
    print(f"  Hydrogen atoms:    {quality.has_hydrogens}")
    print(f"  Disulfide bonds:   {quality.has_ssbonds}")
    print(f"  CONECT records:    {quality.has_conect}")

    # --- Analysis Readiness ---
    print("\n2. ANALYSIS READINESS CHECKS")
    print("-" * 40)

    # is_analysis_ready: single model, no altlocs, not CA-only
    is_ready = quality.is_analysis_ready()
    print(f"Analysis ready: {is_ready}")
    if not is_ready:
        reasons = []
        if quality.has_multiple_models:
            reasons.append("multiple models")
        if quality.has_altlocs:
            reasons.append("alternate locations")
        if quality.has_ca_only:
            reasons.append("CA-only")
        print(f"  Reasons: {', '.join(reasons)}")

    # is_clean: no altlocs, not CA-only
    is_clean = quality.is_clean()
    print(f"Clean structure: {is_clean}")

    # Also available directly on structure
    print(f"\nDirect method calls:")
    print(f"  has_ca_only(): {structure.has_ca_only()}")
    print(f"  has_multiple_models(): {structure.has_multiple_models()}")
    print(f"  has_altlocs(): {structure.has_altlocs()}")
    print(f"  has_hetatm(): {structure.has_hetatm()}")
    print(f"  has_hydrogens(): {structure.has_hydrogens()}")

    # --- Resolution (if available) ---
    print("\n3. RESOLUTION")
    print("-" * 40)

    resolution = structure.get_resolution()
    if resolution is not None:
        print(f"Resolution: {resolution:.2f} A")
    else:
        print("Resolution: Not available (may be NMR or missing REMARK 2)")

    # --- Structure Descriptors ---
    print("\n4. STRUCTURE DESCRIPTORS")
    print("-" * 40)

    # Get all descriptors at once
    descriptors = structure.structure_descriptors()

    print("Size metrics:")
    print(f"  Residues: {descriptors.num_residues}")
    print(f"  Atoms: {descriptors.num_atoms}")

    print("\nGeometric descriptors:")
    print(f"  Radius of gyration: {descriptors.radius_of_gyration:.2f} A")
    print(f"  Max CA distance: {descriptors.max_ca_distance:.2f} A")
    print(f"  Compactness index: {descriptors.compactness_index:.3f}")
    print(f"  CA density: {descriptors.ca_density:.6f}")

    print("\nComposition:")
    print(f"  Glycine ratio: {descriptors.glycine_ratio:.3f}")
    print(f"  Hydrophobic ratio: {descriptors.hydrophobic_ratio:.3f}")
    print(f"  Missing residue ratio: {descriptors.missing_residue_ratio:.3f}")
    print(f"  Secondary structure ratio: {descriptors.secondary_structure_ratio:.3f}")

    print("\nB-factor statistics:")
    print(f"  Mean (all atoms): {descriptors.b_factor_mean:.2f} A^2")
    print(f"  Mean (CA only): {descriptors.b_factor_mean_ca:.2f} A^2")
    print(f"  Min: {descriptors.b_factor_min:.2f} A^2")
    print(f"  Max: {descriptors.b_factor_max:.2f} A^2")
    print(f"  Std: {descriptors.b_factor_std:.2f} A^2")

    # Amino acid composition
    print("\nAmino acid composition (top 5):")
    aa_comp = descriptors.aa_composition
    sorted_aa = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)
    for aa, frac in sorted_aa[:5]:
        print(f"  {aa}: {frac * 100:.1f}%")

    # --- Unified Summary ---
    print("\n5. UNIFIED STRUCTURE SUMMARY")
    print("-" * 40)

    # Get unified summary (combines quality + descriptors)
    summary = structure.summary()

    print("Summary object repr:")
    print(f"  {repr(summary)}")

    # Summary includes both quality and descriptor fields
    print("\nSummary fields:")
    print(f"  num_atoms: {summary.num_atoms}")
    print(f"  num_residues: {summary.num_residues}")
    print(f"  num_chains: {summary.num_chains}")
    print(f"  radius_of_gyration: {summary.radius_of_gyration:.2f} A")
    print(f"  has_altlocs: {summary.has_altlocs}")
    print(f"  is_analysis_ready: {summary.is_analysis_ready()}")

    # --- Dictionary Serialization ---
    print("\n6. DICTIONARY SERIALIZATION")
    print("-" * 40)

    # Quality report to dict
    quality_dict = quality.to_dict()
    print("QualityReport.to_dict() keys:")
    print(f"  {list(quality_dict.keys())}")

    # Descriptors to dict
    desc_dict = descriptors.to_dict()
    print("\nStructureDescriptors.to_dict() keys:")
    print(f"  {list(desc_dict.keys())}")

    # Summary to dict
    summary_dict = summary.to_dict()
    print("\nStructureSummary.to_dict() keys:")
    print(f"  {list(summary_dict.keys())}")

    # --- CSV Export ---
    print("\n7. CSV EXPORT")
    print("-" * 40)

    # Get field names for CSV header
    field_names = summary.field_names()
    print(f"CSV field names ({len(field_names)}):")
    for i in range(0, len(field_names), 5):
        print(f"  {', '.join(field_names[i:i+5])}")

    # Get CSV values
    csv_values = summary.to_csv_values()
    print(f"\nCSV values (first 10):")
    for name, value in zip(field_names[:10], csv_values[:10]):
        print(f"  {name}: {value}")
    print(f"  ... ({len(csv_values) - 10} more values)")

    # --- Comparing Multiple Structures ---
    print("\n8. COMPARING STRUCTURES")
    print("-" * 40)

    # Load second structure
    try:
        structure2 = pdbrust.parse_mmcif_file(MMCIF_FILE)
        summary2 = structure2.summary()

        print(f"{'Metric':<25} {'1UBQ':<15} {'1CRN':<15}")
        print("-" * 55)
        print(f"{'Atoms':<25} {summary.num_atoms:<15} {summary2.num_atoms:<15}")
        print(f"{'Residues':<25} {summary.num_residues:<15} {summary2.num_residues:<15}")
        print(f"{'Chains':<25} {summary.num_chains:<15} {summary2.num_chains:<15}")
        print(f"{'Rg (A)':<25} {summary.radius_of_gyration:<15.2f} {summary2.radius_of_gyration:<15.2f}")
        print(f"{'Max CA dist (A)':<25} {summary.max_ca_distance:<15.2f} {summary2.max_ca_distance:<15.2f}")
        print(f"{'Hydrophobic ratio':<25} {summary.hydrophobic_ratio:<15.3f} {summary2.hydrophobic_ratio:<15.3f}")
        print(f"{'Analysis ready':<25} {str(summary.is_analysis_ready()):<15} {str(summary2.is_analysis_ready()):<15}")
    except Exception as e:
        print(f"Could not load second structure: {e}")

    # --- Practical Workflow ---
    print("\n9. PRACTICAL WORKFLOW: QUALITY FILTERING")
    print("-" * 40)

    print("""
Recommended workflow for analysis:

1. Load structure
2. Get quality report
3. Check is_analysis_ready()
4. If not ready, handle issues:
   - Multiple models: use first model or specify
   - Alternate locs: select primary conformer
   - CA-only: limited analysis possible
5. Get structure summary for metrics
6. Proceed with analysis
""")

    # Example code
    print("Example code:")
    print("""
    structure = pdbrust.parse_pdb_file(path)
    quality = structure.quality_report()

    if not quality.is_analysis_ready():
        if quality.has_multiple_models:
            print("Warning: Multiple models, using first")
        if quality.has_altlocs:
            print("Warning: Alternate locations present")
            # Could filter to primary conformer
        if quality.has_ca_only:
            print("Warning: CA-only structure")
            return

    # Proceed with analysis
    summary = structure.summary()
    print(f"Rg: {summary.radius_of_gyration:.2f} A")
""")

    # --- Summary Table ---
    print("\n10. SUMMARY")
    print("-" * 40)

    print(f"""
Structure: {PDB_FILE}

Quality Assessment:
  Analysis ready: {quality.is_analysis_ready()}
  Clean: {quality.is_clean()}

Key Metrics:
  Atoms: {quality.num_atoms}
  Residues: {quality.num_residues}
  Chains: {quality.num_chains}
  Radius of gyration: {summary.radius_of_gyration:.2f} A
  Max CA distance: {summary.max_ca_distance:.2f} A

Available methods:
  - quality_report(): Full quality assessment
  - structure_descriptors(): All computed metrics
  - summary(): Unified quality + descriptors
  - to_dict(): Dictionary serialization
  - to_csv_values(): CSV export
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
