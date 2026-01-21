#!/usr/bin/env python3
"""
Batch processing example for PDBRust Python bindings.

This example demonstrates:
- Processing multiple PDB/mmCIF files in a loop
- Computing summaries for each structure
- Quality filtering during batch processing
- Exporting results to CSV format
- Error handling for failed parses
- Computing statistics across a dataset
"""

import os
import pdbrust
from pdbrust import StructureSummary

# Directory containing sample PDB files
PDB_DIR = "../../examples/pdb_files"


def main():
    print("=" * 60)
    print("PDBRust Batch Processing Example")
    print("=" * 60)

    # --- Find All PDB/mmCIF Files ---
    print("\n1. DISCOVERING FILES")
    print("-" * 40)

    # Find all parseable files
    pdb_files = []
    for filename in os.listdir(PDB_DIR):
        if filename.endswith(('.pdb', '.cif')):
            pdb_files.append(os.path.join(PDB_DIR, filename))

    print(f"Found {len(pdb_files)} structure files:")
    for f in pdb_files:
        print(f"  - {os.path.basename(f)}")

    # --- Process All Files ---
    print("\n2. PROCESSING FILES")
    print("-" * 40)

    results = []
    failed = []

    for filepath in pdb_files:
        filename = os.path.basename(filepath)
        print(f"Processing {filename}...", end=" ")

        try:
            # Parse structure (auto-detect format)
            structure = pdbrust.parse_structure_file(filepath)

            # Get summary
            summary = structure.summary()
            quality = structure.quality_report()

            results.append({
                'filename': filename,
                'structure': structure,
                'summary': summary,
                'quality': quality
            })

            print(f"OK ({structure.num_atoms} atoms)")

        except Exception as e:
            failed.append({'filename': filename, 'error': str(e)})
            print(f"FAILED: {e}")

    print(f"\nProcessed: {len(results)} successful, {len(failed)} failed")

    # --- Quality Filtering ---
    print("\n3. QUALITY FILTERING")
    print("-" * 40)

    # Filter for analysis-ready structures
    analysis_ready = [r for r in results if r['quality'].is_analysis_ready()]
    print(f"Analysis-ready structures: {len(analysis_ready)}/{len(results)}")

    # Filter for clean structures (no altlocs, not CA-only)
    clean = [r for r in results if r['quality'].is_clean()]
    print(f"Clean structures: {len(clean)}/{len(results)}")

    # Show why structures were filtered out
    not_ready = [r for r in results if not r['quality'].is_analysis_ready()]
    if not_ready:
        print("\nStructures not analysis-ready:")
        for r in not_ready:
            q = r['quality']
            reasons = []
            if q.has_multiple_models:
                reasons.append("multiple models")
            if q.has_altlocs:
                reasons.append("altlocs")
            if q.has_ca_only:
                reasons.append("CA-only")
            print(f"  {r['filename']}: {', '.join(reasons)}")

    # --- Results Summary Table ---
    print("\n4. RESULTS SUMMARY")
    print("-" * 40)

    # Print header
    header = f"{'Filename':<20} {'Atoms':>8} {'Residues':>8} {'Chains':>6} {'Rg':>8} {'Ready':>6}"
    print(header)
    print("-" * len(header))

    for r in results:
        s = r['summary']
        q = r['quality']
        ready = "Yes" if q.is_analysis_ready() else "No"
        print(f"{r['filename']:<20} {s.num_atoms:>8} {s.num_residues:>8} "
              f"{s.num_chains:>6} {s.radius_of_gyration:>8.2f} {ready:>6}")

    # --- Dataset Statistics ---
    print("\n5. DATASET STATISTICS")
    print("-" * 40)

    if results:
        # Compute statistics across all structures
        all_atoms = [r['summary'].num_atoms for r in results]
        all_residues = [r['summary'].num_residues for r in results]
        all_rg = [r['summary'].radius_of_gyration for r in results]

        print(f"Number of structures: {len(results)}")
        print(f"\nAtom counts:")
        print(f"  Min: {min(all_atoms)}")
        print(f"  Max: {max(all_atoms)}")
        print(f"  Mean: {sum(all_atoms) / len(all_atoms):.1f}")
        print(f"  Total: {sum(all_atoms)}")

        print(f"\nResidue counts:")
        print(f"  Min: {min(all_residues)}")
        print(f"  Max: {max(all_residues)}")
        print(f"  Mean: {sum(all_residues) / len(all_residues):.1f}")

        print(f"\nRadius of gyration:")
        print(f"  Min: {min(all_rg):.2f} A")
        print(f"  Max: {max(all_rg):.2f} A")
        print(f"  Mean: {sum(all_rg) / len(all_rg):.2f} A")

        # Count structures with specific features
        with_hetero = sum(1 for r in results if r['quality'].has_hetatm)
        with_ssbonds = sum(1 for r in results if r['quality'].has_ssbonds)
        with_hydrogens = sum(1 for r in results if r['quality'].has_hydrogens)

        print(f"\nFeature counts:")
        print(f"  With HETATM: {with_hetero}/{len(results)}")
        print(f"  With SS bonds: {with_ssbonds}/{len(results)}")
        print(f"  With hydrogens: {with_hydrogens}/{len(results)}")

    # --- CSV Export ---
    print("\n6. CSV EXPORT")
    print("-" * 40)

    # Get CSV header from StructureSummary
    field_names = StructureSummary.field_names()
    csv_header = "filename," + ",".join(field_names)

    print("CSV format:")
    print(f"  {len(field_names) + 1} columns")
    print(f"  Header: filename,{','.join(field_names[:5])}...")

    # Generate CSV content
    csv_lines = [csv_header]
    for r in results:
        values = r['summary'].to_csv_values()
        line = r['filename'] + "," + ",".join(values)
        csv_lines.append(line)

    # Show first few lines
    print("\nCSV content (first 3 rows):")
    for line in csv_lines[:4]:
        # Truncate long lines
        display = line if len(line) < 80 else line[:77] + "..."
        print(f"  {display}")

    # Could write to file:
    # with open("batch_results.csv", "w") as f:
    #     f.write("\n".join(csv_lines))

    # --- Per-Structure Details ---
    print("\n7. DETAILED ANALYSIS (Selected Structures)")
    print("-" * 40)

    # Analyze a few structures in detail
    for r in results[:2]:
        s = r['summary']
        q = r['quality']

        print(f"\n{r['filename']}:")
        print(f"  Size: {s.num_atoms} atoms, {s.num_residues} residues, {s.num_chains} chains")
        print(f"  Geometry: Rg={s.radius_of_gyration:.2f} A, max_dist={s.max_ca_distance:.2f} A")
        print(f"  Composition: hydrophobic={s.hydrophobic_ratio:.2f}, glycine={s.glycine_ratio:.2f}")
        print(f"  Quality: analysis_ready={q.is_analysis_ready()}, clean={q.is_clean()}")

    # --- Error Handling ---
    print("\n8. HANDLING ERRORS")
    print("-" * 40)

    if failed:
        print(f"Failed to parse {len(failed)} files:")
        for f in failed:
            print(f"  {f['filename']}: {f['error']}")
    else:
        print("All files parsed successfully!")

    print("\nRecommended error handling pattern:")
    print("""
    for filepath in pdb_files:
        try:
            structure = pdbrust.parse_structure_file(filepath)
            # Process structure...
        except Exception as e:
            # Log error and continue
            print(f"Failed to parse {filepath}: {e}")
            continue
""")

    # --- Filtering Workflow ---
    print("\n9. COMPLETE FILTERING WORKFLOW")
    print("-" * 40)

    # Example: Find structures suitable for homology comparison
    suitable = []
    for r in results:
        q = r['quality']
        s = r['summary']

        # Apply quality criteria
        if not q.is_analysis_ready():
            continue

        # Apply size criteria
        if s.num_residues < 50:
            continue  # Too small

        # Apply resolution criteria (if available)
        resolution = r['structure'].get_resolution()
        if resolution is not None and resolution > 3.0:
            continue  # Too low resolution

        suitable.append(r)

    print(f"Structures meeting all criteria: {len(suitable)}/{len(results)}")
    print("\nCriteria applied:")
    print("  - Analysis ready (single model, no altlocs)")
    print("  - At least 50 residues")
    print("  - Resolution <= 3.0 A (if available)")

    if suitable:
        print("\nSuitable structures:")
        for r in suitable:
            res = r['structure'].get_resolution()
            res_str = f"{res:.2f} A" if res else "N/A"
            print(f"  {r['filename']}: {r['summary'].num_residues} residues, resolution: {res_str}")

    # --- Summary ---
    print("\n10. SUMMARY")
    print("-" * 40)

    print(f"""
Batch Processing Complete:
  Total files found: {len(pdb_files)}
  Successfully parsed: {len(results)}
  Failed: {len(failed)}

Quality Breakdown:
  Analysis-ready: {len(analysis_ready)}
  Clean: {len(clean)}

Key methods used:
  - pdbrust.parse_structure_file(): Auto-detect format
  - structure.summary(): Get unified summary
  - structure.quality_report(): Get quality assessment
  - summary.to_csv_values(): CSV export
  - StructureSummary.field_names(): CSV header

This workflow is suitable for:
  - Dataset curation
  - Quality control
  - Feature extraction
  - Batch analysis pipelines
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
