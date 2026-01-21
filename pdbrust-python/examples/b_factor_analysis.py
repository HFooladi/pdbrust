#!/usr/bin/env python3
"""
B-factor (temperature factor) analysis example for PDBRust Python bindings.

This example demonstrates:
- B-factor statistics (mean, min, max, std)
- CA-specific B-factor analysis
- Per-residue B-factor profile
- Identifying flexible and rigid regions
- B-factor normalization
- Percentile calculations
"""

import pdbrust

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust B-factor Analysis Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_residues} residues")

    # --- Basic B-factor Statistics ---
    print("\n1. BASIC B-FACTOR STATISTICS")
    print("-" * 40)

    # All atom statistics
    mean_b = structure.b_factor_mean()
    min_b = structure.b_factor_min()
    max_b = structure.b_factor_max()
    std_b = structure.b_factor_std()

    print(f"All atoms:")
    print(f"  Mean B-factor: {mean_b:.2f} A^2")
    print(f"  Min B-factor:  {min_b:.2f} A^2")
    print(f"  Max B-factor:  {max_b:.2f} A^2")
    print(f"  Std deviation: {std_b:.2f} A^2")

    # CA-specific statistics
    mean_b_ca = structure.b_factor_mean_ca()
    print(f"\nCA atoms only:")
    print(f"  Mean B-factor: {mean_b_ca:.2f} A^2")

    # --- Per-Residue B-factor Profile ---
    print("\n2. PER-RESIDUE B-FACTOR PROFILE")
    print("-" * 40)

    profile = structure.b_factor_profile()
    print(f"Total residues in profile: {len(profile)}")

    print(f"\n{'Chain':<6} {'ResSeq':<8} {'Name':<5} {'Mean B':<10} {'Min B':<10} {'Max B':<10}")
    print("-" * 55)

    # Print first 15 residues
    for res in profile[:15]:
        print(f"{res.chain_id:<6} {res.residue_seq:<8} {res.residue_name:<5} "
              f"{res.b_factor_mean:<10.2f} {res.b_factor_min:<10.2f} {res.b_factor_max:<10.2f}")

    print(f"... ({len(profile) - 15} more residues)")

    # --- ResidueBFactor Object ---
    print("\n3. RESIDUEBFACTOR DETAILS")
    print("-" * 40)

    # Get details for first residue
    first_res = profile[0]
    print(f"First residue details:")
    print(f"  Chain ID: {first_res.chain_id}")
    print(f"  Residue seq: {first_res.residue_seq}")
    print(f"  Residue name: {first_res.residue_name}")
    print(f"  Insertion code: {first_res.ins_code}")
    print(f"  Mean B-factor: {first_res.b_factor_mean:.2f} A^2")
    print(f"  Min B-factor: {first_res.b_factor_min:.2f} A^2")
    print(f"  Max B-factor: {first_res.b_factor_max:.2f} A^2")
    print(f"  Atom count: {first_res.atom_count}")

    # Convert to dictionary
    res_dict = first_res.to_dict()
    print(f"\nAs dictionary: {list(res_dict.keys())}")

    # --- Identifying Flexible Residues ---
    print("\n4. FLEXIBLE RESIDUES (HIGH B-FACTORS)")
    print("-" * 40)

    # Find residues with mean B-factor > 30 A^2
    threshold_high = 30.0
    flexible = structure.flexible_residues(threshold_high)

    print(f"Residues with mean B > {threshold_high} A^2: {len(flexible)}")
    if flexible:
        print(f"\n{'Chain':<6} {'ResSeq':<8} {'Name':<5} {'Mean B':<10}")
        print("-" * 32)
        for res in flexible[:10]:
            print(f"{res.chain_id:<6} {res.residue_seq:<8} {res.residue_name:<5} {res.b_factor_mean:<10.2f}")
        if len(flexible) > 10:
            print(f"... and {len(flexible) - 10} more")

    # Find very flexible residues
    threshold_very_high = 40.0
    very_flexible = structure.flexible_residues(threshold_very_high)
    print(f"\nVery flexible (B > {threshold_very_high} A^2): {len(very_flexible)} residues")

    # --- Identifying Rigid Residues ---
    print("\n5. RIGID RESIDUES (LOW B-FACTORS)")
    print("-" * 40)

    # Find residues with mean B-factor < 15 A^2
    threshold_low = 15.0
    rigid = structure.rigid_residues(threshold_low)

    print(f"Residues with mean B < {threshold_low} A^2: {len(rigid)}")
    if rigid:
        print(f"\n{'Chain':<6} {'ResSeq':<8} {'Name':<5} {'Mean B':<10}")
        print("-" * 32)
        for res in rigid[:10]:
            print(f"{res.chain_id:<6} {res.residue_seq:<8} {res.residue_name:<5} {res.b_factor_mean:<10.2f}")
        if len(rigid) > 10:
            print(f"... and {len(rigid) - 10} more")

    # --- B-factor Normalization ---
    print("\n6. B-FACTOR NORMALIZATION (Z-SCORE)")
    print("-" * 40)

    # Normalize B-factors to Z-scores (mean=0, std=1)
    normalized = structure.normalize_b_factors()

    print("Original structure:")
    print(f"  Mean B: {structure.b_factor_mean():.4f}")
    print(f"  Std B:  {structure.b_factor_std():.4f}")

    print("\nNormalized structure:")
    print(f"  Mean B: {normalized.b_factor_mean():.4f} (should be ~0)")
    print(f"  Std B:  {normalized.b_factor_std():.4f} (should be ~1)")

    # Compare profiles
    norm_profile = normalized.b_factor_profile()
    print(f"\nNormalized B-factors (first 5 residues):")
    print(f"  Original -> Normalized")
    for orig, norm in zip(profile[:5], norm_profile[:5]):
        print(f"  {orig.residue_name}{orig.residue_seq}: "
              f"{orig.b_factor_mean:6.2f} -> {norm.b_factor_mean:+6.3f}")

    # --- B-factor Percentile ---
    print("\n7. B-FACTOR PERCENTILE")
    print("-" * 40)

    # Get atoms to check percentiles
    atoms = structure.atoms
    print(f"Total atoms: {len(atoms)}")

    # Check percentile for specific atoms
    atom_serials_to_check = [1, 10, 50, 100]
    for serial in atom_serials_to_check:
        percentile = structure.b_factor_percentile(serial)
        if percentile is not None:
            # Find the atom to get its actual B-factor
            atom = next((a for a in atoms if a.serial == serial), None)
            if atom:
                print(f"Atom {serial:4d} ({atom.name:4s}): "
                      f"B = {atom.temp_factor:6.2f} -> {percentile:5.1f}th percentile")
        else:
            print(f"Atom {serial}: not found")

    # --- Distribution Analysis ---
    print("\n8. B-FACTOR DISTRIBUTION")
    print("-" * 40)

    # Categorize residues by B-factor ranges
    ranges = [
        (0, 15, "Very low (well-ordered)"),
        (15, 25, "Low"),
        (25, 35, "Medium"),
        (35, 50, "High"),
        (50, float('inf'), "Very high (flexible)")
    ]

    print("B-factor distribution by category:")
    for low, high, label in ranges:
        count = sum(1 for r in profile if low <= r.b_factor_mean < high)
        pct = count / len(profile) * 100 if profile else 0
        print(f"  {label:25s}: {count:3d} residues ({pct:5.1f}%)")

    # --- Structural Interpretation ---
    print("\n9. STRUCTURAL INTERPRETATION")
    print("-" * 40)

    # Find highest B-factor residues (potential loop regions)
    sorted_by_b = sorted(profile, key=lambda r: r.b_factor_mean, reverse=True)

    print("Top 5 highest B-factor residues (likely in loops):")
    for res in sorted_by_b[:5]:
        print(f"  {res.chain_id}{res.residue_seq:3d} {res.residue_name}: {res.b_factor_mean:.2f} A^2")

    # Find lowest B-factor residues (likely in core)
    print("\nTop 5 lowest B-factor residues (likely in core):")
    for res in sorted_by_b[-5:]:
        print(f"  {res.chain_id}{res.residue_seq:3d} {res.residue_name}: {res.b_factor_mean:.2f} A^2")

    # --- Summary ---
    print("\n10. SUMMARY")
    print("-" * 40)

    total_residues = len(profile)
    flexible_count = len(flexible)
    rigid_count = len(rigid)
    moderate_count = total_residues - flexible_count - rigid_count

    print(f"""
Structure: {PDB_FILE}
Total residues: {total_residues}

B-factor Statistics:
  Mean:     {mean_b:.2f} A^2
  Std:      {std_b:.2f} A^2
  Range:    {min_b:.2f} - {max_b:.2f} A^2

Flexibility Classification:
  Rigid (B < {threshold_low}):   {rigid_count:3d} residues ({rigid_count/total_residues*100:5.1f}%)
  Moderate:              {moderate_count:3d} residues ({moderate_count/total_residues*100:5.1f}%)
  Flexible (B > {threshold_high}): {flexible_count:3d} residues ({flexible_count/total_residues*100:5.1f}%)

Interpretation:
  Low B-factors indicate well-ordered, rigid regions.
  High B-factors indicate mobile, flexible, or disordered regions.
  Terminal residues and loop regions typically have higher B-factors.
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
