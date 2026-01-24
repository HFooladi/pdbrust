#!/usr/bin/env python3
"""
Ramachandran analysis and backbone dihedral angles example for PDBRust Python bindings.

This example demonstrates:
- Phi/Psi dihedral angle calculation
- Ramachandran region classification
- Ramachandran statistics for structure validation
- Cis peptide bond detection
- Hydrogen bond network analysis
"""

import pdbrust

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Ramachandran & Dihedral Analysis Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms, {structure.num_residues} residues")

    # --- Phi/Psi Dihedral Angles ---
    print("\n1. PHI/PSI DIHEDRAL ANGLES")
    print("-" * 40)

    dihedrals = structure.phi_psi_angles()
    print(f"Total residues with dihedral data: {len(dihedrals)}")

    print(f"\n{'Chain':<6} {'Res#':<6} {'Name':<5} {'Phi':<10} {'Psi':<10} {'Omega':<10} {'Region'}")
    print("-" * 70)

    # Print first 15 residues
    for d in dihedrals[:15]:
        phi_str = f"{d.phi:.1f}" if d.phi is not None else "N/A"
        psi_str = f"{d.psi:.1f}" if d.psi is not None else "N/A"
        omega_str = f"{d.omega:.1f}" if d.omega is not None else "N/A"
        print(f"{d.chain_id:<6} {d.residue_seq:<6} {d.residue_name:<5} "
              f"{phi_str:<10} {psi_str:<10} {omega_str:<10} {d.ramachandran_region}")

    if len(dihedrals) > 15:
        print(f"... ({len(dihedrals) - 15} more residues)")

    # --- Dihedral Angle Reference ---
    print("\n2. DIHEDRAL ANGLE REFERENCE")
    print("-" * 40)

    print("""
Backbone Dihedral Angles:
  Phi (φ): C(i-1)-N(i)-CA(i)-C(i) rotation around N-CA bond
  Psi (ψ): N(i)-CA(i)-C(i)-N(i+1) rotation around CA-C bond
  Omega (ω): CA(i-1)-C(i-1)-N(i)-CA(i) peptide bond planarity

Typical Values:
  Alpha helix: φ ≈ -60°, ψ ≈ -45°
  Beta sheet:  φ ≈ -120°, ψ ≈ +130°
  Trans peptide: ω ≈ ±180°
  Cis peptide:   ω ≈ 0° (rare, mostly proline)
""")

    # --- Ramachandran Regions ---
    print("\n3. RAMACHANDRAN REGION REFERENCE")
    print("-" * 40)

    print("""
Ramachandran Regions:
  Core:     Favored region (most common conformations)
  Allowed:  Allowed region (less common but acceptable)
  Generous: Generously allowed (borderline)
  Outlier:  Unfavored region (potential modeling error)
  Glycine:  Special regions for glycine (no CB)
  Proline:  Special regions for proline (cyclic)
  PrePro:   Residues before proline (restricted)
  Unknown:  Cannot be classified (missing atoms)
""")

    # --- Ramachandran Statistics ---
    print("\n4. RAMACHANDRAN STATISTICS")
    print("-" * 40)

    stats = structure.ramachandran_statistics()

    print(f"Total residues analyzed: {stats.total_residues}")
    print(f"\nRamachandran Classification:")
    print(f"  Favored:  {stats.favored_count:4d} ({stats.favored_fraction * 100:.1f}%)")
    print(f"  Allowed:  {stats.allowed_count:4d} ({stats.allowed_fraction * 100:.1f}%)")
    print(f"  Outliers: {stats.outlier_count:4d} ({stats.outlier_fraction * 100:.1f}%)")

    # Quality assessment
    print(f"\nPeptide Bond Analysis:")
    print(f"  Cis peptide bonds: {stats.cis_peptide_count}")
    print(f"  Cis non-proline:   {stats.cis_nonpro_count}")

    # Interpretation
    print("\nQuality Interpretation:")
    if stats.favored_fraction > 0.95:
        print("  Excellent: >95% in favored regions")
    elif stats.favored_fraction > 0.90:
        print("  Good: >90% in favored regions")
    elif stats.favored_fraction > 0.80:
        print("  Acceptable: >80% in favored regions")
    else:
        print("  Poor: <80% in favored regions - potential issues")

    # --- Ramachandran Outliers ---
    print("\n5. RAMACHANDRAN OUTLIERS")
    print("-" * 40)

    outliers = structure.ramachandran_outliers()
    print(f"Found {len(outliers)} outlier residues")

    if outliers:
        print(f"\n{'Chain':<6} {'Res#':<6} {'Name':<5} {'Phi':<10} {'Psi':<10}")
        print("-" * 40)
        for d in outliers[:10]:
            phi_str = f"{d.phi:.1f}" if d.phi is not None else "N/A"
            psi_str = f"{d.psi:.1f}" if d.psi is not None else "N/A"
            print(f"{d.chain_id:<6} {d.residue_seq:<6} {d.residue_name:<5} "
                  f"{phi_str:<10} {psi_str:<10}")
        if len(outliers) > 10:
            print(f"... ({len(outliers) - 10} more)")

    # --- Cis Peptide Bonds ---
    print("\n6. CIS PEPTIDE BOND DETECTION")
    print("-" * 40)

    cis_residues = [d for d in dihedrals if d.is_cis_peptide()]
    print(f"Residues with cis peptide bond (omega ≈ 0°): {len(cis_residues)}")

    if cis_residues:
        for d in cis_residues:
            omega = d.omega if d.omega is not None else 0
            peptide_type = "cis-Pro" if d.residue_name == "PRO" else "cis non-Pro (unusual!)"
            print(f"  {d.chain_id}{d.residue_seq} {d.residue_name}: omega = {omega:.1f}° ({peptide_type})")
    else:
        print("  All peptide bonds are trans (normal)")

    # --- Hydrogen Bond Network Statistics ---
    print("\n7. HYDROGEN BOND NETWORK ANALYSIS")
    print("-" * 40)

    hbond_stats = structure.hbond_statistics()

    print(f"Total backbone H-bonds: {hbond_stats.total_hbonds}")
    print(f"\nH-bond Classification:")
    print(f"  Intra-helical (α-helix): {hbond_stats.intra_helical}")
    print(f"  Beta-sheet:              {hbond_stats.beta_sheet}")
    print(f"  Turns:                   {hbond_stats.turn}")
    print(f"  Long-range:              {hbond_stats.long_range}")
    print(f"  Inter-chain:             {hbond_stats.inter_chain}")

    print(f"\nH-bond Statistics:")
    print(f"  Mean energy: {hbond_stats.mean_energy:.2f} kcal/mol")
    print(f"  Donor residues: {hbond_stats.donor_residues}")
    print(f"  Acceptor residues: {hbond_stats.acceptor_residues}")

    # --- Per-Residue H-bonds ---
    print("\n8. PER-RESIDUE HYDROGEN BONDS")
    print("-" * 40)

    # Example: query H-bonds for residue 10
    res_hbonds = structure.hbonds_for_residue("A", 10)
    print(f"H-bonds for residue A10:")
    print(f"  Donated: {len(res_hbonds.donated)}")
    print(f"  Accepted: {len(res_hbonds.accepted)}")
    print(f"  Total: {res_hbonds.total()}")
    print(f"  Has H-bonds: {res_hbonds.has_hbonds()}")

    if res_hbonds.donated:
        print("\n  Donated H-bonds:")
        for hb in res_hbonds.donated:
            print(f"    -> {hb.acceptor_chain}{hb.acceptor_resid} {hb.acceptor_resname}: "
                  f"E={hb.energy:.2f} kcal/mol, d={hb.n_o_distance:.2f} A")

    if res_hbonds.accepted:
        print("\n  Accepted H-bonds:")
        for hb in res_hbonds.accepted:
            print(f"    <- {hb.donor_chain}{hb.donor_resid} {hb.donor_resname}: "
                  f"E={hb.energy:.2f} kcal/mol")

    # --- All Mainchain H-bonds Sample ---
    print("\n9. MAINCHAIN H-BOND DETAILS")
    print("-" * 40)

    all_hbonds = structure.mainchain_hbonds()
    print(f"Total mainchain H-bonds: {len(all_hbonds)}")

    if all_hbonds:
        print(f"\n{'Donor':<12} {'Acceptor':<12} {'Energy':<10} {'Type':<15} {'Strong?'}")
        print("-" * 55)

        for hb in all_hbonds[:10]:
            donor = f"{hb.donor_chain}{hb.donor_resid}"
            acceptor = f"{hb.acceptor_chain}{hb.acceptor_resid}"
            strong = "Yes" if hb.is_strong() else "No"
            print(f"{donor:<12} {acceptor:<12} {hb.energy:<10.2f} {hb.hbond_type!r:<15} {strong}")

        if len(all_hbonds) > 10:
            print(f"... ({len(all_hbonds) - 10} more)")

    # --- Summary ---
    print("\n10. STRUCTURE VALIDATION SUMMARY")
    print("-" * 40)

    print(f"""
Structure: {PDB_FILE}
Total residues: {len(dihedrals)}

Ramachandran Analysis:
  Favored:  {stats.favored_fraction * 100:.1f}%
  Allowed:  {stats.allowed_fraction * 100:.1f}%
  Outliers: {stats.outlier_fraction * 100:.1f}% ({stats.outlier_count} residues)

Peptide Bond Geometry:
  Cis peptides: {stats.cis_peptide_count}
  Cis non-proline: {stats.cis_nonpro_count}

Hydrogen Bond Network:
  Total H-bonds: {hbond_stats.total_hbonds}
  Helical: {hbond_stats.intra_helical}, Sheet: {hbond_stats.beta_sheet}
  Mean energy: {hbond_stats.mean_energy:.2f} kcal/mol
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
