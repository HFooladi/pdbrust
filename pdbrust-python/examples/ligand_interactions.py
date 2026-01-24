#!/usr/bin/env python3
"""
Protein-ligand interaction analysis example for PDBRust Python bindings.

This example demonstrates:
- Binding site identification
- Ligand interaction profiling
- Detection of H-bonds, salt bridges, and hydrophobic contacts
- Analyzing all ligands in a structure
"""

import pdbrust

# We'll create a synthetic protein-ligand structure for demonstration
# In practice, you would use a PDB file with a bound ligand


def create_protein_ligand_pdb() -> str:
    """Create a synthetic PDB string with protein and ATP ligand for demonstration."""
    return """HEADER    PROTEIN-LIGAND COMPLEX
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       1.950  -0.700   1.250  1.00 20.00           C
ATOM      6  N   LYS A   2       3.320   1.567   0.000  1.00 20.00           N
ATOM      7  CA  LYS A   2       3.954   2.881   0.000  1.00 20.00           C
ATOM      8  C   LYS A   2       5.464   2.771   0.000  1.00 20.00           C
ATOM      9  O   LYS A   2       6.108   1.726   0.000  1.00 20.00           O
ATOM     10  CB  LYS A   2       3.500   3.700   1.200  1.00 20.00           C
ATOM     11  NZ  LYS A   2       4.500   4.000   2.000  1.00 20.00           N
ATOM     12  N   SER A   3       6.012   3.973   0.000  1.00 20.00           N
ATOM     13  CA  SER A   3       7.445   4.145   0.000  1.00 20.00           C
ATOM     14  C   SER A   3       8.115   2.811   0.000  1.00 20.00           C
ATOM     15  O   SER A   3       7.418   1.800   0.000  1.00 20.00           O
ATOM     16  OG  SER A   3       7.500   5.000   1.000  1.00 20.00           O
ATOM     17  N   ASP A   4       9.400   2.800  -0.300  1.00 20.00           N
ATOM     18  CA  ASP A   4      10.200   1.600  -0.300  1.00 20.00           C
ATOM     19  C   ASP A   4      11.500   1.800   0.400  1.00 20.00           C
ATOM     20  O   ASP A   4      11.600   2.700   1.250  1.00 20.00           O
ATOM     21  OD1 ASP A   4       9.500   0.500  -1.200  1.00 20.00           O
ATOM     22  OD2 ASP A   4      10.000  -0.200  -0.200  1.00 20.00           O
ATOM     23  N   PHE A   5      12.500   1.000   0.100  1.00 20.00           N
ATOM     24  CA  PHE A   5      13.800   1.100   0.700  1.00 20.00           C
ATOM     25  C   PHE A   5      14.700   0.000   0.200  1.00 20.00           C
ATOM     26  O   PHE A   5      14.300  -1.150   0.100  1.00 20.00           O
ATOM     27  CB  PHE A   5      14.400   2.500   0.500  1.00 20.00           C
ATOM     28  CG  PHE A   5      13.500   3.600   1.000  1.00 20.00           C
HETATM   29  PA  ATP A 100       5.000   5.000   1.500  1.00 20.00           P
HETATM   30  O1A ATP A 100       4.000   5.500   2.500  1.00 20.00           O
HETATM   31  O2A ATP A 100       5.500   5.800   0.500  1.00 20.00           O
HETATM   32  O3A ATP A 100       6.000   4.000   2.000  1.00 20.00           O
HETATM   33  N1  ATP A 100       5.500   4.000   0.500  1.00 20.00           N
HETATM   34  C2  ATP A 100       4.500   3.500   0.000  1.00 20.00           C
HETATM   35  C4  ATP A 100       6.000   3.000   1.000  1.00 20.00           C
HETATM   36  C5  ATP A 100       5.000   2.500   1.500  1.00 20.00           C
END
"""


def main():
    print("=" * 60)
    print("PDBRust Protein-Ligand Interaction Analysis Example")
    print("=" * 60)

    # Create and parse the synthetic structure
    pdb_string = create_protein_ligand_pdb()
    structure = pdbrust.parse_pdb_string(pdb_string)
    print(f"Loaded structure: {structure.num_atoms} atoms")

    # Count protein vs ligand atoms
    protein_atoms = [a for a in structure.atoms if a.residue_name in ['ALA', 'LYS', 'SER', 'ASP', 'PHE']]
    ligand_atoms = [a for a in structure.atoms if a.residue_name == 'ATP']
    print(f"  Protein atoms: {len(protein_atoms)}")
    print(f"  Ligand atoms (ATP): {len(ligand_atoms)}")

    # --- Binding Site Identification ---
    print("\n1. BINDING SITE IDENTIFICATION")
    print("-" * 40)

    # Find residues within 6 Angstroms of ATP
    site = structure.binding_site("ATP", 6.0)

    if site:
        print(f"Binding site for {site.ligand_name}:")
        print(f"  Ligand: {site.ligand_name} (Chain {site.ligand_chain}, Residue {site.ligand_resid})")
        print(f"  Distance cutoff: {site.distance_cutoff} Angstroms")
        print(f"  Contact residues: {site.num_residues()}")

        print(f"\n{'Chain':<6} {'Res#':<6} {'Name':<5} {'Min Dist':<10} {'Contacts':<10}")
        print("-" * 45)

        for res in site.contact_residues:
            print(f"{res.chain_id:<6} {res.residue_seq:<6} {res.residue_name:<5} "
                  f"{res.min_distance:<10.2f} {res.num_contacts:<10}")
    else:
        print("No binding site found for ATP")

    # Try different distance cutoffs
    print("\n  Effect of distance cutoff:")
    for cutoff in [4.0, 5.0, 6.0, 8.0]:
        site_test = structure.binding_site("ATP", cutoff)
        if site_test:
            print(f"    {cutoff} A: {site_test.num_residues()} residues")
        else:
            print(f"    {cutoff} A: 0 residues")

    # --- Ligand Interaction Profile ---
    print("\n2. LIGAND INTERACTION PROFILE")
    print("-" * 40)

    profile = structure.ligand_interactions("ATP")

    if profile:
        print(f"Interaction profile for {profile.ligand_name}:")
        print(f"  Chain: {profile.ligand_chain}, Residue ID: {profile.ligand_resid}")
        print(f"\nInteraction Summary:")
        print(f"  Contact residues: {len(profile.contact_residues)}")
        print(f"  Hydrogen bonds: {len(profile.hydrogen_bonds)}")
        print(f"  Salt bridges: {len(profile.salt_bridges)}")
        print(f"  Hydrophobic contacts: {len(profile.hydrophobic_contacts)}")
        print(f"  Total interactions: {profile.total_interactions()}")
        print(f"  Has interactions: {profile.has_interactions()}")
    else:
        print("No interaction profile found for ATP")

    # --- Hydrogen Bond Details ---
    print("\n3. HYDROGEN BOND INTERACTIONS")
    print("-" * 40)

    if profile and profile.hydrogen_bonds:
        print(f"Found {len(profile.hydrogen_bonds)} protein-ligand H-bonds:")
        print(f"\n{'Protein':<20} {'Ligand Atom':<12} {'Distance':<10} {'Donor?'}")
        print("-" * 50)

        for hb in profile.hydrogen_bonds:
            protein_info = f"{hb.protein_resname}{hb.protein_resid}:{hb.protein_atom}"
            donor = "Protein" if hb.is_protein_donor else "Ligand"
            print(f"{protein_info:<20} {hb.ligand_atom:<12} {hb.distance:<10.2f} {donor}")
    else:
        print("No hydrogen bonds detected")
        print("  (H-bonds require donor-acceptor distance <= 3.5 A)")

    # --- Salt Bridge Details ---
    print("\n4. SALT BRIDGE INTERACTIONS")
    print("-" * 40)

    if profile and profile.salt_bridges:
        print(f"Found {len(profile.salt_bridges)} salt bridges:")
        print(f"\n{'Protein':<20} {'Ligand Atom':<12} {'Distance':<10} {'Charge'}")
        print("-" * 50)

        for sb in profile.salt_bridges:
            protein_info = f"{sb.protein_resname}{sb.protein_resid}:{sb.protein_atom}"
            charge = "Protein +" if sb.protein_positive else "Protein -"
            print(f"{protein_info:<20} {sb.ligand_atom:<12} {sb.distance:<10.2f} {charge}")
    else:
        print("No salt bridges detected")
        print("  (Salt bridges: charged atoms within 4.0 A)")
        print("  (Charged residues: LYS-NZ, ARG-NH, ASP-OD, GLU-OE)")

    # --- Hydrophobic Contact Details ---
    print("\n5. HYDROPHOBIC CONTACTS")
    print("-" * 40)

    if profile and profile.hydrophobic_contacts:
        print(f"Found {len(profile.hydrophobic_contacts)} hydrophobic contacts:")
        print(f"\n{'Protein':<20} {'Ligand Atom':<12} {'Distance':<10}")
        print("-" * 45)

        for hc in profile.hydrophobic_contacts[:10]:
            protein_info = f"{hc.protein_resname}{hc.protein_resid}:{hc.protein_atom}"
            print(f"{protein_info:<20} {hc.ligand_atom:<12} {hc.distance:<10.2f}")

        if len(profile.hydrophobic_contacts) > 10:
            print(f"... ({len(profile.hydrophobic_contacts) - 10} more)")
    else:
        print("No hydrophobic contacts detected")
        print("  (Hydrophobic contacts: non-polar atoms within 4.0 A)")

    # --- All Ligand Interactions ---
    print("\n6. ANALYZING ALL LIGANDS")
    print("-" * 40)

    all_profiles = structure.all_ligand_interactions()
    print(f"Found {len(all_profiles)} ligand(s) in structure:")

    for p in all_profiles:
        print(f"\n  {p.ligand_name} (Chain {p.ligand_chain}, Res {p.ligand_resid}):")
        print(f"    Contacts: {len(p.contact_residues)}")
        print(f"    H-bonds: {len(p.hydrogen_bonds)}")
        print(f"    Salt bridges: {len(p.salt_bridges)}")
        print(f"    Hydrophobic: {len(p.hydrophobic_contacts)}")
        print(f"    Total: {p.total_interactions()}")

    # --- Interaction Detection Reference ---
    print("\n7. INTERACTION DETECTION THRESHOLDS")
    print("-" * 40)

    print("""
Distance Thresholds for Interaction Detection:

  Hydrogen Bonds:
    - Distance: <= 3.5 Angstroms (heavy atom distance)
    - Donor atoms: N, O (with attached H)
    - Acceptor atoms: N, O

  Salt Bridges:
    - Distance: <= 4.0 Angstroms
    - Protein positive: LYS (NZ), ARG (NH1, NH2), HIS (ND, NE)
    - Protein negative: ASP (OD1, OD2), GLU (OE1, OE2)

  Hydrophobic Contacts:
    - Distance: <= 4.0 Angstroms
    - Hydrophobic residues: ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO
    - Carbon atoms in non-polar environment
""")

    # --- Practical Usage ---
    print("\n8. PRACTICAL USAGE EXAMPLE")
    print("-" * 40)

    print("""
# For real protein-ligand structures, load from PDB file:

structure = pdbrust.parse_pdb_file("protein_ligand.pdb")

# Find all ligands (HETATM records, excluding HOH/common ions)
profiles = structure.all_ligand_interactions()

# Or analyze a specific ligand
site = structure.binding_site("ATP", 5.0)  # 5 Angstrom cutoff
profile = structure.ligand_interactions("ATP")

# Filter key residues for drug design
if site:
    key_residues = [r for r in site.contact_residues if r.min_distance < 4.0]
    print(f"Key binding site residues: {len(key_residues)}")
""")

    # --- Summary ---
    print("\n9. SUMMARY")
    print("-" * 40)

    if profile:
        print(f"""
Ligand: {profile.ligand_name}
Location: Chain {profile.ligand_chain}, Residue {profile.ligand_resid}

Binding Site ({site.distance_cutoff if site else 'N/A'} A cutoff):
  Contact residues: {site.num_residues() if site else 0}

Interaction Profile:
  Hydrogen bonds:      {len(profile.hydrogen_bonds)}
  Salt bridges:        {len(profile.salt_bridges)}
  Hydrophobic contacts: {len(profile.hydrophobic_contacts)}
  Total interactions:  {profile.total_interactions()}
""")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
