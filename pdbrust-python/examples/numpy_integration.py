#!/usr/bin/env python3
"""
Numpy integration example for PDBRust Python bindings.

This example demonstrates:
- Getting coordinates as numpy arrays
- Distance matrix calculation
- Contact map generation
- Integration with scientific computing workflows
"""

import pdbrust
import numpy as np

# Sample PDB file
PDB_FILE = "../../examples/pdb_files/1UBQ.pdb"


def main():
    print("=" * 60)
    print("PDBRust Numpy Integration Example")
    print("=" * 60)

    # Load structure
    structure = pdbrust.parse_pdb_file(PDB_FILE)
    print(f"Loaded structure: {structure.num_atoms} atoms")

    # --- Coordinate Arrays ---
    print("\n1. COORDINATE ARRAYS")
    print("-" * 40)

    # All atom coordinates
    all_coords = structure.get_coords_array()
    print(f"All atoms: shape = {all_coords.shape}, dtype = {all_coords.dtype}")
    print(f"  First atom: [{all_coords[0, 0]:.3f}, {all_coords[0, 1]:.3f}, {all_coords[0, 2]:.3f}]")

    # CA atom coordinates
    ca_coords = structure.get_ca_coords_array()
    print(f"CA atoms: shape = {ca_coords.shape}")

    # Backbone coordinates (N, CA, C, O)
    bb_coords = structure.get_backbone_coords_array()
    print(f"Backbone atoms: shape = {bb_coords.shape}")

    # Chain-specific CA coordinates
    chains = structure.get_chain_ids()
    if chains:
        chain_ca = structure.get_ca_coords_array(chains[0])
        print(f"Chain {chains[0]} CA atoms: shape = {chain_ca.shape}")

    # --- Distance Matrix ---
    print("\n2. DISTANCE MATRIX")
    print("-" * 40)

    # CA-CA distance matrix
    dist_matrix = structure.distance_matrix_ca()
    print(f"CA distance matrix: shape = {dist_matrix.shape}")
    print(f"  Min distance: {dist_matrix[dist_matrix > 0].min():.2f} Å")
    print(f"  Max distance: {dist_matrix.max():.2f} Å")
    print(f"  Mean distance: {dist_matrix.mean():.2f} Å")

    # Show a small portion of the matrix
    print("\nDistance matrix (first 5x5):")
    print("     ", end="")
    for j in range(min(5, dist_matrix.shape[1])):
        print(f"{j+1:>7}", end="")
    print()
    for i in range(min(5, dist_matrix.shape[0])):
        print(f"{i+1:>4} ", end="")
        for j in range(min(5, dist_matrix.shape[1])):
            print(f"{dist_matrix[i, j]:>7.2f}", end="")
        print()

    # All-atom distance matrix (can be large!)
    # Uncomment if needed:
    # all_dist = structure.distance_matrix()
    # print(f"All-atom distance matrix: shape = {all_dist.shape}")

    # --- Contact Map ---
    print("\n3. CONTACT MAP")
    print("-" * 40)

    # CA contact map with default threshold (8 Å)
    contact_map = structure.contact_map_ca()
    print(f"CA contact map (8 Å): shape = {contact_map.shape}, dtype = {contact_map.dtype}")
    print(f"  Total contacts: {contact_map.sum()}")
    print(f"  Contact density: {contact_map.sum() / contact_map.size:.4f}")

    # Contact map with custom threshold
    contact_map_6 = structure.contact_map_ca(threshold=6.0)
    print(f"CA contact map (6 Å): {contact_map_6.sum()} contacts")

    contact_map_10 = structure.contact_map_ca(threshold=10.0)
    print(f"CA contact map (10 Å): {contact_map_10.sum()} contacts")

    # All-atom contact map
    all_contacts = structure.contact_map(threshold=4.5)
    print(f"All-atom contact map (4.5 Å): {all_contacts.sum()} contacts")

    # --- Numpy Operations ---
    print("\n4. NUMPY OPERATIONS")
    print("-" * 40)

    # Calculate center of mass using numpy
    com = np.mean(all_coords, axis=0)
    print(f"Center of mass (numpy): [{com[0]:.3f}, {com[1]:.3f}, {com[2]:.3f}]")

    # Compare with pdbrust method
    pdbrust_com = structure.get_centroid()
    print(f"Center of mass (pdbrust): [{pdbrust_com[0]:.3f}, {pdbrust_com[1]:.3f}, {pdbrust_com[2]:.3f}]")

    # Radius of gyration calculation
    centered = ca_coords - np.mean(ca_coords, axis=0)
    rg_numpy = np.sqrt(np.mean(np.sum(centered**2, axis=1)))
    rg_pdbrust = structure.radius_of_gyration()
    print(f"\nRadius of gyration:")
    print(f"  Numpy: {rg_numpy:.3f} Å")
    print(f"  PDBRust: {rg_pdbrust:.3f} Å")

    # --- Use Cases ---
    print("\n5. USE CASES")
    print("-" * 40)

    # Use case 1: Find long-range contacts
    n_residues = ca_coords.shape[0]
    long_range = np.zeros_like(contact_map)
    for i in range(n_residues):
        for j in range(i + 12, n_residues):  # At least 12 residues apart
            long_range[i, j] = contact_map[i, j]
            long_range[j, i] = contact_map[i, j]
    print(f"Long-range contacts (>12 residues apart): {long_range.sum() // 2}")

    # Use case 2: Contact order calculation
    if contact_map.sum() > 0:
        contact_indices = np.where(contact_map)
        seq_separations = np.abs(contact_indices[0] - contact_indices[1])
        contact_order = np.mean(seq_separations[seq_separations > 0]) / n_residues
        print(f"Relative contact order: {contact_order:.4f}")

    # Use case 3: Pairwise distances between specific residues
    # Get distances between residue 1 and all others
    distances_from_1 = dist_matrix[0, :]
    print(f"\nDistances from residue 1:")
    print(f"  To residue 10: {distances_from_1[9]:.2f} Å")
    print(f"  To last residue: {distances_from_1[-1]:.2f} Å")

    # --- Machine Learning Ready ---
    print("\n6. MACHINE LEARNING READY")
    print("-" * 40)

    # Contact map as float for ML models
    contact_float = contact_map.astype(np.float32)
    print(f"Contact map as float32: dtype = {contact_float.dtype}")

    # Normalized distance matrix
    dist_normalized = dist_matrix / dist_matrix.max()
    print(f"Normalized distances: range [{dist_normalized.min():.3f}, {dist_normalized.max():.3f}]")

    # Coordinates centered and scaled
    coords_centered = all_coords - np.mean(all_coords, axis=0)
    coords_scale = np.std(coords_centered)
    coords_normalized = coords_centered / coords_scale
    print(f"Normalized coordinates: mean = {np.mean(coords_normalized):.6f}, std = {np.std(coords_normalized):.3f}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
