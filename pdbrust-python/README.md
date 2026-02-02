# PDBRust Python Bindings

High-performance Python bindings for [PDBRust](https://github.com/HFooladi/pdbrust), a Rust library for parsing and analyzing PDB/mmCIF protein structure files.

## Installation

```bash
pip install pdbrust
```

## Development Installation

To build and install from source (useful for development or testing latest changes):

### Prerequisites
- Python 3.9+
- Rust toolchain (1.85.0+)
- [uv](https://github.com/astral-sh/uv) (fast Python package manager)
- [maturin](https://github.com/PyO3/maturin) (Rust-Python build tool)

### Setup

1. Clone the repository and navigate to the Python bindings:
   ```bash
   cd pdbrust-python
   ```

2. Create a virtual environment with uv:
   ```bash
   uv venv
   ```

3. Activate the virtual environment:
   ```bash
   # Linux/macOS
   source .venv/bin/activate

   # Windows
   .venv\Scripts\activate
   ```

4. Install maturin (if not already installed):
   ```bash
   uv pip install maturin
   ```

5. Build and install pdbrust in development mode:
   ```bash
   maturin develop --release
   ```

6. (Optional) Install numpy for array support:
   ```bash
   uv pip install numpy
   ```

## Quick Start

```python
import pdbrust

# Parse a PDB file
structure = pdbrust.parse_pdb_file("protein.pdb")
print(f"Loaded {structure.num_atoms} atoms in {structure.num_chains} chains")

# Get chain IDs
chains = structure.get_chain_ids()
print(f"Chains: {chains}")

# Access atoms
for atom in structure.atoms[:5]:
    print(f"{atom.name} {atom.residue_name}{atom.residue_seq}")
```

## Features

### Parsing

```python
# Parse different formats
structure = pdbrust.parse_pdb_file("protein.pdb")
structure = pdbrust.parse_mmcif_file("protein.cif")
structure = pdbrust.parse_structure_file("protein.ent")  # auto-detect

# Parse gzip-compressed files
structure = pdbrust.parse_gzip_pdb_file("pdb1ubq.ent.gz")

# Parse from string
structure = pdbrust.parse_pdb_string(pdb_content)
```

### Filtering and Cleaning

```python
# Method chaining for clean code
cleaned = structure.remove_ligands().keep_only_chain("A").keep_only_ca()

# Get CA coordinates
ca_coords = structure.get_ca_coords()  # List of (x, y, z) tuples
ca_coords_chain_a = structure.get_ca_coords("A")  # Specific chain

# Cleaning operations
structure.center_structure()
structure.normalize_chain_ids()
structure.reindex_residues()
```

### Structural Descriptors

```python
# Individual metrics
rg = structure.radius_of_gyration()
max_dist = structure.max_ca_distance()
composition = structure.aa_composition()

# All descriptors at once
desc = structure.structure_descriptors()
print(f"Rg: {desc.radius_of_gyration:.2f} A")
print(f"Hydrophobic: {desc.hydrophobic_ratio:.1%}")
```

### Quality Assessment

```python
# Quick checks
if structure.has_altlocs():
    print("Warning: alternate conformations present")

if structure.has_multiple_models():
    print("NMR ensemble detected")

# Full quality report
report = structure.quality_report()
if report.is_analysis_ready():
    print("Structure is ready for analysis")
```

### Writing Files

```python
import pdbrust

structure = pdbrust.parse_pdb_file("input.pdb")

# Write to PDB format
pdbrust.write_pdb_file(structure, "output.pdb")

# Write to mmCIF format
pdbrust.write_mmcif_file(structure, "output.cif")

# Write compressed mmCIF
pdbrust.write_gzip_mmcif_file(structure, "output.cif.gz")

# Get mmCIF as string (useful for web APIs, in-memory processing)
mmcif_string = pdbrust.write_mmcif_string(structure)
```

### Geometry: RMSD and Structure Alignment

```python
from pdbrust import AtomSelection

# Load two structures to compare
structure1 = pdbrust.parse_pdb_file("structure1.pdb")
structure2 = pdbrust.parse_pdb_file("structure2.pdb")

# Calculate RMSD (without alignment)
rmsd = structure1.rmsd_to(structure2)
print(f"RMSD: {rmsd:.3f} Å")

# RMSD with different atom selections
rmsd_ca = structure1.rmsd_to(structure2, AtomSelection.ca_only())      # CA atoms (default)
rmsd_bb = structure1.rmsd_to(structure2, AtomSelection.backbone())    # Backbone (N, CA, C, O)
rmsd_all = structure1.rmsd_to(structure2, AtomSelection.all_atoms())  # All atoms

# Align structures (Kabsch algorithm) - returns aligned structure and result
aligned, result = structure1.align_to(structure2)
print(f"Alignment RMSD: {result.rmsd:.3f} Å ({result.num_atoms} atoms)")

# Per-residue RMSD for flexibility analysis
per_res = structure1.per_residue_rmsd_to(structure2)
for r in per_res:
    if r.rmsd > 2.0:  # Highlight flexible regions
        print(f"Flexible: {r.chain_id}{r.residue_seq} {r.residue_name}: {r.rmsd:.2f} Å")
```

### Numpy Integration

```python
import numpy as np

structure = pdbrust.parse_pdb_file("protein.pdb")

# Get coordinates as numpy arrays
all_coords = structure.get_coords_array()          # Shape: (N_atoms, 3)
ca_coords = structure.get_ca_coords_array()        # Shape: (N_ca, 3)
bb_coords = structure.get_backbone_coords_array()  # Shape: (N_backbone, 3)

# Chain-specific coordinates
chain_a_ca = structure.get_ca_coords_array("A")

# Distance matrix (pairwise CA-CA distances)
dist_matrix = structure.distance_matrix_ca()  # Shape: (N_ca, N_ca)

# Contact map (binary matrix of contacts within threshold)
contact_map = structure.contact_map_ca(threshold=8.0)  # Default: 8 Å

# All-atom versions
all_dist = structure.distance_matrix()
all_contacts = structure.contact_map(threshold=4.5)

# Use with machine learning
print(f"Coords shape: {all_coords.shape}")
print(f"Contact map shape: {contact_map.shape}, contacts: {contact_map.sum()}")
```

### RCSB PDB Integration

```python
from pdbrust import SearchQuery, rcsb_search, download_structure, FileFormat

# Download a structure
structure = download_structure("1UBQ", FileFormat.pdb())

# Download to file directly
pdbrust.download_to_file("1UBQ", "1ubq.pdb", FileFormat.pdb())

# Get as string without saving
pdb_string = pdbrust.download_pdb_string("1UBQ", FileFormat.pdb())

# Search RCSB with various filters
query = (SearchQuery()
    .with_text("kinase")
    .with_organism("Homo sapiens")
    .with_resolution_max(2.0)
    .with_experimental_method(ExperimentalMethod.xray())
    .with_sequence_length_min(100)
    .with_sequence_length_max(500))

results = rcsb_search(query, 10)
print(f"Found {results.total_count} structures")
for pdb_id in results.pdb_ids:
    print(f"  {pdb_id}")
```

### Ligand Pose Quality (PoseBusters-style Checks)

```python
structure = pdbrust.parse_pdb_file("protein_ligand.pdb")

# List all ligands in the structure
ligands = structure.get_ligand_names()
print(f"Ligands: {ligands}")

# Validate a specific ligand
report = structure.ligand_pose_quality("LIG")
if report:
    print(f"Ligand: {report.ligand_name}")
    print(f"Min distance: {report.min_protein_ligand_distance:.2f} Å")
    print(f"Clashes: {report.num_clashes}")
    print(f"Volume overlap: {report.protein_volume_overlap_pct:.1f}%")

    if report.is_geometry_valid:
        print("✓ Pose passes geometry checks")
    else:
        print("✗ Pose fails geometry checks")
        for clash in report.clashes[:3]:
            print(f"  Clash: {clash.protein_residue_name} {clash.protein_atom_name} - "
                  f"{clash.ligand_atom_name}: {clash.distance:.2f}Å")

# Validate all ligands
for report in structure.all_ligand_pose_quality():
    status = "PASS" if report.is_geometry_valid else "FAIL"
    print(f"{report.ligand_name}: {status}")
```

### Additional Structure Methods

```python
# Access sequence from SEQRES records
sequence = structure.get_sequence("A")

# Get residues for a specific chain
residues = structure.get_residues_for_chain("A")  # List of (seq_num, name) tuples

# Access connectivity (CONECT records)
connected = structure.get_connected_atoms(atom_serial=1)

# Get center of mass
centroid = structure.get_centroid()
ca_centroid = structure.get_ca_centroid()

# Translate structure
structure.translate(10.0, 0.0, 0.0)
```

## Running Examples

The `examples/` directory contains Python scripts demonstrating various features.

### Setup
1. Navigate to the pdbrust-python directory and activate your virtual environment:
   ```bash
   cd pdbrust-python
   source .venv/bin/activate  # Linux/macOS
   ```

2. Navigate to the examples directory:
   ```bash
   cd examples
   ```

3. Run any example:
   ```bash
   python basic_usage.py
   python geometry_rmsd.py
   python numpy_integration.py
   ```

### Available Examples
| Example | Description |
|---------|-------------|
| `basic_usage.py` | Parsing, accessing atoms/residues, basic filtering |
| `writing_files.py` | Write PDB/mmCIF files |
| `geometry_rmsd.py` | RMSD calculation, structure alignment |
| `lddt_demo.py` | LDDT calculation (superposition-free) |
| `numpy_integration.py` | Coordinate arrays, distance matrices |
| `rcsb_search.py` | RCSB search queries and downloads |
| `selection_language.py` | PyMOL/VMD-style selection language |
| `secondary_structure.py` | DSSP secondary structure assignment |
| `b_factor_analysis.py` | B-factor statistics and analysis |
| `alphafold_analysis.py` | AlphaFold pLDDT confidence analysis |
| `quality_and_summary.py` | Quality reports, structure summaries |
| `batch_processing.py` | Process multiple files |
| `advanced_filtering.py` | Method chaining, normalization |

> **Note:** Some examples require sample PDB files. You can download test structures from RCSB or use the files in `../examples/pdb_files/`.

## Performance

PDBRust provides **40-260x speedups** over pure Python implementations:

| Operation | Speedup vs Python |
|-----------|-------------------|
| Parsing | 2-3x |
| get_ca_coords | 240x |
| max_ca_distance | 260x |
| radius_of_gyration | 100x |

## Requirements

- Python 3.9-3.13
- No runtime dependencies (Rust code is compiled into the package)

## License

MIT
