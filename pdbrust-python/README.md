# PDBRust Python Bindings

High-performance Python bindings for [PDBRust](https://github.com/HFooladi/pdbrust), a Rust library for parsing and analyzing PDB/mmCIF protein structure files.

## Installation

```bash
pip install pdbrust
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

### RCSB PDB Integration

```python
from pdbrust import SearchQuery, rcsb_search, download_structure, FileFormat

# Download a structure
structure = download_structure("1UBQ", FileFormat.pdb())

# Search RCSB
query = SearchQuery().with_text("kinase").with_organism("Homo sapiens").with_resolution_max(2.0)
results = rcsb_search(query, 10)
print(f"Found {results.total_count} structures")
for pdb_id in results.pdb_ids:
    print(f"  {pdb_id}")
```

## Performance

PDBRust provides **40-260x speedups** over pure Python implementations:

| Operation | Speedup vs Python |
|-----------|-------------------|
| Parsing | 2-3x |
| get_ca_coords | 240x |
| max_ca_distance | 260x |
| radius_of_gyration | 100x |

## Requirements

- Python >= 3.9
- No runtime dependencies (Rust code is compiled into the package)

## License

MIT
