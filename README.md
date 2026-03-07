# Hotspot Selector

Trims protein structures to **surface hotspot regions**: exposed surface residues and the backbone-supporting residues around them.

## Residue classes

| Class | Description | β-factor |
|-------|-------------|----------|
| **Exposed** | Relative SASA ≥ threshold (default 0.25) — definitively on the solvent-accessible surface | **91** |
| **Supporting** | Any atom within 5 Å of an Exposed residue — backbone context for the hotspot | **81** |
| **Other** | Buried interior residues | **49** |

## Outputs

| File | Contents |
|------|----------|
| `<stem>_annotated.pdb` | Full structure; β-factor encodes class for visualisation |
| `<stem>_hotspot.pdb` | Exposed + Supporting residues only — pipeline-ready hotspot |

## Installation

```bash
conda env create -f environment.yml
conda activate hotspot
```

## Usage

```bash
# Basic usage (PDB or CIF)
python hotspot_selector.py structure.pdb
python hotspot_selector.py structure.cif

# Tune thresholds
python hotspot_selector.py structure.pdb --sasa-threshold 0.20 --support-dist 6.0

# Write outputs to a specific directory
python hotspot_selector.py structure.pdb --output-dir ./results
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--sasa-threshold` | `0.25` | Relative SASA cutoff for Exposed classification (0–1) |
| `--support-dist` | `5.0` | Distance (Å) used to identify Supporting residues |
| `--output-dir` | same dir as input | Where to write output files |

## Visualisation

Open `<stem>_annotated.pdb` in **[NanoViewer](https://nanoviewer.xyz)** or **[Mol*](https://molstar.org/viewer/)** and colour by **pLDDT / B-factor**:

- **Blue / 91** → Exposed surface residues
- **Cyan / 81** → Supporting residues  
- **Orange / 49** → Other (buried) residues

This gives an immediate visual overview of which parts of the structure will be retained in the hotspot output.

## How it works

1. Loads the structure with BioPython (auto-detects PDB vs CIF format).
2. Runs **ShrakeRupley SASA** calculation (BioPython built-in, no external dependencies).
3. Computes relative SASA per residue against empirical maximum values (Tien et al. 2013).
4. Labels Exposed residues → builds a neighbour search → labels Supporting residues.
5. Writes the two output files.
