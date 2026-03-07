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

# Anchor / surface-walk mode: select a surface patch around known residues
python hotspot_selector.py structure.pdb --anchor A:42
python hotspot_selector.py structure.pdb --anchor A:42 A:43 B:10 --surface-radius 30.0
```

### Anchor / surface-walk mode

Pass `--anchor CHAIN:RESNUM [...]` to restrict the output to a local surface patch around known residue(s).

Instead of Euclidean distance (which would select residues through the protein interior), the tool builds a **surface graph** — nodes are Exposed residues, edges connect Cα atoms within `--graph-step` Å — and runs **Dijkstra** from the anchor(s). Only exposed residues within `--surface-radius` Å of *surface-path* distance are selected, along with their Supporting neighbours.

**Cross-gap solvent edges** (`--cross-gap`, default 20 Å): allows the walk to jump across solvent-filled cavities and clefts — useful for curved, concave, or hollow surfaces where the around-the-surface path is much longer than the straight-line distance. Each candidate long edge is validated by sampling the Cα–Cα segment and rejecting it if any sample point is within `--probe-radius` Å of a protein atom (i.e., it passes through the interior). Set `--cross-gap` equal to `--graph-step` to disable.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--sasa-threshold` | `0.25` | Relative SASA cutoff for Exposed classification (0–1) |
| `--support-dist` | `5.0` | Distance (Å) used to identify Supporting residues |
| `--output-dir` | same dir as input | Where to write output files |
| `--anchor` | *(off)* | Activate surface-walk mode: one or more `CHAIN:RESNUM` anchors, e.g. `--anchor A:42 B:10` |
| `--surface-radius` | `25.0` | (Anchor mode) Max surface-path distance (Å) from anchor(s) |
| `--graph-step` | `10.0` | (Anchor mode) Max Cα–Cα distance (Å) for a surface graph edge |
| `--cross-gap` | `20.0` | (Anchor mode) Max Cα–Cα for solvent-crossing edges; set equal to `--graph-step` to disable |
| `--probe-radius` | `2.5` | (Anchor mode) Min clearance (Å) from any protein atom for a cross-gap edge |

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
