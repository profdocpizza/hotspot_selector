# Hotspot Selector

Trims protein structures to **surface hotspot regions**: exposed surface residues and the backbone-supporting residues around them.

## Residue classes

| Class | Description | β-factor |
|-------|-------------|----------|
| **Exposed** | Relative SASA ≥ threshold (default 0.10) — definitively on the solvent-accessible surface | **91** |
| **Supporting** | Any atom within 4 Å of an Exposed residue — backbone context for the hotspot | **81** |
| **Other** | Buried interior residues | **49** |

## Outputs

| File | Contents |
|------|----------|
| `<stem>_annotated.pdb` | Full structure; β-factor encodes class for visualisation |
| `<stem>_hotspot.pdb` | Exposed + Supporting residues only — pipeline-ready hotspot |
| `<stem>_hotspot_residue_index.txt` | Compact hotspot residue index string (for example `A96-116,A144-155,B1-7`) |

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
python hotspot_selector.py structure.pdb --anchor A:42 --surface-radius 30.0 --max_residues 400
```

### Example run (2XRP, anchors E:308 and A:108)

```bash
# Fetch 2XRP into examples/ (one-time setup)
mkdir -p examples outputs
curl -fsSL https://files.rcsb.org/download/2XRP.pdb -o examples/2XRP.pdb

# Run hotspot selection around anchor residues E:308 and A:108
python hotspot_selector.py examples/2XRP.pdb --anchor E:308 A:108 --output-dir ./outputs

# Inspect generated files (includes residue index .txt for easy copying)
ls -1 ./outputs/2XRP_*
cat ./outputs/2XRP_hotspot_residue_index.txt
```

### Anchor / surface-walk mode

Pass `--anchor CHAIN:RESNUM [...]` to restrict the output to a local surface patch around known residue(s).

Instead of Euclidean distance (which would select residues through the protein interior), the tool builds a **surface graph** — nodes are Exposed residues, edges connect Cα atoms within `--graph-step` Å — and runs **Dijkstra** from the anchor(s). Only exposed residues within `--surface-radius` Å of *surface-path* distance are selected, along with their Supporting neighbours.

**Cross-gap solvent edges** — controlled by `--graph-step`. The walk only considers edges between exposed residues whose Cα–Cα distance is within `--graph-step` (default 20 Å). Each candidate edge is validated: the Cα–Cα segment is sampled every ~1 Å and checked against a KD-tree of **buried (Other) residue atoms only** — if any sample point is within `--probe-radius` Å of a buried atom, the edge is rejected (it cuts through the protein core). Surface backbone atoms are correctly ignored. To allow hops across larger solvent gaps (cavities, clefts), simply increase `--graph-step`.

If `--max_residues` is set (anchor mode), the tool evaluates the selected hotspot size (Exposed + Supporting residues) and, when needed, automatically decreases the effective `--surface-radius` until the selection is within the cap.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--sasa-threshold` | `0.1` | Relative SASA cutoff for Exposed classification (0–1) |
| `--support-dist` | `4.0` | Distance (Å) used to identify Supporting residues |
| `--output-dir` | same dir as input | Where to write output files |
| `--anchor` | *(off)* | Activate surface-walk mode: one or more `CHAIN:RESNUM` anchors, e.g. `--anchor A:42 B:10` |
| `--surface-radius` | `25.0` | (Anchor mode) Max surface-path distance (Å) from anchor(s) |
| `--max_residues` | `None` | (Anchor mode) Cap hotspot size (Exposed+Supporting); if exceeded, `--surface-radius` is reduced automatically until the cap is met |
| `--graph-step` | `20.0` | (Anchor mode) Max Cα–Cα for a graph edge; increase to hop across larger solvent gaps |
| `--probe-radius` | `2.5` | (Anchor mode) Min clearance (Å) from any buried protein atom; edges cutting through buried core are rejected |

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
5. Writes the two PDB outputs and a hotspot residue index `.txt` file.
