# Hotspot Selector

Trims protein structures to **surface hotspot regions**: exposed surface residues and the backbone-supporting residues around them.

This tool is useful for:
- **Protein design** — focus computational effort on surface-exposed hotspots
- **Binding site analysis** — identify and visualize accessible surface patches
- **Interface design** — select regions for mutagenesis or optimization
- **Binder design** — scaffold generation around specific surface anchors

## Residue classes

Hotspot Selector classifies all residues into three categories based on solvent accessibility and spatial proximity:

| Class | Description | β-factor | Molecular Role |
|-------|-------------|----------|-----------------|
| **Exposed** | Relative SASA ≥ threshold (default 0.10) — definitively on the solvent-accessible surface | **91** | Directly accessible to solvent; primary targets for design or binding |
| **Supporting** | Any atom within 4 Å of an Exposed residue — backbone context for the hotspot | **81** | Provide structural context and support interactions; included to preserve backbone geometry |
| **Other** | Buried interior residues without Exposed neighbors | **49** | Buried protein core; typically unchanged in design workflows |

**SASA Background:** Relative SASA (Solvent-Accessible Surface Area) is computed using the **ShrakeRupley algorithm** against empirical maximum values per residue type (Tien et al. 2013). A residue with SASA ≥ 10% of its theoretical maximum is considered Exposed.

## Outputs

| File | Contents |
|------|----------|
| `<stem>_annotated.pdb` | Full structure; β-factor encodes class for visualisation |
| `<stem>_hotspot.pdb` | Exposed + Supporting residues only — pipeline-ready hotspot |
| `<stem>_hotspot_residue_index.txt` | Compact hotspot residue index string (for example `A96-116,A144-155,B1-7`) |

## Installation

**Requirements:**
- Python 3.8+
- Conda or Mamba (for environment management)

**Setup (first time only):**

```bash
# Clone and navigate to the repository
cd hotspot_selector

# Create conda environment with dependencies
conda env create -f environment.yml
conda activate hotspot
```

**Verify installation:**

```bash
python hotspot_selector.py --help
```

**Troubleshooting:**
- If `conda activate hotspot` doesn't work, try: `source activate hotspot` (for older conda versions)
- For environment conflicts, use `mamba` instead of `conda` (faster, more reliable)
- Note: The `examples/` directory is provided for testing only and is excluded from version control

## Usage

### Basic usage (minimal options)

Select the entire accessible surface:

```bash
# Load a PDB file (auto-detects format)
python hotspot_selector.py structure.pdb

# Or a CIF file
python hotspot_selector.py structure.cif

# Write outputs to a specific directory
python hotspot_selector.py structure.pdb --output-dir ./results
```

### Tuning sensitivity (global mode)

Adjust what counts as "exposed" or "supporting":

```bash
# Stricter Exposed threshold (10% → 20% SASA)
python hotspot_selector.py structure.pdb --sasa-threshold 0.20

# Larger Supporting shell (4 Å → 6 Å)
python hotspot_selector.py structure.pdb --support-dist 6.0

# Both together
python hotspot_selector.py structure.pdb --sasa-threshold 0.20 --support-dist 6.0
```

### Anchor / surface-walk mode (local patches)

Select a localized surface patch around one or more anchor residues:

```bash
# Single anchor
python hotspot_selector.py structure.pdb --anchor A:42

# Multiple anchors
python hotspot_selector.py structure.pdb --anchor A:42 A:43 B:10 --surface-radius 30.0

# With hotspot size cap (auto-reduces radius if needed)
python hotspot_selector.py structure.pdb --anchor A:42 --max_residues 200
```

### Example walkthrough (2XRP from PDB)

```bash
# Create directories (first time only)
mkdir -p examples outputs

# Download a test structure from RCSB PDB
curl -fsSL https://files.rcsb.org/download/2XRP.pdb -o examples/2XRP.pdb

# Run hotspot selection with a single anchor (B:409)
python hotspot_selector.py examples/2XRP.pdb --anchor B:409 --output-dir ./outputs

# List outputs
ls -lh ./outputs/2XRP_*

# View the residue index (compact format for pipelines)
cat ./outputs/2XRP_hotspot_residue_index.txt
```

**Expected output for B:409:**
```
Exposed residues: 64
Supporting residues: 147
Output files:
  - 2XRP_annotated.pdb (2.2M) — full structure with β-factors
  - 2XRP_hotspot.pdb (128K) — hotspot only
  - 2XRP_hotspot_residue_index.txt → A2-3,A102-104,...,B389-431
```

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

### How anchor / surface-walk mode works

**The Problem:** Euclidean distance selects residues through the protein interior, which is not useful for surface design. You need *surface-path* distance instead.

**The Solution:** Build a surface-only graph and use Dijkstra's algorithm to find true geodesic distance on the surface.

**Algorithm:**
1. **Build surface graph**: Nodes = all Exposed residues. Edges connect residues whose Cα atoms are within `--graph-step` Å (default 20 Å).
2. **Validate edges**: For each candidate edge, sample the Cα–Cα segment every ~1 Å and check against a KD-tree of buried (Other) atom positions. If any sample point is within `--probe-radius` Å of a buried atom, **reject the edge** (it cuts through the protein core). Backbone atoms of surface residues are correctly ignored.
3. **Run Dijkstra**: Find shortest *surface-path* distance from anchor residue(s) to all other Exposed residues.
4. **Select residues**: Keep all Exposed residues within `--surface-radius` Å of surface-path distance, plus their Supporting neighbors.
5. **Enforce size cap** (optional): If `--max_residues` is set, automatically reduce `--surface-radius` until the selection fits.

**Why it matters:** Surface pathways can avoid protein interior even when Euclidean distance is small. This is essential for binding site analysis and binder design, where you want realistic surface adjacency.

<p align="center">
  <img src="images/hotspot.png" alt="Hotspot surface patch visualization" width="500"/>
</p>

<p align="center">
  <em>Surface patch around an anchor residue. Exposed residues (blue), Supporting residues (cyan), and buried core (orange).</em>
</p>

## Visualisation

Open `<stem>_annotated.pdb` in **[NanoViewer](https://nanoviewer.xyz)** or **[Mol*](https://molstar.org/viewer/)** and colour by **pLDDT / B-factor** to see the classification:

### Color mapping table

| β-factor | Class | Color | Hex Code | Meaning |
|----------|-------|-------|----------|---------|
| **91** | **Exposed** | <span style="color: #0066ff">■</span> Blue | `#0066ff` | Directly accessible to solvent; primary design targets |
| **81** | **Supporting** | <span style="color: #10cff1">■</span> Cyan | `#10cff1` | Backbone context; structural support for hotspot |
| **49** | **Other** | <span style="color: #ff8c00">■</span> Orange | `#ff8c00` | Buried protein interior; typically preserved |

**Visualization tips:**
- The annotated structure shows the full protein with residue classes encoded as β-factors
- Use the color scale to quickly identify which residues will be retained in the hotspot
- The hotspot output (`*_hotspot.pdb`) contains only Exposed + Supporting residues

<p align="center">
  <img src="images/hotpsot_annotated.png" alt="Annotated structure with color-coded residue classes" width="500"/>
</p>

<p align="center">
  <em>Annotated PDB structure colored by residue class (β-factor). Blue = Exposed (design targets), Cyan = Supporting (backbone), Orange = Other (buried core).</em>
</p>

## How it works

The pipeline consists of four main steps:

### Step 1: Load structure
Loads the protein structure using **BioPython**, automatically detecting PDB or CIF format.

### Step 2: Compute SASA
Runs the **ShrakeRupley solvent-accessible surface area (SASA)** algorithm, built into BioPython.
- No external dependencies required
- Computes atomic SASA, then aggregates to per-residue values
- Fast (~1–2 seconds for typical structures)

### Step 3: Classify residues
Computes **relative SASA** per residue by dividing by empirical maximum values (from Tien et al. 2013):
- Residues with relative SASA ≥ `--sasa-threshold` (default 0.10) → **Exposed**
- All other residues → queued for Supporting check

### Step 4: Find Supporting residues
For each Exposed residue, queries all atoms within `--support-dist` Å (default 4 Å):
- Any residue with an atom in this radius → **Supporting**
- Provides structural context and backbone geometry
- Preserves local geometry for design workflows

### Step 5: Write outputs
Generates three output files:
- `*_annotated.pdb` — Full structure with β-factors encoding residue class (91/81/49)
- `*_hotspot.pdb` — Hotspot only (Exposed + Supporting residues)
- `*_hotspot_residue_index.txt` — Compact residue range string (e.g., `A2-4,A53-55,...,B406-412`)
