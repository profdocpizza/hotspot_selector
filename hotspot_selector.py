#!/usr/bin/env python3
"""
hotspot_selector.py — Trim protein structures to surface hotspot regions.

Classifies residues as:
  Exposed   — surface residues with relative SASA above threshold
  Supporting — residues with any atom within `support_dist` Å of an exposed residue
  Other      — buried residues

Outputs:
  <stem>_annotated.pdb  — full structure; beta-factor encodes class (91/81/49)
  <stem>_hotspot.pdb    — Exposed + Supporting residues only
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

# Maximum accessible surface area (Å²) per residue type in extended Gly-X-Gly context
# Values from Tien et al. 2013 (empirical, all-atom)
MAX_SASA = {
    "ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0, "CYS": 148.0,
    "GLN": 214.0, "GLU": 214.0, "GLY": 97.0,  "HIS": 216.0, "ILE": 195.0,
    "LEU": 191.0, "LYS": 230.0, "MET": 203.0, "PHE": 228.0, "PRO": 154.0,
    "SER": 143.0, "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0,
}
DEFAULT_MAX_SASA = 200.0  # fallback for non-standard residues

BFACTOR_EXPOSED    = 91.0
BFACTOR_SUPPORTING = 81.0
BFACTOR_OTHER      = 49.0


def load_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def compute_sasa(structure):
    sr = ShrakeRupley()
    sr.compute(structure, level="R")  # stores sasa on each residue


def get_relative_sasa(residue) -> float:
    resname = residue.get_resname().strip()
    max_sasa = MAX_SASA.get(resname, DEFAULT_MAX_SASA)
    return residue.sasa / max_sasa


def classify_residues(structure, sasa_threshold: float, support_dist: float):
    """Return dicts mapping residue full_id -> class label."""
    all_residues = [
        res for model in structure
        for chain in model
        for res in chain
        if res.id[0] == " "  # skip HETATMs (waters, ligands)
    ]

    exposed_set = set()
    for res in all_residues:
        rel = get_relative_sasa(res)
        if rel >= sasa_threshold:
            exposed_set.add(res.full_id)

    # Build coordinate array for exposed atoms
    exposed_atoms = [
        atom for model in structure
        for chain in model
        for res in chain
        if res.full_id in exposed_set
        for atom in res.get_atoms()
    ]
    if exposed_atoms:
        exposed_coords = np.array([a.get_vector().get_array() for a in exposed_atoms])
    else:
        exposed_coords = np.empty((0, 3))

    classification = {}
    for res in all_residues:
        if res.full_id in exposed_set:
            classification[res.full_id] = "Exposed"
            continue

        # Check if any atom of this residue is within support_dist of any exposed atom
        res_coords = np.array([a.get_vector().get_array() for a in res.get_atoms()])
        if len(res_coords) > 0 and len(exposed_coords) > 0:
            diffs = res_coords[:, np.newaxis, :] - exposed_coords[np.newaxis, :, :]
            dists = np.sqrt((diffs ** 2).sum(axis=2))
            if dists.min() <= support_dist:
                classification[res.full_id] = "Supporting"
                continue

        classification[res.full_id] = "Other"

    return classification


def annotate_bfactors(structure, classification: dict):
    """Set b-factor on every atom according to residue classification."""
    for model in structure:
        for chain in model:
            for res in chain:
                label = classification.get(res.full_id, "Other")
                bf = {
                    "Exposed": BFACTOR_EXPOSED,
                    "Supporting": BFACTOR_SUPPORTING,
                }.get(label, BFACTOR_OTHER)
                for atom in res.get_atoms():
                    atom.set_bfactor(bf)


class HotspotSelect(Select):
    """Select only Exposed and Supporting residues."""
    def __init__(self, classification: dict):
        self.classification = classification

    def accept_residue(self, residue):
        label = self.classification.get(residue.full_id, "Other")
        return label in ("Exposed", "Supporting")


def save_pdb(structure, path: Path, select=None):
    io = PDBIO()
    io.set_structure(structure)
    if select is not None:
        io.save(str(path), select)
    else:
        io.save(str(path))


def print_summary(classification: dict, out_annotated: Path, out_hotspot: Path):
    counts = {"Exposed": 0, "Supporting": 0, "Other": 0}
    for label in classification.values():
        counts[label] += 1

    print("\n── Hotspot Selector Results ──────────────────────────────────")
    print(f"  Exposed    residues : {counts['Exposed']:>5}  (β = {BFACTOR_EXPOSED:.0f})")
    print(f"  Supporting residues : {counts['Supporting']:>5}  (β = {BFACTOR_SUPPORTING:.0f})")
    print(f"  Other      residues : {counts['Other']:>5}  (β = {BFACTOR_OTHER:.0f})")
    print(f"  Total               : {sum(counts.values()):>5}")
    print()
    print(f"  Annotated PDB  → {out_annotated}")
    print(f"  Hotspot PDB    → {out_hotspot}")
    print()
    print("── Visualisation ─────────────────────────────────────────────")
    print("  Open the annotated PDB in NanoViewer (https://nanoviewer.xyz)")
    print("  or Mol* (https://molstar.org/viewer/) and colour by pLDDT /")
    print("  B-factor to see Exposed (91), Supporting (81), Other (49).")
    print("──────────────────────────────────────────────────────────────\n")


def main():
    parser = argparse.ArgumentParser(
        description="Trim protein structures to surface hotspot regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Input structure file (.pdb or .cif)")
    parser.add_argument(
        "--sasa-threshold", type=float, default=0.25,
        help="Relative SASA threshold to classify a residue as Exposed (0–1)",
    )
    parser.add_argument(
        "--support-dist", type=float, default=5.0,
        help="Distance (Å) within which a residue is considered Supporting",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for output files (default: same as input)",
    )
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Error: input file not found: {args.input}")

    out_dir = args.output_dir if args.output_dir else args.input.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.input.stem
    out_annotated = out_dir / f"{stem}_annotated.pdb"
    out_hotspot   = out_dir / f"{stem}_hotspot.pdb"

    print(f"Loading structure: {args.input}")
    structure = load_structure(args.input)

    print("Computing SASA (ShrakeRupley)…")
    compute_sasa(structure)

    print(f"Classifying residues (SASA threshold={args.sasa_threshold}, support dist={args.support_dist} Å)…")
    classification = classify_residues(structure, args.sasa_threshold, args.support_dist)

    annotate_bfactors(structure, classification)
    save_pdb(structure, out_annotated)
    save_pdb(structure, out_hotspot, select=HotspotSelect(classification))

    print_summary(classification, out_annotated, out_hotspot)


if __name__ == "__main__":
    main()
