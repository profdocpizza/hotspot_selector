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

Anchor / surface-walk mode (--anchor):
  Restricts output to exposed residues reachable from anchor residue(s) by walking
  *along the surface* (Dijkstra on a Cα graph of Exposed residues only), plus their
  Supporting neighbours.  Paths cannot shortcut through buried interior.
"""

import argparse
import heapq
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


def build_surface_graph(exposed_residues: list, graph_step: float) -> dict:
    """
    Build an adjacency graph over Exposed residues using Cα–Cα distances.

    Nodes  : residue full_id
    Edges  : pair of exposed residues whose Cα–Cα distance ≤ graph_step
    Returns: dict { full_id -> list of (neighbour_full_id, distance) }
    """
    # Collect Cα positions; fall back to first atom if no Cα present
    nodes = []
    for res in exposed_residues:
        if "CA" in res:
            coord = res["CA"].get_vector().get_array()
        else:
            coord = next(iter(res.get_atoms())).get_vector().get_array()
        nodes.append((res.full_id, coord))

    adj = {fid: [] for fid, _ in nodes}
    coords = np.array([c for _, c in nodes])
    ids = [fid for fid, _ in nodes]

    # Pairwise distances (O(n²) but surface residues are a fraction of the total)
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dists = np.sqrt((diff ** 2).sum(axis=2))

    n = len(ids)
    for i in range(n):
        for j in range(i + 1, n):
            d = dists[i, j]
            if d <= graph_step:
                adj[ids[i]].append((ids[j], float(d)))
                adj[ids[j]].append((ids[i], float(d)))

    return adj


def surface_walk(adj: dict, anchor_ids: set, surface_radius: float) -> set:
    """
    Dijkstra from the set of anchor residues through the surface graph.

    Returns the set of residue full_ids reachable within `surface_radius` Å
    of accumulated surface-path distance.
    """
    # Multi-source Dijkstra: initialise heap with all anchors at distance 0
    dist = {fid: float("inf") for fid in adj}
    heap = []
    for aid in anchor_ids:
        if aid in dist:
            dist[aid] = 0.0
            heapq.heappush(heap, (0.0, aid))

    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u]:
            continue
        for v, w in adj.get(u, []):
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                heapq.heappush(heap, (nd, v))

    return {fid for fid, d in dist.items() if d <= surface_radius}


def parse_anchor(token: str) -> tuple:
    """Parse 'CHAIN:RESNUM' or 'CHAIN:RESNUM:ICODE' into (chain, resnum, icode)."""
    parts = token.split(":")
    if len(parts) < 2:
        raise argparse.ArgumentTypeError(
            f"Anchor '{token}' must be in CHAIN:RESNUM format, e.g. A:42"
        )
    chain = parts[0]
    try:
        resnum = int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Residue number in '{token}' must be an integer"
        )
    icode = parts[2] if len(parts) > 2 else " "
    return (chain, resnum, icode)


def resolve_anchor_full_ids(structure, anchors: list) -> set:
    """Map parsed anchor tuples to BioPython full_ids, with helpful errors."""
    resolved = set()
    for chain_id, resnum, icode in anchors:
        found = False
        for model in structure:
            if chain_id not in model:
                continue
            chain = model[chain_id]
            res_id = (" ", resnum, icode)
            if res_id in chain:
                resolved.add(chain[res_id].full_id)
                found = True
                break
        if not found:
            sys.exit(
                f"Error: anchor residue {chain_id}:{resnum} not found in structure."
            )
    return resolved


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


def print_summary(classification: dict, out_annotated: Path, out_hotspot: Path, anchor_mode: bool = False):
    counts = {"Exposed": 0, "Supporting": 0, "Other": 0}
    for label in classification.values():
        counts[label] += 1

    mode_tag = " [anchor/surface-walk]" if anchor_mode else ""
    print(f"\n── Hotspot Selector Results{mode_tag} ──────────────────────────────")
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
    # Anchor / surface-walk mode
    parser.add_argument(
        "--anchor", nargs="+", metavar="CHAIN:RESNUM", default=None,
        help=(
            "Activate surface-walk mode: one or more anchor residues as CHAIN:RESNUM "
            "(e.g. --anchor A:42 A:43). Only Exposed residues reachable from the "
            "anchor(s) along the surface (within --surface-radius) are kept."
        ),
    )
    parser.add_argument(
        "--surface-radius", type=float, default=25.0,
        help="(Anchor mode) Max surface-path distance (Å) from anchor(s)",
    )
    parser.add_argument(
        "--graph-step", type=float, default=10.0,
        help="(Anchor mode) Max Cα–Cα distance (Å) for a surface graph edge",
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

    anchor_mode = args.anchor is not None
    if anchor_mode:
        anchors = [parse_anchor(t) for t in args.anchor]
        anchor_full_ids = resolve_anchor_full_ids(structure, anchors)

        # Verify anchors are Exposed
        for fid in anchor_full_ids:
            if classification.get(fid) != "Exposed":
                label = classification.get(fid, "unknown / HETATM")
                print(
                    f"  Warning: anchor {fid} is classified as '{label}', not Exposed. "
                    "It will still be used as a surface-walk start point."
                )

        # Build surface graph over all Exposed residues
        exposed_residues = [
            res for model in structure
            for chain in model
            for res in chain
            if classification.get(res.full_id) == "Exposed"
        ]
        print(
            f"Building surface graph ({len(exposed_residues)} exposed residues, "
            f"edge step ≤ {args.graph_step} Å)…"
        )
        adj = build_surface_graph(exposed_residues, args.graph_step)

        print(f"Surface walk from {len(anchor_full_ids)} anchor(s), radius ≤ {args.surface_radius} Å…")
        reachable_exposed = surface_walk(adj, anchor_full_ids, args.surface_radius)

        # Rebuild classification: keep only reachable exposed + their supporting neighbours
        # For Supporting, recompute from scratch restricted to reachable exposed atoms
        reachable_exposed_atoms = [
            atom for model in structure
            for chain in model
            for res in chain
            if res.full_id in reachable_exposed
            for atom in res.get_atoms()
        ]
        if reachable_exposed_atoms:
            reachable_coords = np.array(
                [a.get_vector().get_array() for a in reachable_exposed_atoms]
            )
        else:
            reachable_coords = np.empty((0, 3))

        restricted = {}
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    if res.full_id in reachable_exposed:
                        restricted[res.full_id] = "Exposed"
                        continue
                    # Supporting: any atom within support_dist of a reachable exposed atom
                    res_coords = np.array(
                        [a.get_vector().get_array() for a in res.get_atoms()]
                    )
                    if len(res_coords) > 0 and len(reachable_coords) > 0:
                        diffs = res_coords[:, np.newaxis, :] - reachable_coords[np.newaxis, :, :]
                        dists = np.sqrt((diffs ** 2).sum(axis=2))
                        if dists.min() <= args.support_dist:
                            restricted[res.full_id] = "Supporting"
                            continue
                    restricted[res.full_id] = "Other"

        classification = restricted
        print(
            f"  Selected {sum(1 for v in classification.values() if v=='Exposed')} exposed "
            f"and {sum(1 for v in classification.values() if v=='Supporting')} supporting "
            "residues in anchor neighbourhood."
        )

    annotate_bfactors(structure, classification)
    save_pdb(structure, out_annotated)
    save_pdb(structure, out_hotspot, select=HotspotSelect(classification))

    print_summary(classification, out_annotated, out_hotspot, anchor_mode=anchor_mode)


if __name__ == "__main__":
    main()
