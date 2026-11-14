#!/usr/bin/env python3
"""Epitope: does the design hit the same patch on HER2 as trastuzumab does?

This is the second of the three axes the lab measures, and the only one we can
compute without relying on a confidence score. We do not look at how sure the
model is; we look at which HER2 residues actually sit next to the antibody in
the predicted structure, and compare the sets.

Two numbers per design, and the second is what makes the first readable:

  overlap_vs_tz_same_model   the design's epitope against trastuzumab's epitope,
                             predicted by the SAME model with the same seed
  tz_vs_crystal              trastuzumab's predicted epitope against the one in
                             the crystal structure (1N8Z)

The second number is the noise floor. If the model cannot reproduce
trastuzumab's own epitope, a deviation in a design means nothing at all.

    python scripts/epitope.py --runs runs --out results/epitope
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

CUTOFF = 4.5          # Å, heavy atom to heavy atom
AB_CHAINS = ("A", "B")        # light chain, heavy chain
AG_CHAINS = ("C", "D")        # HER2 domain IV, split at the unresolved loop

# Chains C and D are numbered 1.. in the predicted files, but are HER2 residues
# 489.. and 591.. in the crystal. Without this we compare two different numberings.
AG_OFFSET = {"C": 488, "D": 590}


def read_cif_atoms(path: Path) -> dict[str, list[tuple[int, np.ndarray]]]:
    """chain -> [(residue number, coordinate)] for all heavy atoms."""
    out: dict[str, list[tuple[int, np.ndarray]]] = {}
    in_atoms = False
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            if line.strip() in ("#", "") and out:
                break
            continue
        f = line.split()
        if len(f) < 13 or f[2] == "H":
            continue
        chain, resi = f[6], f[8]
        try:
            xyz = np.array([float(f[10]), float(f[11]), float(f[12])])
            resi_i = int(resi)
        except ValueError:
            continue
        out.setdefault(chain, []).append((resi_i, xyz))
    return out


def read_pdb_atoms(path: Path, chains) -> dict[str, list[tuple[int, np.ndarray]]]:
    out: dict[str, list[tuple[int, np.ndarray]]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        ch = line[21]
        if ch not in chains or line[76:78].strip() == "H":
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        out.setdefault(ch, []).append((int(line[22:26]), xyz))
    return out


def contact_residues(atoms, ab_chains, ag_chains, offset=None,
                     cutoff: float = CUTOFF) -> set[int]:
    """HER2 residues with a heavy atom closer than `cutoff` to the antibody."""
    offset = offset or {}
    ab = [xyz for c in ab_chains for _, xyz in atoms.get(c, [])]
    if not ab:
        return set()
    ab_arr = np.asarray(ab)
    hits: set[int] = set()
    for c in ag_chains:
        entries = atoms.get(c, [])
        if not entries:
            continue
        resis = np.array([r for r, _ in entries])
        coords = np.asarray([x for _, x in entries])
        # in blocks, so a 2500 x 1000 matrix is not built all at once
        for start in range(0, len(coords), 512):
            block = coords[start:start + 512]
            d = np.linalg.norm(block[:, None, :] - ab_arr[None, :, :], axis=-1)
            near = d.min(axis=1) <= cutoff
            for r in resis[start:start + 512][near]:
                hits.add(int(r) + offset.get(c, 0))
    return hits


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def find_runs(runs_dir: Path) -> dict[str, Path]:
    found = {}
    for run in sorted(runs_dir.glob("*-fold-*")):
        manifest, pred = run / "run.json", run / "predictions"
        if not manifest.exists() or not pred.is_dir():
            continue
        name = json.loads(manifest.read_text()).get("name", run.name)
        found[name.replace("fold-", "", 1)] = pred
    return found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, default=Path("runs"))
    ap.add_argument("--panel", type=Path, default=Path("runs/panel/panel.csv"))
    ap.add_argument("--crystal", type=Path, default=Path("data/1n8z.pdb"))
    ap.add_argument("--out", type=Path, default=Path("results/epitope"))
    ap.add_argument("--cutoff", type=float, default=CUTOFF)
    args = ap.parse_args()

    with args.panel.open(newline="") as fh:
        panel = {r["id"]: r for r in csv.DictReader(fh)}

    crystal: set[int] = set()
    if args.crystal.exists():
        # in 1N8Z, A is the light chain, B the heavy chain, C HER2 — and C is author-
        # numbered, i.e. already at 489.. , so no offset here
        crystal = contact_residues(read_pdb_atoms(args.crystal, ("A", "B", "C")),
                                   ("A", "B"), ("C",), cutoff=args.cutoff)

    rows = []
    for model, pred_dir in sorted(find_runs(args.runs).items()):
        epitopes: dict[str, set[int]] = {}
        for cif in sorted(pred_dir.glob("*/seed_*/*_model.cif")):
            qid = cif.parent.parent.name
            epitopes[qid] = contact_residues(read_cif_atoms(cif), AB_CHAINS,
                                             AG_CHAINS, AG_OFFSET, args.cutoff)
        tz = epitopes.get("trastuzumab")
        for qid, ep in sorted(epitopes.items()):
            rows.append({
                "model": model,
                "id": qid,
                "group": panel.get(qid, {}).get("group", "?"),
                "n_epitope_residues": len(ep),
                "overlap_vs_tz_same_model": (round(jaccard(ep, tz), 3)
                                             if tz is not None else ""),
                "overlap_vs_crystal": (round(jaccard(ep, crystal), 3)
                                       if crystal else ""),
            })

    if not rows:
        raise SystemExit(f"found no predicted structures under {args.runs}")

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "epitope.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"crystal epitope (1N8Z): {len(crystal)} HER2 residues")
    print(f"{len(rows)} rows -> {args.out/'epitope.csv'}\n")
    print(f'{"model":22s} {"tz vs crystal":>16s} {"design vs tz, median":>22s} '
          f'{"lowest":>9s}')
    for model in sorted({r["model"] for r in rows}):
        mine = [r for r in rows if r["model"] == model]
        tzrow = next((r for r in mine if r["id"] == "trastuzumab"), None)
        des = [r["overlap_vs_tz_same_model"] for r in mine
               if r["group"] == "design" and r["overlap_vs_tz_same_model"] != ""]
        print(f'{model:22s} {str(tzrow["overlap_vs_crystal"]) if tzrow else "-":>16s} '
              f'{(round(float(np.median(des)), 3) if des else "-")!s:>22s} '
              f'{(min(des) if des else "-")!s:>9s}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
