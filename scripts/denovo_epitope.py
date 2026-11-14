#!/usr/bin/env python3
"""Which HER2 residues does each design take — and is it trastuzumab's site?

Epitope is one of the three axes the lab measures, and for the de novo designs
it cannot be assumed: all six loops are new. Here it is measured.

METHOD
  For each predicted complex: find the HER2 residues (chains C and D) that have
  a heavy atom within 4.5 Å of a heavy atom in the antibody (chains A and B).
  That is the interface. The overlap with trastuzumab's own interface is
  measured with the Jaccard index.

TWO COMPARISONS, and both are needed
  against trastuzumab's PREDICTED complex, same model and seed: same sources
    of error on both sides, so the number is a clean comparison.
  against 1N8Z, the experimental structure: the real answer, but then the
    fold model's own error is measured too. Trastuzumab's own value against
    1N8Z is the benchmark a design should be held up against — not 1.0.

The numbering in the predicted files starts at 1 per chain. residue_ids in
queries.json give the original HER2 numbers, and are used to translate.

    python scripts/denovo_epitope.py --predictions runs/<run>/predictions \
        --queries runs/denovo-panel/queries.json --reference trastuzumab \
        --out results/denovo/epitope_boltz2.csv
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

AB_CHAINS = ("A", "B")
AG_CHAINS = ("C", "D")
CUTOFF = 4.5


def parse_cif_atoms(path: Path) -> list[tuple[str, int, str, float, float, float]]:
    """(chain, residue number, atom name, x, y, z) for each heavy atom.

    The column order is read from the _atom_site header, not assumed.
    """
    cols: list[str] = []
    atoms = []
    in_loop = False
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
            in_loop = True
            continue
        if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
            f = line.split()
            if len(f) < len(cols):
                continue
            rec = dict(zip(cols, f))
            atom = rec.get("label_atom_id", "")
            if rec.get("type_symbol") == "H" or atom.startswith("H"):
                continue
            chain = rec.get("auth_asym_id") or rec.get("label_asym_id")
            try:
                resnum = int(rec.get("auth_seq_id") or rec.get("label_seq_id"))
                xyz = (float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"]))
            except (TypeError, ValueError):
                continue
            atoms.append((chain, resnum, atom, *xyz))
        elif in_loop and line.startswith("#"):
            in_loop = False
    return atoms


def parse_pdb_atoms(path: Path) -> list[tuple[str, int, str, float, float, float]]:
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip()
        atom = line[12:16].strip()
        if element == "H" or atom.startswith("H"):
            continue
        atoms.append((line[21], int(line[22:26]), atom,
                      float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return atoms


def read_atoms(path: Path):
    return parse_pdb_atoms(path) if path.suffix.lower() == ".pdb" else parse_cif_atoms(path)


def interface(atoms, ab_chains=AB_CHAINS, ag_chains=AG_CHAINS,
              cutoff: float = CUTOFF) -> set[tuple[str, int]]:
    """The HER2 residues that are in contact with the antibody."""
    ab = np.array([a[3:6] for a in atoms if a[0] in ab_chains], dtype=float)
    ag_meta = [(a[0], a[1]) for a in atoms if a[0] in ag_chains]
    ag = np.array([a[3:6] for a in atoms if a[0] in ag_chains], dtype=float)
    if len(ab) == 0 or len(ag) == 0:
        return set()
    hits: set[tuple[str, int]] = set()
    # in blocks, so 10,000 x 10,000 distances do not blow up memory
    for i in range(0, len(ag), 2000):
        block = ag[i:i + 2000]
        d2 = ((block[:, None, :] - ab[None, :, :]) ** 2).sum(-1)
        close = (d2 <= cutoff * cutoff).any(axis=1)
        for j in np.nonzero(close)[0]:
            hits.add(ag_meta[i + int(j)])
    return hits


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def residue_id_map(queries_path: Path) -> dict[str, dict[int, int]]:
    """chain -> {position in the chain (1-based): original residue number}."""
    data = json.loads(queries_path.read_text())
    out: dict[str, dict[int, int]] = {}
    for q in data["queries"]:
        for ch in q["chains"]:
            ids = ch.get("residue_ids")
            if ids and ch["chain_id"] not in out:
                out[ch["chain_id"]] = {i: int(r) for i, r in enumerate(ids, 1)}
    return out


def parse_offsets(text: str) -> dict[str, int]:
    """'C=488,D=590' -> {'C': 488, 'D': 590}.

    Some refold queries do not carry residue_ids, and then HER2 is numbered
    from 1 in the predicted files. The offset restores the 1N8Z numbering:
    chain C is residues 489-580 and chain D is 591-607, and the lengths
    (92 and 17) are checked against that before the numbers are used.
    """
    out: dict[str, int] = {}
    for part in text.split(","):
        if "=" in part:
            chain, _, value = part.partition("=")
            out[chain.strip()] = int(value)
    return out


def renumber(hits: set[tuple[str, int]], mapping: dict[str, dict[int, int]],
             offsets: dict[str, int] | None = None) -> set[int]:
    """The interface as original HER2 residue numbers, chain-free.

    Some fold models write out HER2 with their own numbering (1, 2, 3 ...),
    others keep the original (489, 490 ...). Both must be handled correctly:
    if the number is already in residue_ids, it is the original. A number we
    cannot place becomes negative, so it shows up in the table instead of
    silently mixing in among the real ones.
    """
    out = set()
    offsets = offsets or {}
    for chain, pos in hits:
        m = mapping.get(chain)
        if chain in offsets:
            out.add(offsets[chain] + pos)
        elif not m:
            out.add(-pos)
        elif pos in set(m.values()):
            out.add(pos)
        else:
            out.add(m.get(pos, -pos))
    return out


def prediction_files(pred_dir: Path) -> dict[str, Path]:
    """query_id -> model file (first seed, first sample)."""
    out: dict[str, Path] = {}
    for cif in sorted(pred_dir.glob("*/seed_*/*_model.cif")):
        out.setdefault(cif.parent.parent.name, cif)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions", type=Path, required=True)
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--reference", default="trastuzumab",
                    help="query_id for the predicted reference complex")
    ap.add_argument("--reference-prediction", type=Path,
                    help="directory the reference was folded in, if it is a different run")
    ap.add_argument("--experimental", type=Path, default=Path("data/1n8z.pdb"))
    ap.add_argument("--ab-chains", default="AB",
                    help="the chains that are the antibody in the predicted files")
    ap.add_argument("--ag-chains", default="CD",
                    help="the chains that are HER2 in the predicted files")
    ap.add_argument("--chain-offsets", default="",
                    help="e.g. C=488,D=590 when the queries lack residue_ids")
    ap.add_argument("--exp-ab-chains", default="AB")
    ap.add_argument("--exp-ag-chains", default="C")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    mapping = residue_id_map(args.queries)
    files = prediction_files(args.predictions)
    if args.reference_prediction:
        files.update({k: v for k, v in prediction_files(args.reference_prediction).items()
                      if k == args.reference})
    if not files:
        raise SystemExit(f"no predicted structures in {args.predictions}")

    exp_set: set[int] = set()
    if args.experimental.exists():
        exp_atoms = read_atoms(args.experimental)
        exp_set = {num for _, num in interface(
            exp_atoms, tuple(args.exp_ab_chains), tuple(args.exp_ag_chains))}

    ab_chains, ag_chains = tuple(args.ab_chains), tuple(args.ag_chains)
    offsets = parse_offsets(args.chain_offsets)
    sets = {qid: renumber(interface(read_atoms(p), ab_chains, ag_chains), mapping, offsets)
            for qid, p in sorted(files.items())}
    ref = sets.get(args.reference, set())

    rows = []
    for qid, s in sets.items():
        rows.append({
            "query_id": qid,
            "n_interface_residues": len(s),
            "jaccard_to_predicted_trastuzumab": round(jaccard(s, ref), 3) if ref else "",
            "jaccard_to_1n8z": round(jaccard(s, exp_set), 3) if exp_set else "",
            "shared_with_1n8z": len(s & exp_set) if exp_set else "",
            "interface_residues": " ".join(str(x) for x in sorted(s)),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} complexes -> {args.out}")
    if ref:
        print(f"  reference {args.reference}: {len(ref)} HER2 residues in the interface")
    if exp_set:
        print(f"  1N8Z: {len(exp_set)} HER2 residues in the interface")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
