#!/usr/bin/env python3
"""Build queries.json for the whole panel, in the same form that refold-prepare produces.

The basis is a real `prepared/holo/queries.json` from a refold-prepare run.
From it, chain A (VL), C and D (HER2 domain IV, split at the unresolved
loop) are inherited with their residue_ids. Only chain B (VH) is swapped
out, once per panel entry. The template (template.cif) omits the ten CDR-H3
residues, so the same template applies to all the variants.

    python scripts/fold_queries.py --panel runs/panel/panel.csv \
        --reference runs/panel/reference_queries.json \
        --out runs/panel/queries.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def build(reference: dict, panel: list[dict]) -> dict:
    ref_q = reference["queries"][0]
    by_id = {c["chain_id"]: c for c in ref_q["chains"]}
    missing = {"A", "B", "C", "D"} - set(by_id)
    if missing:
        raise ValueError(f"the reference is missing chain {sorted(missing)}")

    vh_len = len(by_id["B"]["sequence"])
    queries = []
    for row in panel:
        vh = row["vh_sequence"]
        if len(vh) != vh_len:
            raise ValueError(
                f"{row['id']}: VH is {len(vh)} residues, the template expects {vh_len}"
            )
        chains = []
        for cid in ("A", "B", "C", "D"):
            c = dict(by_id[cid])
            if cid == "B":
                c["sequence"] = vh
            elif cid == "A":
                c["sequence"] = row["vl_sequence"]
            chains.append(c)
        queries.append({"query_id": row["id"], "chains": chains})
    return {"queries": queries}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--groups", default="", help="comma-separated list, e.g. baseline,calibration")
    args = ap.parse_args()

    with args.panel.open(newline="") as fh:
        panel = list(csv.DictReader(fh))
    if args.groups:
        want = set(args.groups.split(","))
        panel = [r for r in panel if r["group"] in want]
    if not panel:
        raise SystemExit("the panel is empty")

    reference = json.loads(args.reference.read_text())
    out = build(reference, panel)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out))
    print(f"{len(out['queries'])} queries -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
