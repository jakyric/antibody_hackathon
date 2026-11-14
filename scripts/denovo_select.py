#!/usr/bin/env python3
"""Pick the few designs we actually submit, using a rule that is written in code.

Only a handful get made in the lab, so every pick is expensive. The rule:

  1. only designs that passed the gate (locally and in the official checker)
  2. only designs where the epitope holds: the overlap with trastuzumab's own
     HER2 residues is at least as large as that of the worst of the anchors, so
     we do not submit something that binds somewhere else
  3. sort by mean percentile of interface ipTM across the
     folding models
  4. spread the picks: at most one design per CDR-H3 cluster (clustered on
     edit distance), so we do not submit four variants of the same thing

Each pick is written out with its justification, so the order can be verified.

    python scripts/denovo_select.py --scores results/denovo/scores.csv \
        --candidates results/denovo/candidates_pass.csv --n 4 \
        --out results/denovo/selection.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.gate import levenshtein  # noqa: E402

CLUSTER_DISTANCE = 4     # CDR-H3s closer than this count as the same idea


def num(value, default=float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def cluster(ids_seqs: list[tuple[str, str]], cutoff: int = CLUSTER_DISTANCE) -> dict[str, int]:
    """Simple greedy clustering on edit distance. The input order decides
    which one becomes the first of a cluster, so pass them in ranked order."""
    centres: list[tuple[str, int]] = []
    out: dict[str, int] = {}
    for ident, seq in ids_seqs:
        for centre, cid in centres:
            if levenshtein(seq, centre) < cutoff:
                out[ident] = cid
                break
        else:
            cid = len(centres)
            centres.append((seq, cid))
            out[ident] = cid
    return out


def select(rows: list[dict], n: int, epitope_floor: float) -> list[dict]:
    ranked = [r for r in rows if r.get("group") == "denovo"]
    ranked.sort(key=lambda r: -num(r.get("binding_percentile_mean"), -1))

    clusters = cluster([(r["design_id"], r.get("cdrh3", "")) for r in ranked])
    picked: list[dict] = []
    used: set[int] = set()
    for r in ranked:
        ej = num(r.get("epitope_jaccard"))
        if ej == ej and ej < epitope_floor:
            r["decision"] = f"not picked: epitope overlap {ej} below the anchor floor {epitope_floor}"
            continue
        cid = clusters[r["design_id"]]
        if cid in used:
            r["decision"] = f"not picked: same CDRH3 cluster ({cid}) as a higher-ranked design"
            continue
        if len(picked) >= n:
            r["decision"] = "not picked: outside the top n"
            continue
        used.add(cid)
        r["decision"] = (f"picked: rank {len(picked) + 1}, binding percentile "
                         f"{r.get('binding_percentile_mean')}, cluster {cid}")
        r["cdrh3_cluster"] = cid
        picked.append(r)
    for r in ranked:
        r.setdefault("cdrh3_cluster", clusters[r["design_id"]])
    return picked


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, required=True,
                    help="candidates_pass.csv: the sequences and their origin are taken from here")
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--epitope-floor", type=float, default=-1.0,
                    help="smallest allowed overlap with trastuzumab's epitope; "
                         "negative means the floor is read from the anchors")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    with args.scores.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    with args.candidates.open(newline="") as fh:
        cands: dict[str, dict] = {}
        for r in csv.DictReader(fh):
            # design_id is "<run-id>-<example_id>"; the refold tables use
            # only example_id, so both must look up the same row
            cands[r["design_id"]] = r
            cands.setdefault(r["design_id"].split("-", 1)[-1], r)

    floor = args.epitope_floor
    if floor < 0:
        anchor_vals = [num(r.get("epitope_jaccard")) for r in rows
                       if r.get("group") in ("anchor", "baseline")]
        anchor_vals = [v for v in anchor_vals if v == v]
        floor = min(anchor_vals) if anchor_vals else 0.0
        print(f"the epitope floor is set by the anchors: {floor}")

    picked = select(rows, args.n, floor)
    if not picked:
        raise SystemExit("no design made it through the selection")

    out_rows = []
    for r in picked:
        c = cands.get(r["design_id"], {})
        out_rows.append({
            "design_id": r["design_id"],
            "rank": len(out_rows) + 1,
            "vh_sequence": c.get("vh_sequence", ""),
            "vl_sequence": c.get("vl_sequence", ""),
            "cdrh3": r.get("cdrh3", ""),
            "weight": c.get("weight", ""),
            "seed": c.get("seed", ""),
            "run_id": c.get("run_id", ""),
            "identity_cdr_pct": c.get("identity_cdr_pct", ""),
            "cdrh3_min_dist_provided": c.get("cdrh3_min_dist_provided", ""),
            "binding_percentile_mean": r.get("binding_percentile_mean", ""),
            "binding_percentile_spread": r.get("binding_percentile_spread", ""),
            "epitope_jaccard": r.get("epitope_jaccard", ""),
            "ddg_vs_trastuzumab": r.get("ddg_vs_trastuzumab", ""),
            "cdrh3_cluster": r.get("cdrh3_cluster", ""),
            "decision": r.get("decision", ""),
        })
    missing = [r["design_id"] for r in out_rows if not r["vh_sequence"]]
    if missing:
        raise SystemExit(f"missing sequence for {missing}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    print(f"{len(out_rows)} designs selected -> {args.out}")
    for r in out_rows:
        print(f"  {r['rank']}. {r['design_id']} {r['cdrh3']} {r['decision']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
