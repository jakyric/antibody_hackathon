#!/usr/bin/env python3
"""Rank the VHH designs on the refold scores from several folding models.

refold-score gives the organizers' own numbers per design: `pass` with a reason,
`min_ipae` (lower is better), `iptm` and `binder_plddt`, plus how far the
binding domain moved from the pose it was designed in.

THE RULE
  Each model gives a percentile for min_ipae (inverted, so high is good). The
  mean over the models is the binding axis. The models are not pooled: the
  percentile is computed within a single model run. If the two models
  disagree, it shows in the spread, and then we do not tell the designs apart.

  `pass` from the oracle is not a ranking but a threshold: a design that no
  model lets pass is still listed, but flagged.

The output has the same shape as scripts/denovo_rank.py, so
scripts/denovo_select.py can select from it.

    python scripts/denovo_vhh_rank.py --metrics boltz2=runs/<k>/refold/scored-boltz2/refold_metrics.csv \
        --epitope boltz2=results/denovo/epitope_vhh_boltz2.csv --out results/denovo/vhh
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


def pairs(items: list[str]) -> dict[str, list[Path]]:
    """model -> all the files for that model. One file per design run."""
    out: dict[str, list[Path]] = {}
    for item in items:
        name, _, path = item.partition("=")
        if not path:
            raise SystemExit(f"write model=path, not {item!r}")
        out.setdefault(name, []).append(Path(path))
    return out


def run_prefix(path: Path) -> str:
    """The first 15 characters of the run id in the path, or an empty string.

    WHY THIS IS NEEDED: each design run numbers its examples from zero, so
    `her2_vhh_denovo_0` exists in ALL the runs. Without the run in the key,
    different designs merge into one, and the ranking becomes nonsense.
    The prefix is the same one design_id uses in candidates_pass.csv.
    """
    m = re.search(r"\d{8}T\d{6}Z", str(path))
    return m.group(0)[:15] if m else ""


def num(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def read_metrics(path: Path, prefix: str = "") -> dict[str, dict]:
    """design_id -> the row from refold_metrics.csv.

    design_id is "<run id>-<example_id>" when the prefix is known, so designs
    from different runs do not collide. Several seeds per design: keep the
    best min_ipae, and count how many rows the design has.
    """
    out: dict[str, dict] = {}
    with path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            qid = f"{prefix}-{rec['query_id']}" if prefix else rec["query_id"]
            rec["design_id"] = qid
            cur = out.get(qid)
            if cur is None or num(rec.get("min_ipae")) < num(cur.get("min_ipae")):
                rec["_n"] = (cur or {}).get("_n", 0) + 1
                out[qid] = rec
            else:
                cur["_n"] = cur.get("_n", 1) + 1
    return out


def read_all(paths: list[Path]) -> dict[str, dict]:
    """All the runs for one model, merged with the run in the key."""
    merged: dict[str, dict] = {}
    for p in paths:
        merged.update(read_metrics(p, run_prefix(p)))
    return merged


def percentiles(values: dict[str, float], higher_is_better: bool) -> dict[str, float]:
    items = [(k, v) for k, v in values.items() if not math.isnan(v)]
    order = sorted(items, key=lambda kv: kv[1], reverse=not higher_is_better)
    n = len(order)
    # low min_ipae should give a high percentile. A single design has nothing to
    # be ranked against; it gets 100 and must be read together with models_scored.
    if n == 1:
        return {order[0][0]: 100.0}
    return {k: round(100 * i / (n - 1), 1) for i, (k, _) in enumerate(order)} if n else {}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", action="append", required=True, metavar="MODEL=CSV")
    ap.add_argument("--epitope", action="append", default=[], metavar="MODEL=CSV")
    ap.add_argument("--candidates", type=Path,
                    help="candidates_pass.csv, for CDR3 and origin in the table")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    metrics = {m: read_all(ps) for m, ps in pairs(args.metrics).items()}
    epi: dict[str, dict] = {}
    for m, ps in pairs(args.epitope).items():
        rows: dict[str, dict] = {}
        for p in ps:
            prefix = run_prefix(p)
            with p.open(newline="") as fh:
                for r in csv.DictReader(fh):
                    key = f"{prefix}-{r['query_id']}" if prefix else r["query_id"]
                    rows[key] = r
        epi[m] = rows

    cands: dict[str, dict] = {}
    if args.candidates:
        with args.candidates.open(newline="") as fh:
            for r in csv.DictReader(fh):
                # design_id is "<run id>-<example_id>"; refold uses example_id
                cands[r["design_id"]] = r
                cands.setdefault(r["design_id"].split("-", 1)[-1], r)

    pct = {m: percentiles({q: num(r.get("min_ipae")) for q, r in rows.items()},
                          higher_is_better=False)
           for m, rows in metrics.items()}

    ids = sorted({q for rows in metrics.values() for q in rows})
    rows_out = []
    for qid in ids:
        rec = cands.get(qid, {})
        row = {"design_id": qid, "group": "denovo",
               "cdrh3": rec.get("cdrh3", ""),
               "weight": rec.get("weight", ""), "seed": rec.get("seed", ""),
               "run_id": rec.get("run_id", ""),
               "kind": rec.get("kind", "vhh")}
        pcts, passes = [], []
        for model in sorted(metrics):
            m = metrics[model].get(qid, {})
            row[f"{model}_min_ipae"] = m.get("min_ipae", "")
            row[f"{model}_iptm"] = m.get("iptm", "")
            row[f"{model}_binder_plddt"] = m.get("binder_plddt", "")
            row[f"{model}_rmsd_to_designed_pose"] = m.get("target_aligned_binder_rmsd", "")
            row[f"{model}_oracle_pass"] = m.get("pass", "")
            p = pct[model].get(qid)
            row[f"{model}_percentile"] = p if p is not None else ""
            if p is not None:
                pcts.append(p)
            if str(m.get("pass", "")).lower() == "true":
                passes.append(model)
            if model in epi:
                e = epi[model].get(qid, {})
                row[f"{model}_epitope_jaccard"] = e.get("jaccard_to_1n8z", "")
        row["models_scored"] = len(pcts)
        row["oracle_pass_models"] = ";".join(passes)
        row["binding_percentile_mean"] = round(float(np.mean(pcts)), 1) if pcts else ""
        row["binding_percentile_spread"] = round(float(max(pcts) - min(pcts)), 1) if len(pcts) > 1 else ""
        epi_model = sorted(epi)[0] if epi else ""
        row["epitope_jaccard"] = row.get(f"{epi_model}_epitope_jaccard", "") if epi_model else ""
        rows_out.append(row)

    rows_out.sort(key=lambda r: -(r["binding_percentile_mean"] if r["binding_percentile_mean"] != "" else -1))
    for i, r in enumerate(rows_out, 1):
        r["rank"] = i

    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / "scores.csv"
    with dest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
        w.writeheader()
        w.writerows(rows_out)
    n_pass = sum(1 for r in rows_out if r["oracle_pass_models"])
    print(f"{len(rows_out)} designs ranked -> {dest}")
    print(f"  {n_pass} designs pass the oracle in at least one model")
    for r in rows_out[:5]:
        print(f"  {r['rank']}. {r['design_id']} percentile {r['binding_percentile_mean']} "
              f"min_ipae " + ", ".join(f"{m}={r.get(f'{m}_min_ipae', '')}" for m in sorted(metrics)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

