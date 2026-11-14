#!/usr/bin/env python3
"""What the measured data say about each of our designs, before any model gets a say.

Two questions:

1. Does the design touch a position where *all* measured single mutants lost badly?
   That is known counter-evidence, not a model score, and it should weigh more
   heavily than either of our two ranking axes.
2. How close to trastuzumab is the design on BLOSUM45? The baseline we have to beat.

    python scripts/precedent_report.py                     # results/precedent.csv
    python scripts/precedent_report.py --table results/ranked_submission.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10, cdrh3_from_vh, load_spr  # noqa: E402
from abx.precedent import (  # noqa: E402
    POSITION_LABELS,
    blosum45_similarity,
    design_evidence,
    position_precedent,
    single_mutant_table,
    wildtype_kd,
)

FIELDS = [
    "design_id", "cdrh3", "n_mutations", "mutations", "blosum45",
    "her2_contacts_touched", "n_her2_contacts_touched", "measured_singles",
    "worst_measured_fold", "against_measured_precedent", "structural_risk",
    "n_risk_flags",
]


def spr_rows() -> list[dict[str, str]]:
    spr = load_spr()
    return [
        {"cdrh3": str(c), "kd_nM": "" if k != k else str(k), "binding_result": str(r)}
        for c, k, r in zip(spr.cdrh3, spr.kd_nM, spr.binding_result, strict=True)
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", type=Path, default=ROOT / "results/ranked_submission.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "results/precedent.csv")
    args = ap.parse_args()

    rows = spr_rows()
    singles = single_mutant_table(rows)
    wt = wildtype_kd(rows)

    print(f"trastuzumab, median measured KD: {wt:.4f} nM")
    print(f"measured single mutants with quantified KD: {len(singles)}")
    print("\nposition  n measured   loss vs wild type")
    for i, label in enumerate(POSITION_LABELS):
        folds = position_precedent(singles, i)
        if folds:
            span = f"{min(folds):.0f}–{max(folds):.0f}x" if len(folds) > 1 else f"{folds[0]:.0f}x"
            verdict = "  ← all lost at least 10x" if min(folds) >= 10 else ""
            print(f"  {label:<6} {len(folds):>5}   {span}{verdict}")
        else:
            print(f"  {label:<6} {0:>5}   no measured single mutant")

    with args.table.open(newline="") as fh:
        designs = list(csv.DictReader(fh))

    out_rows = []
    skipped = []
    for d in designs:
        # The ranking table has a cdrh3 column; the submission table has only
        # the full VH. Both must be readable, so the report follows the submission.
        cdrh3 = d.get("cdrh3") or cdrh3_from_vh(d.get("heavy_sequence", ""))
        name = d.get("design_id") or d.get("sequence_id", "")
        if not cdrh3 or len(cdrh3) != len(TZ_CDRH3_10):
            skipped.append(f"{name} ({cdrh3 or 'no CDR-H3 found'})")
            continue
        out_rows.append({"design_id": name, **design_evidence(cdrh3, singles)})
    if skipped:
        print(f"\nskipped {len(skipped)}: {', '.join(skipped)}")
    # Trastuzumab itself as a reference row: it should have zero of everything.
    out_rows.append({"design_id": "trastuzumab", **design_evidence(TZ_CDRH3_10, singles)})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(out_rows)

    flagged = [r for r in out_rows if r["n_risk_flags"]]
    print(f"\n{len(out_rows) - 1} designs read from {args.table.name}")
    print(f"BLOSUM45 against trastuzumab, highest possible: {blosum45_similarity(TZ_CDRH3_10)}")
    print(f"designs with at least one risk flag: {len(flagged)}")
    for r in sorted(out_rows[:-1], key=lambda r: -r["n_risk_flags"]):
        if r["n_risk_flags"]:
            print(f"  {r['design_id']:<12} {r['cdrh3']}  blosum45={r['blosum45']:>3}  "
                  f"{r['against_measured_precedent']}{' | ' if r['against_measured_precedent'] and r['structural_risk'] else ''}{r['structural_risk']}")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
