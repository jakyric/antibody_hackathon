#!/usr/bin/env python3
"""Build the folding panel for the de novo designs: the designs themselves, plus the anchors.

The anchors are what make the scores readable. Without them an ipTM is just a number.
  baseline     trastuzumab itself, measured KD 0.13 nM
  anchor-*     three SPR variants with measured KD, from strongest to weakest
All are folded in the same batch, with the same template, model and seed as the designs, so
the comparison is within a single run and not across runs.

The panel is written in the same format that scripts/fold_panel.py uses, so
scripts/fold_queries.py can build queries.json from it.

    python scripts/denovo_panel.py --candidates results/denovo/candidates_pass.csv \
        --out runs/denovo-panel/panel.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10, load_spr  # noqa: E402
from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

FIELDS = ["id", "group", "cdrh3", "kd_nM", "vh_sequence", "vl_sequence"]


def anchors(n: int = 3) -> list[dict]:
    """n measured variants evenly spaced over log10(KD), plus trastuzumab.

    Only loops as long as trastuzumab's: the template omits exactly ten residues
    where CDR-H3 goes, so a shorter loop would have folded a different geometry.
    """
    spr = load_spr()
    keep = spr.has_kd
    uniq: dict[str, list] = {}
    for c, k, h, l in zip(spr.cdrh3[keep], spr.kd_nM[keep], spr.vh[keep], spr.vl[keep]):
        if len(str(c)) != len(TZ_CDRH3_10) or str(c) == TZ_CDRH3_10:
            continue
        uniq.setdefault(str(c), [[], str(h), str(l)])[0].append(float(k))

    items = sorted(((c, float(np.median(v[0])), v[1], v[2]) for c, v in uniq.items()),
                   key=lambda t: t[1])
    picks = [items[round(i * (len(items) - 1) / (n - 1))] for i in range(n)] if n > 1 else items[:1]

    rows = [{
        "id": "trastuzumab", "group": "baseline", "cdrh3": TZ_CDRH3_10,
        "kd_nM": "0.1345", "vh_sequence": TRASTUZUMAB_VH, "vl_sequence": TRASTUZUMAB_VL,
    }]
    for i, (cdrh3, kd, vh, vl) in enumerate(picks):
        rows.append({"id": f"anchor-{i:02d}", "group": "anchor", "cdrh3": cdrh3,
                     "kd_nM": f"{kd:.4g}", "vh_sequence": vh, "vl_sequence": vl})
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--anchors", type=int, default=3)
    ap.add_argument("--top-n", type=int, default=0,
                    help="fold only the N best by the design model's own numbers")
    args = ap.parse_args(argv)

    with args.candidates.open(newline="") as fh:
        cands = list(csv.DictReader(fh))
    if not cands:
        raise SystemExit("no candidates in the table")

    def design_score(rec: dict) -> tuple[float, float]:
        """The design model's own numbers: does it hit the hotspots, and did it
        hold the sequence constraints? This is not a binding prediction, just a
        reasonable order in which to spend GPU time."""
        def num(key: str) -> float:
            try:
                return float(rec.get(key, "") or "nan")
            except ValueError:
                return float("nan")
        h, s = num("hotspot_obedience"), num("sequence_adherence")
        return (-1e9 if h != h else h, -1e9 if s != s else s)

    if args.top_n:
        cands = sorted(cands, key=design_score, reverse=True)[:args.top_n]

    rows = anchors(args.anchors)
    vh_len = len(TRASTUZUMAB_VH)
    vl_len = len(TRASTUZUMAB_VL)
    skipped = 0
    for rec in cands:
        vh, vl = rec["vh_sequence"].strip(), rec["vl_sequence"].strip()
        if len(vh) != vh_len or len(vl) != vl_len:
            skipped += 1      # the template fits one length; another would fold something else
            continue
        rows.append({"id": rec["design_id"], "group": "denovo", "cdrh3": rec.get("cdrh3", ""),
                     "kd_nM": "nan", "vh_sequence": vh, "vl_sequence": vl})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_design = sum(1 for r in rows if r["group"] == "denovo")
    print(f"{len(rows)} rows -> {args.out} ({n_design} designs, {len(rows) - n_design} anchors)")
    if skipped:
        print(f"  {skipped} designs dropped: VH/VL has a different length than the template")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
