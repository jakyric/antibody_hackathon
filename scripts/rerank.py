"""Re-rank the designs on two axes instead of one.

WHY THIS EXISTS
  The original ranking in `results/designs.csv` rests on the binding proxy:
  closeness to trastuzumab, which is the only signal that transferred to
  measured SPR affinity, and transferred weakly (Spearman +0.23). We now have
  a second axis from the organisers' own liability reference, and it is clean:
  it compares each design against a molecule we know survives manufacturing.

HOW THE TWO ARE COMBINED
  Not by a weighted sum. Weights we cannot calibrate would invent precision
  we do not have: no measurement tells us that one liability is worth N
  mutations of edit distance. We use rank aggregation instead -- each design's
  position on each axis, added. It needs no units and no weights, and it
  rewards being good on both axes rather than maximal on one.

  Sorting strictly by developability was the first thing we tried, and it was
  wrong. It promoted designs from binding rank 15 and 18 into the top four on
  the strength of a liability count. Liability motifs speak to chemical
  degradation and expression risk; they are NOT a thermostability prediction,
  and thermostability is the judged axis. Letting them override the one signal
  we have for a judged property over-trusts the proxy.

  Designs that introduce a high-severity motif trastuzumab does not carry are
  demoted below everything else regardless of rank sum. That is the FragH:DP
  group. Several look neutral on the count because they dropped trastuzumab's
  DG and added DP; that trade is not neutral. DG is demonstrably tolerated in
  a licensed product, DP is risk we would introduce with no such evidence.
  Same nominal severity, very different standing.

WHAT THIS IS NOT
  Not an affinity prediction, and not a thermostability prediction. Both of
  those are measured by the lab after the event. This orders candidates by
  how good a bet they look on two computational proxies, nothing more.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def tier_of(row: dict) -> int:
    introduces_high = bool(row.get("h3_gained_high", "").strip())
    d_high = int(row["d_n_high"])
    if introduces_high:
        return 3
    return 1 if d_high < 0 else 2


TIER_NOTE = {
    1: "removes a high-severity liability without introducing one",
    2: "developability unchanged against trastuzumab",
    3: "introduces a high-severity motif trastuzumab does not carry",
}


def developability_rank(rows: list[dict]) -> dict[str, int]:
    """Order on the developability axis alone: fewer high-severity motifs
    first, then fewer liabilities overall."""
    ordered = sorted(
        rows,
        key=lambda r: (tier_of(r), int(r["d_n_high"]), int(r["n_total"]),
                       r["sequence_id"]),
    )
    return {r["sequence_id"]: i for i, r in enumerate(ordered, 1)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--liabilities", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    with args.designs.open(newline="", encoding="utf-8") as fh:
        designs = list(csv.DictReader(fh))
    with args.liabilities.open(newline="", encoding="utf-8") as fh:
        liab = {r["sequence_id"]: r for r in csv.DictReader(fh)}

    dev_rank = developability_rank(
        [l for sid, l in liab.items() if sid != "tz_reference"]
    )

    merged = []
    for d in designs:
        l = liab.get(d["design_id"])
        if l is None:
            raise SystemExit(f"no liability row for {d['design_id']}")
        tier = tier_of(l)
        merged.append({
            "design_id": d["design_id"],
            "binding_rank": int(d["rank"]),
            "dev_rank": dev_rank[d["design_id"]],
            "rank_sum": int(d["rank"]) + dev_rank[d["design_id"]],
            "dev_tier": tier,
            "dev_note": TIER_NOTE[tier],
            "n_liabilities": l["n_total"],
            "d_high": l["d_n_high"],
            "h3_gained": l["h3_gained"],
            "h3_lost": l["h3_lost"],
            "cdrh3": d["cdrh3"],
            "n_mut_from_trastuzumab": d["n_mutations_from_trastuzumab"],
            "min_mut_from_provided": d["min_mutations_from_provided"],
            "model_percentile": d["model_percentile"],
            "scfv_sequence": d["scfv_sequence"],
            "vh_sequence": d["vh_sequence"],
            "vl_sequence": d["vl_sequence"],
        })

    # Designs introducing new high-severity risk go last whatever their sum.
    merged.sort(key=lambda r: (r["dev_tier"] == 3, r["rank_sum"],
                               r["binding_rank"]))
    for i, r in enumerate(merged, 1):
        r["rank"] = i
        r["rank_change"] = r["binding_rank"] - i

    fields = ["rank", "design_id", "binding_rank", "dev_rank", "rank_sum",
              "rank_change", "dev_tier", "dev_note", "n_liabilities", "d_high",
              "h3_gained", "h3_lost", "cdrh3", "n_mut_from_trastuzumab",
              "min_mut_from_provided", "model_percentile", "vh_sequence",
              "vl_sequence", "scfv_sequence"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(merged)

    print(f"{'new':>3}  {'design':<12} {'bind':>5} {'dev':>4} {'sum':>4} "
          f"{'move':>5} {'tier':>5}  note")
    for r in merged:
        move = f"{r['rank_change']:+d}" if r["rank_change"] else "-"
        print(f"{r['rank']:>3}  {r['design_id']:<12} {r['binding_rank']:>5} "
              f"{r['dev_rank']:>4} {r['rank_sum']:>4} {move:>5} "
              f"{r['dev_tier']:>5}  {r['dev_note']}")

    print(f"\ntop 4 (the designs every team is guaranteed to get tested):")
    for r in merged[:4]:
        print(f"  {r['rank']}. {r['design_id']}  CDRH3 {r['cdrh3']}  "
              f"was rank {r['binding_rank']}")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
