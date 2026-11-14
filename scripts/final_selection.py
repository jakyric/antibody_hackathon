"""Choose the final candidates, on the basis that survived the audit.

WHY THIS REPLACES THE EARLIER RANKING
  `results/RANKING.md` ordered the twenty by aggregating a binding rank and a
  developability rank. The audit in `results/audit/` removed the first of
  those. All twenty designs sit at exactly 4 CDRH3 mutations, so edit distance
  -- the only scorer that tracks measured affinity -- is constant across them
  and orders nothing. Inside that stratum no scorer we have beats chance, and
  the term that actually dominated the old formula, `conservative_fraction`,
  has a negative point estimate against affinity in every subset tested.

  So we stop pretending to rank on affinity. What remains is real:

  1. developability, measured against trastuzumab with the organisers' own
     liability reference;
  2. diversity, which is the correct response to not being able to pick a
     winner.

SELECTION RULE
  - Exclude every design that introduces a high-severity motif trastuzumab
    does not carry. That is the FragH:DP group, nine designs. The trade they
    make -- drop trastuzumab's tolerated DG, add a novel acid-labile DP -- is
    not neutral, and it is the one developability judgement we are confident
    in.
  - Order the rest by developability, best first.
  - Report pairwise CDRH3 distances, so "diverse" is a number and not a hope.

WHAT WE DO NOT CLAIM
  The order within the selected set is NOT an affinity prediction. We cannot
  rank these against each other on binding and we say so. The order reflects
  developability and nothing else.
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 99


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", type=Path, default=Path("results/ranked_submission.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n", type=int, default=10)
    args = ap.parse_args(argv)

    with args.ranked.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    clean = [r for r in rows if int(r["dev_tier"]) != 3]
    dropped = [r for r in rows if int(r["dev_tier"]) == 3]

    # Best developability first; dev_rank already encodes tier then counts.
    clean.sort(key=lambda r: int(r["dev_rank"]))
    chosen = clean[: args.n]

    print(f"{len(rows)} designs -> excluded {len(dropped)} for introducing "
          f"a high-severity motif -> {len(clean)} eligible -> selected {len(chosen)}\n")

    print(f"{'#':>2}  {'design':<12} {'CDRH3':<12} {'liab':>4} {'d_hi':>5} "
          f"{'gained':<10} {'lost':<14} tier")
    for i, r in enumerate(chosen, 1):
        print(f"{i:>2}  {r['design_id']:<12} {r['cdrh3']:<12} "
              f"{r['n_liabilities']:>4} {int(r['d_high']):>+5} "
              f"{(r['h3_gained'] or '-'):<10} {(r['h3_lost'] or '-'):<14} "
              f"{r['dev_tier']}")

    # ---- diversity, as a number ------------------------------------------
    loops = [r["cdrh3"] for r in chosen]
    dists = [hamming(a, b) for a, b in combinations(loops, 2)]
    print(f"\npairwise CDRH3 distance among the {len(chosen)} selected: "
          f"min {min(dists)}, median {sorted(dists)[len(dists)//2]}, max {max(dists)}")
    close = [(a["design_id"], b["design_id"], hamming(a["cdrh3"], b["cdrh3"]))
             for a, b in combinations(chosen, 2)
             if hamming(a["cdrh3"], b["cdrh3"]) <= 3]
    if close:
        print("pairs within 3 mutations of each other (same bet twice):")
        for x, y, d in close:
            print(f"  {x} / {y}: {d}")
    else:
        print("no two selected designs are within 3 mutations of each other")

    # position-level variety: how many distinct residues appear per position
    print("\nresidue variety per CDRH3 position across the selection:")
    for p in range(len(loops[0])):
        seen = sorted({l[p] for l in loops})
        print(f"  pos {p+1:>2}: {len(seen)} distinct  {''.join(seen)}")

    fields = list(chosen[0].keys())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i, r in enumerate(chosen, 1):
            r = dict(r)
            r["rank"] = i
            w.writerow(r)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
