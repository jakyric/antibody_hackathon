"""Pick a diverse handful from the liability-clean, novelty-passing pool.

WHY DIVERSITY FIRST AND NOT THE SCORE
  Every candidate in this pool sits at exactly four substitutions from
  trastuzumab, which is the regime where our binding proxy has no
  demonstrated skill: inside a fixed mutation count, edit distance is
  constant and the trained model beat nothing. Ordering thousands of
  candidates by a score that cannot separate them is false precision.

  What we can defend is spread. Four guaranteed lab slots are four
  independent bets, and two designs that differ in one residue are one bet
  written twice. So the selection maximises the minimum pairwise distance
  between chosen loops, and uses the proxy only to break ties and to choose
  the seed. If the proxy carries signal we could not demonstrate, starting
  from its best candidate keeps that upside; if it carries none, we have
  still bought spread.

USAGE
  python scripts/clean_select.py --n 6 --out results/pool/selected.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.gate import hamming, levenshtein  # noqa: E402
from abx.pool import CONTACT_POSITIONS, is_clean, is_clean_outside  # noqa: E402
from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402


def greedy_spread(loops: list[str], n: int) -> list[str]:
    """Farthest-point selection: each pick is the one furthest from the rest.

    `loops` arrives in proxy order, so the first pick is the proxy's best and
    every later tie is broken towards it.
    """
    chosen = [loops[0]]
    while len(chosen) < n and len(chosen) < len(loops):
        best, best_d = None, -1
        for cand in loops:
            if cand in chosen:
                continue
            d = min(levenshtein(cand, c) for c in chosen)
            if d > best_d:
                best, best_d = cand, d
        chosen.append(best)
    return chosen


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pool", type=Path, default=Path("results/pool/clean_pool.csv"))
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--consider", type=int, default=400,
                    help="how far down the proxy order the spread may reach")
    ap.add_argument("--out", type=Path, default=Path("results/pool/selected.csv"))
    ap.add_argument(
        "--contact-safe",
        action="store_true",
        help="the pool keeps the four HER2-contacting residues, so the Trp "
             "oxidation liability on W99 is expected and is not a reason to drop",
    )
    ap.add_argument(
        "--require-residue",
        action="append",
        default=[],
        metavar="POS:AA",
        help="keep only loops with this residue at this 1-based position, e.g. "
             "10:D. Use for a position the measured data says not to move.",
    )
    args = ap.parse_args(argv)

    with args.pool.open(newline="", encoding="utf-8") as fh:
        pool = list(csv.DictReader(fh))
    print(f"pool: {len(pool)} loops, novel")

    # The pool was enumerated under an earlier, shorter motif list; the official checker
    # then showed that list was missing DT. Re-applying the current filter
    # here means the widened list governs the picks without re-running the
    # novelty check, which is the expensive step and is unaffected by it.
    before = len(pool)
    if args.contact_safe:
        pool = [r for r in pool
                if is_clean_outside(r["cdrh3"], CONTACT_POSITIONS)]
    else:
        pool = [r for r in pool if is_clean(r["cdrh3"])]
    print(f"      {len(pool)} still clean under the widened motif list "
          f"({before - len(pool)} dropped)")
    if not pool:
        raise SystemExit("the widened motif list emptied the pool")

    for spec in args.require_residue:
        pos, aa = spec.split(":")
        i = int(pos) - 1
        was = len(pool)
        pool = [r for r in pool if r["cdrh3"][i] == aa]
        print(f"      {len(pool)} keep {aa} at position {pos} "
              f"({was - len(pool)} dropped)")
        if not pool:
            raise SystemExit(f"requiring {spec} emptied the pool")

    by_loop = {r["cdrh3"]: r for r in pool}
    loops = [r["cdrh3"] for r in pool][: args.consider]
    chosen = greedy_spread(loops, args.n)

    pair = [levenshtein(a, b) for i, a in enumerate(chosen) for b in chosen[i + 1:]]
    print(f"\nchosen {len(chosen)}; pairwise edit distance "
          f"min {min(pair)}, median {sorted(pair)[len(pair) // 2]}, max {max(pair)}")

    rows = []
    for i, loop in enumerate(chosen, 1):
        src = by_loop[loop]
        rows.append(
            {
                "design_id": f"tz-clean-{i:02d}",
                "cdrh3": loop,
                "n_mutations_from_trastuzumab": hamming(loop, TZ_CDRH3_10),
                "min_edit_distance_to_provided": src["min_edit_distance_to_provided"],
                "nearest_provided": src["nearest_provided"],
                "proxy_rank": src["proxy_rank"],
                "proxy_score": src["proxy_score"],
                "vh_sequence": TRASTUZUMAB_VH.replace(TZ_CDRH3_10, loop),
                "vl_sequence": TRASTUZUMAB_VL,
            }
        )
        print(f"  {rows[-1]['design_id']}  {loop}  proxy rank {src['proxy_rank']:>5}"
              f"  novelty {src['min_edit_distance_to_provided']}")

    for r in rows:
        if TZ_CDRH3_10 in r["vh_sequence"] or r["cdrh3"] not in r["vh_sequence"]:
            raise SystemExit(f"{r['design_id']}: the loop was not grafted correctly")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
