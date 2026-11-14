"""Rank the liability-clean picks on the two axes that still separate them.

All of these designs carry the same liability profile -- every one drops all
three of trastuzumab's CDR-H3 motifs and adds none -- so developability,
which ordered the old list, is now a constant and orders nothing. What is
left is:

  epitope   the predicted HER2 contact set against trastuzumab's, same model,
            same seed, same template. This is a structural measurement of the
            axis the lab judges, and the parent reproduces its own crystal
            epitope exactly in this run, so the comparison has a floor.

  proxy     the trained FACS model. It has no demonstrated skill inside a
            fixed mutation count, which is where every candidate sits, so it
            is used as the second of two ranks and never alone.

Combined by rank aggregation, as before: positions added, no weights, no
invented precision. A design whose predicted epitope has drifted materially
from the parent's is dropped rather than ranked, because an optimisation-track
design that no longer competes with trastuzumab for the site has lost the
axis outright.

USAGE
  python scripts/clean_rank.py --min-epitope 0.94 --out results/pool/ranked.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.gate import levenshtein  # noqa: E402


def average_ranks(values: list[float], higher_is_better: bool) -> list[float]:
    """Ranks from 1, with ties sharing the average of the places they fill.

    Averaging matters here: three designs with an identical epitope overlap
    must not be ordered by the accident of their row order, and giving them
    all the same rank keeps the second axis as the thing that separates them.
    """
    order = sorted(values, reverse=higher_is_better)
    out = []
    for v in values:
        places = [i + 1 for i, x in enumerate(order) if x == v]
        out.append(sum(places) / len(places))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selected", type=Path,
                    default=Path("results/pool/selected8.csv"))
    ap.add_argument("--epitope", type=Path,
                    default=Path("results/pool/epitope/epitope.csv"))
    ap.add_argument("--min-epitope", type=float, default=0.94)
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--out", type=Path, default=Path("results/pool/ranked.csv"))
    args = ap.parse_args(argv)

    with args.selected.open(newline="", encoding="utf-8") as fh:
        designs = {r["design_id"]: r for r in csv.DictReader(fh)}
    with args.epitope.open(newline="", encoding="utf-8") as fh:
        epi = {r["id"]: r for r in csv.DictReader(fh)}

    baseline = epi.get("trastuzumab")
    if baseline is None or float(baseline["overlap_vs_crystal"]) < 0.99:
        raise SystemExit(
            "the parent does not reproduce its own crystal epitope in this run; "
            "without that floor a design's overlap means nothing"
        )
    print(f"floor ok: trastuzumab reproduces its crystal epitope at "
          f"{baseline['overlap_vs_crystal']}")

    rows = []
    for did, d in designs.items():
        e = epi.get(did)
        if e is None:
            print(f"  {did}: not folded, left out")
            continue
        overlap = float(e["overlap_vs_tz_same_model"])
        if overlap < args.min_epitope:
            print(f"  {did}: epitope {overlap} below {args.min_epitope}, dropped")
            continue
        rows.append({**d, "epitope": overlap, "proxy_rank_in_pool": int(d["proxy_rank"])})

    if len(rows) < args.n:
        raise SystemExit(f"only {len(rows)} designs survive the epitope floor")

    e_ranks = average_ranks([r["epitope"] for r in rows], higher_is_better=True)
    p_ranks = average_ranks([-r["proxy_rank_in_pool"] for r in rows],
                            higher_is_better=True)
    for r, er, pr in zip(rows, e_ranks, p_ranks):
        r["epitope_rank"] = er
        r["proxy_rank_here"] = pr
        r["rank_sum"] = er + pr
    rows.sort(key=lambda r: (r["rank_sum"], r["proxy_rank_in_pool"]))
    rows = rows[: args.n]

    print(f"\n{'rank':>4}  {'id':<13} {'cdrh3':<12} {'epitope':>7} {'proxy':>7} "
          f"{'sum':>5}")
    for i, r in enumerate(rows, 1):
        print(f"{i:>4}  {r['design_id']:<13} {r['cdrh3']:<12} "
              f"{r['epitope']:>7} {r['proxy_rank_in_pool']:>7} {r['rank_sum']:>5.1f}")

    loops = [r["cdrh3"] for r in rows]
    d = [levenshtein(a, b) for i, a in enumerate(loops) for b in loops[i + 1:]]
    print(f"\npairwise edit distance: min {min(d)}, "
          f"median {sorted(d)[len(d) // 2]}, max {max(d)}")

    cols = ["design_id", "cdrh3", "n_mutations_from_trastuzumab",
            "min_edit_distance_to_provided", "nearest_provided", "proxy_rank",
            "proxy_score", "epitope", "rank_sum", "vh_sequence", "vl_sequence"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r[c] for c in cols})
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
