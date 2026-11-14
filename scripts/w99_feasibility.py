#!/usr/bin/env python3
"""Are there legal designs that keep W99, and what does it cost?

BACKGROUND
  All six submitted designs replace W99. This is not an oversight:
  W99 *is* trastuzumab's TrpOx liability, and `src/abx/pool.py` removes it
  on purpose. But W99 is also a HER2 contact, and among the 150 measured
  variants in the legal region the sub-nM rate is 6/22 for those that keep W99
  versus 0/12 for those that replace it, at exactly four mutations (Fisher
  one-sided p = 0.055; see results/PRECEDENT.md).

  So the replacement is a real trade-off, not a mistake: developability
  against measured affinity. This script answers what has to be settled before
  the trade-off can be made at all: whether the other side exists.

THE QUESTION
  Are there loops at four mutations that keep W99 (and preferably the other
  three contacts), introduce no new liability motifs, and lie at least three
  mutations from all 525,063 provided CDR-H3s?

  Four mutations, not three: three give 10/13 = 76.9% identity on the official
  checker's 13-residue numbering, and the requirement is below 70%.

    python scripts/w99_feasibility.py --out results/w99_pool.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.gate import ProvidedIndex, count_unpaired_cysteines  # noqa: E402
from abx.pool import (  # noqa: E402
    CONTACT_POSITIONS,
    enumerate_contact_safe,
    is_clean,
)
from abx.precedent import blosum45_similarity, design_evidence  # noqa: E402
from abx.precedent import single_mutant_table  # noqa: E402


def spr_rows() -> list[dict[str, str]]:
    from abx.data import load_spr  # noqa: PLC0415

    spr = load_spr()
    return [
        {"cdrh3": str(c), "kd_nM": "" if k != k else str(k), "binding_result": str(r)}
        for c, k, r in zip(spr.cdrh3, spr.kd_nM, spr.binding_result, strict=True)
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mut", type=int, default=4)
    ap.add_argument("--want", type=int, default=300,
                    help="how many passed candidates we need")
    ap.add_argument("--out", type=Path, default=ROOT / "results/w99_pool.csv")
    args = ap.parse_args()

    print("building index of the provided CDR-H3s ...")
    idx = ProvidedIndex()
    print(f"  {idx.n:,} sequences\n")

    # The novelty check costs ~1.9 ms per sequence, so a full scan of all the
    # enumerations takes hours. It is not necessary either: we do not need the
    # whole legal space, only to know that it exists and to pull the best ones
    # out of it. So the candidates are sorted first, and the check is run down
    # the list until we have enough, with the observed share reported so the
    # size of the space can be estimated honestly.
    loops = enumerate_contact_safe(args.n_mut, frozen=CONTACT_POSITIONS)
    label = "all four HER2 contacts frozen"
    if not loops:
        loops = enumerate_contact_safe(args.n_mut, frozen=frozenset({0}))
        label = "W99 frozen only"
    print(f"enumeration, {label}: {len(loops):,} loops without new motifs")

    loops = [s for s in loops if count_unpaired_cysteines(s) == 0]
    loops.sort(key=lambda s: -blosum45_similarity(s))
    print(f"  {len(loops):,} without unpaired cysteines, sorted by BLOSUM45\n")

    keep, tested = [], 0
    for s in loops:
        tested += 1
        if idx.is_novel(s):
            keep.append(s)
            if len(keep) >= args.want:
                break
    rate = len(keep) / tested
    print(f"novelty check: {len(keep):,} passed out of {tested:,} tried "
          f"({100 * rate:.1f}%)")
    print(f"estimated size of the legal space: ~{int(rate * len(loops)):,} "
          f"of {len(loops):,}")

    if not keep:
        print("\nNo legal loop keeps W99. The trade-off is settled by the requirements.")
        return 0

    print(f"\nwriting the {len(keep):,} from \"{label}\"")
    singles = single_mutant_table(spr_rows())
    rows = []
    for s in keep:
        ev = design_evidence(s, singles)
        rows.append({
            "cdrh3": s,
            "blosum45": ev["blosum45"],
            "mutations": ev["mutations"],
            "n_her2_contacts_touched": ev["n_her2_contacts_touched"],
            "motif_free": int(is_clean(s)),
            "against_measured_precedent": ev["against_measured_precedent"],
            "min_distance_to_provided": idx.min_distance(s),
        })
    rows.sort(key=lambda r: (-int(r["min_distance_to_provided"]), -int(r["blosum45"])))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    motif_free = sum(r["motif_free"] for r in rows)
    print(f"  of them, {motif_free:,} are fully motif-free "
          f"(the rest keep TrpOx:W, which is the price of holding on to W99)")
    print(f"  margin to provided data: "
          f"{min(int(r['min_distance_to_provided']) for r in rows)}–"
          f"{max(int(r['min_distance_to_provided']) for r in rows)} mutations")
    print(f"\ntrastuzumab itself: BLOSUM45 {blosum45_similarity(TZ_CDRH3_10)}")
    print("the ten best (largest margin, then most similar to trastuzumab):")
    for r in rows[:10]:
        print(f"  {r['cdrh3']}  bl45={r['blosum45']:>3}  "
              f"distance={r['min_distance_to_provided']}  {r['mutations']}")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
