"""Merge both tracks into the single ranked CSV the organisers asked for.

THE TWO DECISIONS THIS FILE ENCODES

1. ONE FILE, NOT TWO.
   The channel instruction says to send "a CSV file", singular, named after the
   team. The task description says a team submits one ranked list of up to 20
   designs, and that the official checker decides the track per row -- a design is never
   moved between tracks, it is classified. So a single file is what the process
   expects, and two files risk one being read as a correction of the other.

2. THE ORDER ALTERNATES BETWEEN THE TRACKS, STARTING 2 + 2.
   This is the decision that actually costs something, so here is the whole
   reasoning.

   The prize is split equally between the two tracks and each is ranked
   separately. Only the top N of this one list is tested, and N is 4 unless we
   place. Four optimisation designs on top therefore means nothing tested in
   the de novo track, and a guaranteed zero in a track worth half the pot.

   The counter-argument is that our de novo VHHs are the longer shot: they are
   designed from scratch and we have no evidence any of them binds. That is
   true. But it is true for every team's de novo entries, and the tier is
   scored relative to trastuzumab, not relative to perfection. A track where
   everyone's designs are speculative is a track with a low bar to place in.
   Forfeiting it to add a third and fourth optimisation design -- which sit at
   the same four mutations as the first two, in a regime where no scorer we
   have can separate them -- buys very little.

   So: two from each track in the guaranteed four.

   Within the optimisation track the guaranteed slots go to STABILISED
   designs. Under the brief's rule a property counts only if it beats the
   baseline ("ties and regressions are hollow"), so an unmodified parent can
   only win thermostability by luck, while a stabilised design has a
   literature-backed reason to win it. That makes the stabilised rows the
   better bet for A- and B-tier, which is where this competition will actually
   be decided.

   The two unmodified parents follow immediately at positions 5 and 6, so that
   placing anywhere at all completes the matched pairs and makes the R66G
   result readable. The two stabilised designs with no parent in the list go
   last: they cannot support a controlled comparison, so they are the weakest
   rows to spend a slot on.

USAGE
  python scripts/build_delivery.py   -> results/delivery/fabsolutely-binding.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The order is written out by hand, because a rank is a judgement and a sort
# key would hide it. Each entry is (source track, sequence_id, why this slot).
ORDER = [
    # ── the four that are tested whatever happens ──────────────────────
    ("opt", "tz-clean-02-stab",
     "optimisation, stabilised. Epitope 0.95 against trastuzumab with the "
     "crystal control at 1.00, and the best BLOSUM45 similarity of the "
     "designs we folded cleanly."),
    ("denovo", "team5-denovo-vhh-01",
     "de novo. The design both folding models agree on: best averaged "
     "percentile, 2.7 A refold to its designed pose, passes the organisers' "
     "oracle under OpenFold3."),
    ("opt", "tz-clean-03-stab",
     "optimisation, stabilised, a loop distinct from tz-clean-02. Epitope "
     "1.00, the highest of the six."),
    ("denovo", "team5-denovo-vhh-05",
     "de novo. Best single Boltz-2 interface score and the largest epitope "
     "overlap, with a 1.4 A refold; the two folding models disagree about it, "
     "which is why it is second and not first."),

    # ── 5-6: complete the matched pairs, so any placement makes R66G readable
    ("opt", "tz-clean-02",
     "CONTROL: tz-clean-02 without R66G. Completes the pair with slot 1 so the "
     "lab result separates the loop from the mutation."),
    ("opt", "tz-clean-03",
     "CONTROL: tz-clean-03 without R66G. Completes the pair with slot 3."),

    # ── 7-10 ───────────────────────────────────────────────────────────
    ("opt", "tz-clean-01-stab",
     "optimisation, stabilised. Highest BLOSUM45 similarity to trastuzumab of "
     "all six, but the lowest epitope retention (0.85), so it ranks below the "
     "two above rather than above them."),
    ("denovo", "team5-denovo-vhh-04",
     "de novo. Second-largest epitope overlap (0.56)."),
    ("opt", "tz-clean-07-stab",
     "optimisation, stabilised, fourth distinct loop. Epitope 0.90."),
    ("denovo", "team5-denovo-vhh-03",
     "de novo. Tight pLDDT (94) but a loose 8.2 A refold."),

    # ── 11-16 ──────────────────────────────────────────────────────────
    ("opt", "tz-clean-01",
     "CONTROL: tz-clean-01 without R66G, completing the pair with slot 7."),
    ("opt", "tz-clean-07",
     "CONTROL: tz-clean-07 without R66G, completing the pair with slot 9."),
    ("denovo", "team5-denovo-vhh-02",
     "de novo. Moderate scores, kept for spread across the designed loops."),
    ("denovo", "team5-denovo-vhh-06",
     "de novo. The weakest interface scores of the six."),
    ("opt", "tz-clean-06-stab",
     "optimisation, stabilised, epitope 1.00 -- but no unmodified parent is in "
     "this list, so testing it alone cannot say what R66G did."),
    ("opt", "tz-clean-04-stab",
     "optimisation, stabilised, epitope 0.95, and likewise unpaired."),
]

FIELDS = ["sequence_id", "kind", "heavy_sequence", "light_sequence",
          "design_description"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--optimisation", type=Path,
                    default=ROOT / "results/SUBMISSION.csv")
    ap.add_argument("--denovo", type=Path,
                    default=ROOT / "results/denovo/denovo_submission.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results/delivery/fabsolutely-binding.csv")
    args = ap.parse_args(argv)

    src: dict[str, dict] = {}
    for track, path in (("opt", args.optimisation), ("denovo", args.denovo)):
        for r in csv.DictReader(path.open(newline="", encoding="utf-8")):
            src[r["sequence_id"]] = {"_track": track, **r}

    missing = [sid for _, sid, _ in ORDER if sid not in src]
    if missing:
        raise SystemExit(f"not in either source file: {', '.join(missing)}")
    unplaced = sorted(set(src) - {sid for _, sid, _ in ORDER})
    if unplaced:
        raise SystemExit(
            "every design must be given a rank or deliberately dropped; "
            f"unplaced: {', '.join(unplaced)}"
        )

    rows = []
    for rank, (track, sid, why) in enumerate(ORDER, 1):
        r = src[sid]
        if r["_track"] != track:
            raise SystemExit(f"{sid} is {r['_track']}, the plan says {track}")
        rows.append({
            "sequence_id": sid,
            "kind": r["kind"],
            "heavy_sequence": r["heavy_sequence"],
            "light_sequence": r.get("light_sequence", ""),
            # The organisers read this column. Lead with the rank and the
            # reason, keep whatever the track's own builder wrote after it.
            "design_description": f"Rank {rank}. {why} | {r['design_description']}",
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    n_opt = sum(1 for t, _, _ in ORDER if t == "opt")
    n_dn = len(ORDER) - n_opt
    print(f"{len(rows)} rows -> {args.out}")
    print(f"  optimisation {n_opt}, de novo {n_dn}")
    top4 = [t for t, _, _ in ORDER[:4]]
    print(f"  guaranteed four: {top4.count('opt')} optimisation, "
          f"{top4.count('denovo')} de novo")
    for i, (t, sid, _) in enumerate(ORDER, 1):
        mark = "  <- guaranteed" if i <= 4 else ""
        print(f"  {i:>2}. {sid:<24} {t}{mark}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
