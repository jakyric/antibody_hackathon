"""Exhaustive novelty check: every submitted design against EVERY provided CDRH3.

WHY THIS EXISTS
  `scripts/build_reference.py` builds the official checker reference from a 4,579-row
  subset of the provided data: every sequence within 4 mutations of a design,
  plus random controls. The argument for the shortcut is sound -- a sequence
  more than 4 mutations away cannot set a design's minimum distance -- but it
  is an argument, not a measurement, and it rests on two assumptions that have
  never been tested:

    1. that every provided CDRH3 is exactly 10 residues, so that Hamming
       distance is defined for all of them (`build_reference.py` assigns
       distance 99 to any row of a different length, silently dropping it);
    2. that no provided sequence of a DIFFERENT length sits within 3 edits of
       a design, which Hamming cannot see at all.

  This script makes no such assumption. It compares all ten submitted designs
  against all 525,063 provided CDRH3s, by Hamming where lengths match and by
  Levenshtein where they do not, and reports the true minimum.

  It is the evidence that the gate covers the whole dataset.

USAGE
  python scripts/verify_novelty_full.py --out results/gate/novelty_full.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abx.data import TZ_CDRH3_10, provided_cdrh3s  # noqa: E402
from abx.gate import MIN_MUTATIONS_FROM_PROVIDED, levenshtein  # noqa: E402
from abx.sequences import CDRH3_PREFIX  # noqa: E402

SUFFIX = "YWGQGT"


def cdrh3_from_heavy(heavy: str) -> str:
    """The 10-residue loop in the provided tables' convention.

    The heavy chain reads ...YYCSR <10-mer> Y WGQGT..., and the trailing Y is
    not part of the `cdrh3` column even though Kabat counts it. Slicing on the
    anchors keeps us in the tables' convention, which is the one the novelty
    rule is measured in.
    """
    start = heavy.index(CDRH3_PREFIX) + len(CDRH3_PREFIX)
    end = heavy.index(SUFFIX, start)
    return heavy[start:end]


def edit_distances_same_length(query: str, rows: np.ndarray) -> np.ndarray:
    """Levenshtein distance from `query` to every row of an (N, L) matrix.

    WHY NOT HAMMING
      Equal length does not mean the alignment is positional. `YPGRYFYTYD`
      and `YPAGRFYTYD` are both 10-mers differing at four positions, but one
      is the other with a residue inserted and another deleted: two edits.
      The official checker aligns, so it sees 2 where a positional comparison sees 4 --
      and 2 is below the three-mutation floor. Scoring these by Hamming
      reports a design as novel that the official checker rejects, and a rejected
      design is never made.

    The dynamic-programming table is walked one query residue at a time with
    the whole column held as a vector over the N provided sequences, so the
    cost is (len(query) x L) numpy operations rather than N python loops.
    """
    n_rows, width = rows.shape
    q = np.frombuffer(query.encode(), dtype=np.uint8)
    prev = np.tile(np.arange(width + 1, dtype=np.int16), (n_rows, 1))
    cur = np.empty_like(prev)
    for i in range(1, len(q) + 1):
        cur[:, 0] = i
        for j in range(1, width + 1):
            cost = (rows[:, j - 1] != q[i - 1]).astype(np.int16)
            np.minimum(prev[:, j] + 1, cur[:, j - 1] + 1, out=cur[:, j])
            np.minimum(cur[:, j], prev[:, j - 1] + cost, out=cur[:, j])
        prev, cur = cur, prev
    return prev[:, width]


def min_distance_exhaustive(
    query: str, same_len: np.ndarray, other_len: list[str]
) -> tuple[int, str]:
    """True minimum distance of `query` to every provided sequence.

    `same_len` is an (N, L) uint8 matrix of the provided sequences whose length
    equals len(query); those are compared by edit distance, vectorised.
    `other_len` holds the rest, compared by Levenshtein one at a time. Returns
    the minimum and the sequence that attains it.
    """
    best, best_seq = len(query) + 99, ""
    if same_len.size:
        dists = edit_distances_same_length(query, same_len)
        i = int(dists.argmin())
        best, best_seq = int(dists[i]), same_len[i].tobytes().decode()
    for seq in other_len:
        # A length difference already costs |len(a)-len(b)| edits, so anything
        # that far off cannot beat the current best; skip it.
        if abs(len(seq) - len(query)) >= best:
            continue
        d = levenshtein(query, seq)
        if d < best:
            best, best_seq = d, seq
    return best, best_seq


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submission", type=Path, default=Path("results/SUBMISSION.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/gate/novelty_full.csv"))
    ap.add_argument(
        "--required",
        type=int,
        default=MIN_MUTATIONS_FROM_PROVIDED,
        help="minimum mutations the official checker enforces (3)",
    )
    args = ap.parse_args(argv)

    provided = sorted(provided_cdrh3s())
    lengths = Counter(len(p) for p in provided)
    print(f"provided CDRH3s: {len(provided)} unique")
    print(f"length distribution: {dict(sorted(lengths.items()))}")

    with args.submission.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    # Trastuzumab itself is the negative control: its loop IS in the provided
    # data, so it must come back at distance 0. A run where it does not is a
    # broken comparison, not a passing design.
    rows.append(
        {
            "sequence_id": "trastuzumab_CONTROL",
            "heavy_sequence": "",
            "_cdrh3": TZ_CDRH3_10,
        }
    )

    queries = [
        (r["sequence_id"], r.get("_cdrh3") or cdrh3_from_heavy(r["heavy_sequence"]))
        for r in rows
    ]

    by_len: dict[int, list[str]] = {}
    for p in provided:
        by_len.setdefault(len(p), []).append(p)
    packed = {
        L: np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)
        for L, seqs in by_len.items()
    }

    out_rows = []
    failures = []
    for sid, cdrh3 in queries:
        same = packed.get(len(cdrh3), np.empty((0, len(cdrh3)), dtype=np.uint8))
        other = [p for L, seqs in by_len.items() if L != len(cdrh3) for p in seqs]
        dist, nearest = min_distance_exhaustive(cdrh3, same, other)
        ok = dist >= args.required
        if sid.endswith("_CONTROL"):
            ok = dist == 0  # the control must collide
        out_rows.append(
            {
                "sequence_id": sid,
                "cdrh3": cdrh3,
                "min_distance_to_provided": dist,
                "nearest_provided": nearest,
                "required": args.required,
                "passes": "yes" if ok else "NO",
            }
        )
        print(f"  {sid:24s} {cdrh3:12s} min={dist}  nearest={nearest}  {'ok' if ok else 'FAIL'}")
        if not ok:
            failures.append(sid)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nwrote {args.out}")

    if failures:
        print(f"FAILED: {', '.join(failures)}")
        return 1
    print("all rows behaved as required (designs >= %d, control at 0)" % args.required)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
