"""Build a reference database for the official checker from the provided HER2 data.

WHY THIS EXISTS
  The novelty rule says a design must differ from every CDRH3 in the two
  provided tables (525,063 unique sequences). But the reference bundled with
  the official checker does not contain them: in our gate run every design's nearest known
  CDRH3 came back as trastuzumab itself, at distance 4. The rule that actually
  constrains us has therefore never been checked with the organisers' tool.

  Our own gate puts the nearest provided sequence at distance 3, exactly the
  minimum the official checker requires. Zero margin. If the official checker were to read a CDRH3
  differently than we do -- it numbers 13 residues (SRWGGDGFYAMDY) where we
  work on the 10-residue core -- a design we believe is at 3 could land at 2
  and be rejected. A rejected design is not tested at all.

WHAT IT BUILDS
  reference-ingest wants assay CSVs with a heavy chain column, because it
  reads the CDR3 off the heavy chain with immunum. The provided tables hold
  bare 10-residue loops, so each one is grafted into trastuzumab's VH.

  Numbering half a million heavy chains is the slow step, and we do not need
  it: only sequences NEAR a design can determine that design's minimum
  distance. Everything at distance > 4 is irrelevant to the verdict. So the
  reference holds every provided sequence within --radius of any design, plus
  a random sample as a control that the ingest path works on ordinary rows.

  This is a faithful test of the binding constraint, not a shortcut around it:
  if no provided sequence sits closer than 3 from a design, that design's
  minimum distance is 3 whether or not the other 500,000 rows are present.
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from abx.data import provided_cdrh3s, vh_from_cdrh3  # noqa: E402
from abx.gate import levenshtein  # noqa: E402
from abx.sequences import TRASTUZUMAB_VL  # noqa: E402
from verify_novelty_full import edit_distances_same_length  # noqa: E402


def designs_cdrh3(path: Path) -> list[tuple[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return [(r["design_id"], r["cdrh3"]) for r in csv.DictReader(fh)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--radius",
        type=int,
        default=4,
        help="keep provided sequences within this many mutations of a design",
    )
    ap.add_argument("--sample", type=int, default=2000, help="random control rows")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    designs = designs_cdrh3(args.designs)
    provided = sorted(provided_cdrh3s())
    print(f"{len(provided)} provided CDRH3s, {len(designs)} designs")

    # Distances are EDIT distances, not positional ones, even where the two
    # loops are the same length. Two 10-mers that differ at three positions
    # can be one insertion and one deletion apart, which the official checker sees and a
    # Hamming comparison does not; a reference built on Hamming can therefore
    # leave out the very neighbour that decides a design's verdict.
    by_len: dict[int, list[str]] = {}
    for p in provided:
        by_len.setdefault(len(p), []).append(p)
    packed = {
        L: np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)
        for L, seqs in by_len.items()
    }

    near: dict[str, int] = {}
    for did, d in designs:
        best = 99
        for L, seqs in by_len.items():
            if L == len(d):
                dists = edit_distances_same_length(d, packed[L])
            else:
                if abs(L - len(d)) > args.radius and abs(L - len(d)) >= best:
                    continue
                dists = np.array([levenshtein(d, s) for s in seqs])
            best = min(best, int(dists.min()))
            for i in np.flatnonzero(dists <= args.radius):
                seq = seqs[int(i)]
                near[seq] = min(near.get(seq, 99), int(dists[i]))
        print(f"  {did}: nearest provided sequence is {best} mutation(s) away")

    rng = random.Random(args.seed)
    control = rng.sample(provided, min(args.sample, len(provided)))
    keep = dict(near)
    for p in control:
        keep.setdefault(p, -1)

    print(f"\nreference holds {len(keep)} sequences "
          f"({len(near)} within {args.radius} of a design, "
          f"{len(keep) - len(near)} random controls)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sequence_id", "heavy_sequence", "light_sequence", "cdrh3_10mer"])
        for i, (loop, dist) in enumerate(sorted(keep.items())):
            tag = f"near{dist}" if dist >= 0 else "control"
            w.writerow([f"provided_{tag}_{i:06d}", vh_from_cdrh3(loop),
                        TRASTUZUMAB_VL, loop])

    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
