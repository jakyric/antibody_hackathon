"""Search the liability-clean CDR-H3 space, and keep only what passes novelty.

WHY THIS EXISTS
  Our twenty designs were sampled from a position-weight matrix and inspected
  for liabilities afterwards. Nine of them introduce `FragH:DP`, a
  high-severity acid-labile motif trastuzumab does not carry, and only four
  reduce the high-severity count at all. A published comparison set scored
  with the same tool had nine of ten reducing it and none introducing FragH.

  The difference is the search, not the scoring. `abx.pool` enumerates the
  whole space of loops that drop all three of trastuzumab's own CDR-H3 motifs
  and add none, at exactly the four substitutions the novelty rule needs.
  This script scores that space and keeps the candidates that survive the
  rule we cannot negotiate.

THE ORDER OF THE FILTERS, AND WHY
  1. enumerate       complete, no seed, no sampling
  2. binding proxy   the trained FACS model, used to ORDER not to select --
                     it has no demonstrated skill inside a fixed mutation
                     count, which is exactly where every candidate here sits
  3. novelty         edit distance >= 3 from all 525,063 provided loops, by
                     the same vectorised check that judges the submission
  4. the official checker: the authority, run separately on the survivors

  Novelty is last of the cheap steps because it is the expensive one: it is
  the full provided set per candidate. Scoring first and cutting to a few
  thousand is what makes it affordable. Nothing is selected on the proxy --
  the cut is wide enough that the ranking inside it is not load-bearing.

USAGE
  python scripts/clean_pool.py --top 4000 --out results/pool/clean_pool.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from abx.data import TZ_CDRH3_10, provided_cdrh3s  # noqa: E402
from abx.gate import MIN_MUTATIONS_FROM_PROVIDED  # noqa: E402
from abx.model import load as load_model  # noqa: E402
from abx.model import onehot  # noqa: E402
from abx.pool import enumerate_clean, enumerate_contact_safe  # noqa: E402
from verify_novelty_full import edit_distances_same_length  # noqa: E402


def score_in_batches(net, loops: list[str], batch: int = 20000) -> np.ndarray:
    out = np.empty(len(loops), dtype=float)
    for start in range(0, len(loops), batch):
        chunk = loops[start : start + batch]
        out[start : start + len(chunk)] = net.logits(onehot(chunk))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-mut", type=int, default=4)
    ap.add_argument("--top", type=int, default=4000,
                    help="candidates carried into the novelty check")
    ap.add_argument("--out", type=Path, default=Path("results/pool/clean_pool.csv"))
    ap.add_argument(
        "--contact-safe",
        action="store_true",
        help="leave the four HER2-contacting residues alone (W99, G103, F104, "
             "Y105), which keeps the Trp oxidation liability and the paratope",
    )
    args = ap.parse_args(argv)

    loops = (enumerate_contact_safe(args.n_mut) if args.contact_safe
             else enumerate_clean(args.n_mut))
    kind = "contact-safe" if args.contact_safe else "liability-clean"
    print(f"{kind} loops at {args.n_mut} substitutions: {len(loops)}")
    if not loops:
        raise SystemExit("nothing to score")

    net = load_model()
    scores = score_in_batches(net, loops)
    order = np.argsort(-scores)
    keep = order[: args.top]
    print(f"carrying the top {len(keep)} into the novelty check "
          f"(proxy score {scores[keep].min():.4f} .. {scores[keep].max():.4f})")

    provided = sorted(provided_cdrh3s())
    by_len: dict[int, list[str]] = {}
    for p in provided:
        by_len.setdefault(len(p), []).append(p)
    packed = {
        L: np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)
        for L, seqs in by_len.items()
    }
    print(f"provided loops: {len(provided)} "
          f"(lengths {dict(sorted(Counter(len(p) for p in provided).items()))})")

    # The parent must come back at 0. A novelty check that calls trastuzumab
    # novel is measuring something else.
    control = min(
        int(edit_distances_same_length(TZ_CDRH3_10, packed[len(TZ_CDRH3_10)]).min()),
        99,
    )
    if control != 0:
        raise SystemExit(f"control failed: trastuzumab came back at {control}, not 0")
    print("control ok: trastuzumab collides with the provided data at 0")

    rows = []
    for n, i in enumerate(keep, 1):
        loop = loops[int(i)]
        best, nearest = 99, ""
        for L, seqs in by_len.items():
            if L == len(loop):
                d = edit_distances_same_length(loop, packed[L])
                j = int(d.argmin())
                if int(d[j]) < best:
                    best, nearest = int(d[j]), seqs[j]
            else:
                for s in seqs:
                    if abs(len(s) - len(loop)) >= best:
                        continue
                    from abx.gate import levenshtein

                    dd = levenshtein(loop, s)
                    if dd < best:
                        best, nearest = dd, s
        if best >= MIN_MUTATIONS_FROM_PROVIDED:
            rows.append(
                {
                    "cdrh3": loop,
                    "proxy_score": round(float(scores[int(i)]), 6),
                    "proxy_rank": n,
                    "min_edit_distance_to_provided": best,
                    "nearest_provided": nearest,
                }
            )
        if n % 500 == 0:
            print(f"  checked {n}/{len(keep)}, kept {len(rows)}")

    print(f"\n{len(rows)} of {len(keep)} clear the novelty rule "
          f"(>= {MIN_MUTATIONS_FROM_PROVIDED} edits from every provided loop)")
    dist = Counter(r["min_edit_distance_to_provided"] for r in rows)
    print(f"distance to the nearest provided loop: {dict(sorted(dist.items()))}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
