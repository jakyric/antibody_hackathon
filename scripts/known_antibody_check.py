"""Check our designs against the reference the official checker actually ships.

WHY THIS EXISTS
  The official checker's `team-validate --reference-db <ours>` **overrides** the
  bundled reference -- the help says so in one line that is easy to read past.
  Every validation we ran therefore measured our designs against the provided
  HER2 data and nothing else.

  The bundled reference holds 610,109 CDRH3s: the 524,308 from the FACS table,
  and also uspto (38,308), plabdab (38,197), pdb (5,943), plabdab_nano,
  therasabdab, swissprot and several target sets. Those ~84,000 known and
  patented loops are part of the rule -- a design must differ from *every
  known* CDRH3 -- and we had never measured against them.

  `novelty-report` leaves `is_novel_enough` null whether or not the bundle is
  passed, and the parquet directory is not accepted as a reference bundle, so
  the tool would not answer the question. This does, directly, on the file the
  tool ships.

WHAT IT COMPARES
  The official checker numbers CDRH3 with IMGT and gets 13 residues for a
  trastuzumab variant: SR + the 10-mer + Y. The reference column is in the same
  convention, so the design is padded to match before measuring. Distances are
  edit distances, because that is what the official checker uses -- the mistake
  that nearly cost us a design once already.

USAGE
  python scripts/known_antibody_check.py --table results/SUBMISSION.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.gate import levenshtein  # noqa: E402
from abx.sequences import CDRH3_PREFIX  # noqa: E402
from verify_novelty_full import edit_distances_same_length  # noqa: E402

BUNDLED = Path("reference_data/official_checker_cdrh3.parquet")
SUFFIX = "YWGQGT"
REQUIRED = 3  # mutations a trastuzumab variant needs, per the official checker


def imgt_cdrh3(heavy: str) -> str:
    """The 13-residue CDRH3 the official checker numbers, from a heavy chain."""
    start = heavy.index(CDRH3_PREFIX) + len(CDRH3_PREFIX)
    loop = heavy[start : heavy.index(SUFFIX, start)]
    return "SR" + loop + "Y"


def nearest(query: str, by_len: dict, packed: dict, meta: dict) -> tuple[int, str, str]:
    best, best_seq, best_src = len(query) + 99, "", ""
    for L, seqs in by_len.items():
        if L == len(query):
            d = edit_distances_same_length(query, packed[L])
            i = int(d.argmin())
            if int(d[i]) < best:
                best, best_seq, best_src = int(d[i]), seqs[i], meta[L][i]
        else:
            if abs(L - len(query)) >= best:
                continue
            for s, src in zip(seqs, meta[L]):
                if abs(len(s) - len(query)) >= best:
                    continue
                dd = levenshtein(query, s)
                if dd < best:
                    best, best_seq, best_src = dd, s, src
    return best, best_seq, best_src


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, default=Path("results/SUBMISSION.csv"))
    ap.add_argument("--reference", type=Path, default=BUNDLED)
    ap.add_argument("--out", type=Path,
                    default=Path("results/gate/known_antibody_check.csv"))
    args = ap.parse_args(argv)

    ref = pl.read_parquet(args.reference)
    print(f"bundled reference: {ref.height} CDRH3s")
    print("sources: " + ", ".join(
        f"{s} {n}" for s, n in
        ref["source"].value_counts().sort("count", descending=True)
        .head(8).iter_rows()))

    seqs_all = ref["cdrh3"].to_list()
    srcs_all = ref["source"].to_list()
    by_len: dict[int, list[str]] = {}
    meta: dict[int, list[str]] = {}
    for s, src in zip(seqs_all, srcs_all):
        by_len.setdefault(len(s), []).append(s)
        meta.setdefault(len(s), []).append(src)
    packed = {
        L: np.frombuffer("".join(v).encode(), dtype=np.uint8).reshape(len(v), L)
        for L, v in by_len.items()
    }

    with args.table.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    # Trastuzumab is the control. Its own loop is in this reference, so it must
    # come back at 0. A run where it does not is measuring the wrong thing.
    queries = [(r["sequence_id"], imgt_cdrh3(r["heavy_sequence"])) for r in rows]
    queries.append(("trastuzumab_CONTROL", "SR" + TZ_CDRH3_10 + "Y"))

    out, failures = [], []
    for sid, q in queries:
        d, seq, src = nearest(q, by_len, packed, meta)
        ok = (d == 0) if sid.endswith("_CONTROL") else (d >= REQUIRED)
        mark = "ok" if ok else "FAIL"
        print(f"  {sid:<24} {q}  min={d}  nearest={seq} [{src}]  {mark}")
        out.append({"sequence_id": sid, "imgt_cdrh3": q, "min_edit_distance": d,
                    "nearest": seq, "nearest_source": src, "passes": ok})
        if not ok:
            failures.append(sid)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {args.out}")

    non_facs = Counter(
        r["nearest_source"] for r in out if not r["sequence_id"].endswith("_CONTROL")
    )
    print(f"nearest-source breakdown: {dict(non_facs)}")
    if failures:
        print(f"FAILED: {', '.join(failures)}")
        return 1
    print("every design is at least 3 edits from every known CDRH3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
