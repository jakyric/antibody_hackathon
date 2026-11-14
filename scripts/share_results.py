#!/usr/bin/env python3
"""Write the shareable result tables to the shared repo.

Why a script and not copying by hand: the branch is public. Which
columns go out should be a rule that can be read and tested, not something
someone remembered to remove that one time.

THE RULE
  Numbers about methods are shared. Sequences are not.
  Anything that can assemble a candidate — CDR-H3, VH, VL, scFv — is removed.
  Design IDs stay: "tz-opt-05 got ipTM 0.936" is verifiable
  for the team and useless to another team.

    python scripts/share_results.py                 # write to the clone
    python scripts/share_results.py --check         # only say what would happen
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Columns that can assemble a candidate. Removed no matter which file they
# show up in, so a new column in a new table doesn't slip out with it.
SEQUENCE_COLUMNS = {
    "cdrh3", "vh_sequence", "vl_sequence", "scfv_sequence",
    "heavy_sequence", "light_sequence", "sequence", "sequences",
    "binder_sequence", "designed_sequences", "vh", "vl", "mutations",
}

# (source in results/, destination in the clone). Only these. A file that is
# not listed here is not shared, no matter how harmless it looks.
TABLES = [
    ("folding/calibration.csv", "results/folding/calibration.csv"),
    ("folding/agreement.csv", "results/folding/agreement.csv"),
    ("folding/vs_trastuzumab.csv", "results/folding/vs_trastuzumab.csv"),
    ("folding/scores.csv", "results/folding/scores.csv"),
    ("epitope/epitope.csv", "results/epitope/epitope.csv"),
]


def strip_sequences(rows: list[dict], fields: list[str]) -> tuple[list[dict], list[str]]:
    keep = [f for f in fields if f.lower() not in SEQUENCE_COLUMNS]
    dropped = [f for f in fields if f not in keep]
    return [{k: r.get(k, "") for k in keep} for r in rows], dropped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=ROOT / "results")
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    clone = Path(json.loads((ROOT / "scripts" / "share_map.json").read_text())
                 ["clone"].replace("~", os.path.expanduser("~")))

    written = 0
    for src_rel, dest_rel in TABLES:
        src = args.results / src_rel
        if not src.exists():
            print(f"  skipping {src_rel} (does not exist)")
            continue
        with src.open(newline="") as fh:
            reader = csv.DictReader(fh)
            rows, fields = list(reader), list(reader.fieldnames or [])
        rows, dropped = strip_sequences(rows, fields)
        note = f"  removed {dropped}" if dropped else ""
        print(f"{'would write' if args.check else 'wrote'} {dest_rel}"
              f"  {len(rows)} rows{note}")
        if args.check:
            continue
        dest = clone / dest_rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else fields)
            w.writeheader()
            w.writerows(rows)
        written += 1
    if not args.check:
        print(f"\n{written} tables to {clone}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
