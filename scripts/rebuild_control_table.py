"""Rebuild the negative-control table: the submission plus trastuzumab itself.

Trastuzumab must FAIL the novelty rule. Keeping it in a table of its own means
every gate run has something in it that we know the answer to, so a run where
everything passes can be told apart from a run where the check is switched off.

The file has to be rebuilt whenever the submission changes, or the control run
silently validates an old set of designs.

USAGE
  python scripts/rebuild_control_table.py --submission results/SUBMISSION.csv --out SUB_with_tz.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submission", type=Path, default=Path("results/SUBMISSION.csv"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(args.submission.open(newline="", encoding="utf-8")))
    fields = list(rows[0])
    control = {k: "" for k in fields}
    control.update({
        "sequence_id": "tz_reference",
        "kind": "paired",
        "heavy_sequence": TRASTUZUMAB_VH,
        "light_sequence": TRASTUZUMAB_VL,
        "design_description": "NEGATIVE CONTROL: trastuzumab itself. This row must fail.",
    })
    rows.append(control)

    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows ({len(rows)-1} designs + trastuzumab) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
