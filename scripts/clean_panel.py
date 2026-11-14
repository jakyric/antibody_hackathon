"""Panel for folding the liability-clean picks, with trastuzumab in the run.

The epitope axis is a comparison inside one run, so the parent is folded with
the designs, by the same model with the same seed and the same template. The
run is otherwise identical to the one that produced our existing numbers, so
the new loops and the old ones can be read on one scale.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--selected", type=Path, default=Path("results/pool/selected.csv"))
ap.add_argument("--out", type=Path, default=Path("candidates/clean_panel.csv"))
args = ap.parse_args()

with args.selected.open(newline="", encoding="utf-8") as fh:
    designs = list(csv.DictReader(fh))

rows = [
    {
        "id": "trastuzumab",
        "group": "baseline",
        "cdrh3": TZ_CDRH3_10,
        "kd_nM": "nan",
        "vh_sequence": TRASTUZUMAB_VH,
        "vl_sequence": TRASTUZUMAB_VL,
    }
]
for d in designs:
    rows.append(
        {
            "id": d["design_id"],
            "group": "design",
            "cdrh3": d["cdrh3"],
            "kd_nM": "nan",
            "vh_sequence": d["vh_sequence"],
            "vl_sequence": d["vl_sequence"],
        }
    )

args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open("w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} entries ({len(designs)} designs + the parent) -> {args.out}")
