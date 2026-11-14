#!/usr/bin/env python3
"""Build the panel we fold: a calibration split with measured KD, and our designs.

The point of the calibration split: we know KD for 182 variants. If a folding
model's ipTM / interface-PAE ranks those 182 correctly, we can rank designs on
it. If not, we cannot, and then we have measured that instead of assuming it.

The panel is the same for all folding models, so the numbers can be compared.

    python scripts/fold_panel.py --out runs/panel
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abx.data import TZ_CDRH3_10, load_spr  # noqa: E402
from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

N_CALIB = 24  # calibration variants, evenly spread over log KD


def calibration_panel(n: int = N_CALIB) -> list[dict]:
    """n SPR variants evenly spread over log10(KD), plus trastuzumab itself."""
    spr = load_spr()
    keep = spr.has_kd
    cdrh3 = spr.cdrh3[keep]
    kd = spr.kd_nM[keep]
    vh = spr.vh[keep]
    vl = spr.vl[keep]

    # one row per unique CDRH3 (median KD where it was measured several times)
    uniq: dict[str, list] = {}
    for c, k, h, l in zip(cdrh3, kd, vh, vl):
        uniq.setdefault(str(c), [[], str(h), str(l)])[0].append(float(k))
    # The template we fold against leaves out exactly ten residues where CDR-H3 sits.
    # The SPR table also has shorter loops; they do not fit the template and must go,
    # otherwise we would be comparing two different geometries. The number that
    # drop out is shown in the log.
    n_all = len(uniq)
    uniq = {c: v for c, v in uniq.items() if len(c) == len(TZ_CDRH3_10)}
    if n_all != len(uniq):
        print(f"  dropped {n_all - len(uniq)} variants with CDR-H3 "
              f"not equal to {len(TZ_CDRH3_10)} residues (does not fit the template)")
    rows = [
        {"cdrh3": c, "kd_nM": float(np.median(v[0])), "vh": v[1], "vl": v[2]}
        for c, v in uniq.items()
    ]
    rows.sort(key=lambda r: r["kd_nM"])

    # evenly in rank across the whole range, deterministic
    idx = np.unique(np.linspace(0, len(rows) - 1, n).round().astype(int))
    panel = []
    for rank, i in enumerate(idx):
        r = rows[i]
        panel.append(
            {
                "id": f"calib-{rank:02d}",
                "group": "calibration",
                "cdrh3": r["cdrh3"],
                "kd_nM": r["kd_nM"],
                "vh_sequence": r["vh"],
                "vl_sequence": r["vl"],
            }
        )
    # trastuzumab as the baseline, always included, always first
    tz = [r for r in rows if r["cdrh3"] == TZ_CDRH3_10]
    panel = [
        {
            "id": "trastuzumab",
            "group": "baseline",
            "cdrh3": TZ_CDRH3_10,
            "kd_nM": tz[0]["kd_nM"] if tz else float("nan"),
            "vh_sequence": TRASTUZUMAB_VH,
            "vl_sequence": TRASTUZUMAB_VL,
        }
    ] + [p for p in panel if p["cdrh3"] != TZ_CDRH3_10]
    return panel


def design_panel(designs_csv: Path) -> list[dict]:
    out = []
    with designs_csv.open(newline="") as fh:
        for row in csv.DictReader(fh):
            out.append(
                {
                    "id": row["design_id"],
                    "group": "design",
                    "cdrh3": row["cdrh3"],
                    "kd_nM": float("nan"),
                    "vh_sequence": row["vh_sequence"],
                    "vl_sequence": row["vl_sequence"],
                }
            )
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--n-calib", type=int, default=N_CALIB)
    ap.add_argument("--no-designs", action="store_true")
    args = ap.parse_args()

    panel = calibration_panel(args.n_calib)
    if not args.no_designs and args.designs.exists():
        panel += design_panel(args.designs)

    args.out.mkdir(parents=True, exist_ok=True)
    cols = ["id", "group", "cdrh3", "kd_nM", "vh_sequence", "vl_sequence"]
    with (args.out / "panel.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(panel)
    (args.out / "panel.json").write_text(json.dumps(panel, indent=1))

    kds = [p["kd_nM"] for p in panel if p["kd_nM"] == p["kd_nM"]]
    print(f"{len(panel)} entries -> {args.out/'panel.csv'}")
    print(f"  calibration: {sum(p['group']=='calibration' for p in panel)}"
          f"  designs: {sum(p['group']=='design' for p in panel)}")
    print(f"  KD range: {min(kds):.3g} .. {max(kds):.3g} nM")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
