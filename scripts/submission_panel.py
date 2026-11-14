"""Build a folding panel from whatever is currently in the submission table.

WHY THIS EXISTS
  The epitope column is one of the three judged axes, and it is the only one
  we can put a number on before the lab runs. That number comes from folding
  the design against HER2 and comparing its contact set with trastuzumab's.

  The danger is drift: the submission gets rebuilt, the design ids and the
  loops change, and the epitope numbers left in results/ silently describe
  molecules we are no longer submitting. That happened once already. This
  script reads the submission as it stands now, so the panel can never be
  about an older generation of designs.

  Only the distinct heavy chains are folded. A `-stab` row differs from its
  parent by one light-chain residue, so folding both doubles the cost for a
  comparison that `scripts/r66g_analyse.py` already covers; the epitope value
  is carried by the loop.

USAGE
  python scripts/submission_panel.py --out runs/subpanel/panel.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

PREFIX, SUFFIX = "YYCSR", "YWGQGT"


def cdrh3(heavy: str) -> str:
    start = heavy.index(PREFIX) + len(PREFIX)
    return heavy[start:heavy.index(SUFFIX, start)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submission", type=Path, default=Path("results/SUBMISSION.csv"))
    ap.add_argument("--out", type=Path, default=Path("runs/subpanel/panel.csv"))
    args = ap.parse_args(argv)

    rows = list(csv.DictReader(args.submission.open(newline="", encoding="utf-8")))

    # Trastuzumab first: every epitope number is an overlap against it, so a
    # panel without it produces comparisons with nothing to compare to.
    panel = [{
        "id": "trastuzumab", "group": "baseline", "cdrh3": cdrh3(TRASTUZUMAB_VH),
        "kd_nM": "0.1345", "vh_sequence": TRASTUZUMAB_VH, "vl_sequence": TRASTUZUMAB_VL,
    }]

    seen: dict[str, str] = {}
    for r in rows:
        vh = r["heavy_sequence"].strip()
        if vh in seen:
            continue
        base = r["sequence_id"].replace("-stab", "")
        seen[vh] = base
        panel.append({
            "id": base, "group": "design", "cdrh3": cdrh3(vh), "kd_nM": "",
            "vh_sequence": vh, "vl_sequence": r["light_sequence"].strip(),
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0]))
        w.writeheader()
        w.writerows(panel)

    print(f"{len(rows)} submitted rows -> {len(panel)-1} distinct heavy chains "
          f"+ trastuzumab = {len(panel)} queries")
    for p in panel:
        print(f"  {p['id']:20s} {p['group']:9s} {p['cdrh3']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
