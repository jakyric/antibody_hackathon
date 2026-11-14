"""Positive control: does the official checker's novelty gate actually read our reference?

WHY THIS EXISTS
  A gate that passes everything looks identical to a gate that is switched
  off. After building a reference from the provided HER2 data, the official
  checker still reported trastuzumab as the nearest known CDRH3 for every
  design, at distance 4 -- even though our own index finds provided sequences
  at distance 3. Either the official checker measures something we do not, or
  the reference was not consulted. "All 20 pass" means nothing until we know
  which.

  So we hand it a design that MUST be rejected: one mutation away from a
  sequence that is in the reference, and far from trastuzumab. If that design
  passes, the reference is not being used and every green result built on it
  is void.

  This is the same logic as including trastuzumab itself in the submission
  table. A control that must fail is the only way to prove a filter is on.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abx.data import provided_cdrh3s, vh_from_cdrh3  # noqa: E402
from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.gate import hamming  # noqa: E402
from abx.sequences import TRASTUZUMAB_VL  # noqa: E402

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--reference-csv",
        type=Path,
        required=True,
        help="the CSV that was ingested, so the control is built from a row it holds",
    )
    args = ap.parse_args(argv)

    with args.reference_csv.open(newline="", encoding="utf-8") as fh:
        in_reference = [r["cdrh3_10mer"] for r in csv.DictReader(fh)]
    provided = provided_cdrh3s()

    rows = []
    # Pick reference loops that are far from trastuzumab, so the only thing
    # that can reject the control is the provided-data comparison.
    far = [s for s in in_reference if hamming(s, TZ_CDRH3_10) >= 6]
    if not far:
        raise SystemExit("no reference loop far enough from trastuzumab")

    for n, base in enumerate(far[:3]):
        # One substitution away from `base`, which IS in the reference.
        for pos in range(len(base)):
            made = False
            for aa in ALPHABET:
                if aa == base[pos]:
                    continue
                cand = base[:pos] + aa + base[pos + 1:]
                if cand in provided:
                    continue  # must be novel vs the data, except for `base`
                if hamming(cand, TZ_CDRH3_10) < 4:
                    continue
                rows.append({
                    "sequence_id": f"control_d1_{n}_{pos}",
                    "kind": "paired",
                    "heavy_sequence": vh_from_cdrh3(cand),
                    "light_sequence": TRASTUZUMAB_VL,
                    "design_description": (
                        f"POSITIVE CONTROL: 1 mutation from provided {base}, "
                        f"{hamming(cand, TZ_CDRH3_10)} from trastuzumab. "
                        "MUST be rejected when the reference is in use."
                    ),
                })
                made = True
                break
            if made:
                break

    if not rows:
        raise SystemExit("could not construct a control")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"{r['sequence_id']}: {r['design_description']}")
    print(f"\nwrote {len(rows)} control(s) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
