"""Build the submission table for the official checker from our ranked designs.

WHY THIS EXISTS
  `results/designs.csv` is our working table: it carries the ranking, the
  reasoning and the intermediate numbers. The organisers' validator wants
  exactly five columns and nothing else. Keeping the two apart means we can
  change our own bookkeeping without breaking the submission format.

FORMAT
  The official checker splits rows on `kind` (vhh / scfv / paired) and validates each
  partition against its own schema. Trastuzumab is a paired antibody, so the
  optimisation designs are emitted as `paired` by default (heavy = VH,
  light = VL). `--kind scfv` emits the single-chain form instead, with the
  linker-joined sequence in heavy_sequence and an empty light_sequence, so we
  can check both against the real validator rather than guessing which one it
  wants.

  The trastuzumab reference row is included on purpose when --with-reference
  is passed: it MUST fail the novelty gate. A gate run where trastuzumab
  passes is a broken gate, not a good result.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

TZ_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TZ_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
SCFV_LINKER = "GGGGSGGGGSGGGGS"

COLUMNS = [
    "sequence_id",
    "kind",
    "heavy_sequence",
    "light_sequence",
    "design_description",
]


def rows_from_designs(path: Path, kind: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            vh = rec["vh_sequence"].strip()
            vl = rec["vl_sequence"].strip()
            if kind == "paired":
                heavy, light = vh, vl
            elif kind == "scfv":
                heavy, light = rec["scfv_sequence"].strip(), ""
            else:
                raise SystemExit(f"unknown kind: {kind}")

            # Keep the description short and factual; it is read by a human.
            desc = (
                f"optimisation; CDRH3 {rec['cdrh3']}; "
                f"{rec['n_mutations_from_trastuzumab']} mutations from trastuzumab; "
                f"nearest provided CDRH3 {rec['min_mutations_from_provided']} mutations away"
            )
            out.append(
                {
                    "sequence_id": rec["design_id"],
                    "kind": kind,
                    "heavy_sequence": heavy,
                    "light_sequence": light,
                    "design_description": desc,
                }
            )
    return out


def reference_row(kind: str) -> dict[str, str]:
    if kind == "paired":
        heavy, light = TZ_VH, TZ_VL
    else:
        heavy, light = TZ_VH + SCFV_LINKER + TZ_VL, ""
    return {
        "sequence_id": "tz_reference",
        "kind": kind,
        "heavy_sequence": heavy,
        "light_sequence": light,
        "design_description": "trastuzumab itself: a control that MUST fail the novelty gate",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--kind", choices=["paired", "scfv"], default="paired")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--with-reference",
        action="store_true",
        help="prepend trastuzumab itself as a negative control",
    )
    args = ap.parse_args(argv)

    rows = rows_from_designs(args.designs, args.kind)
    if args.with_reference:
        rows.insert(0, reference_row(args.kind))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {len(rows)} rows ({args.kind}) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
