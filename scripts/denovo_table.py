#!/usr/bin/env python3
"""Write the submission table for the de novo track, with provenance in the description.

The official checker wants exactly five columns: sequence_id, kind, heavy_sequence,
light_sequence, design_description. Our designs are Fvs, and are delivered as
scFv (VH + (G4S)x3 + VL in heavy_sequence, empty light_sequence) — one of the three
legal formats in the de novo track.

design_description carries the traceability: model, weight, seed, run id, and the
novelty numbers the selection rests on. All in English; the table is public.

    python scripts/denovo_table.py --designs results/denovo/candidates_pass.csv \
        --out results/denovo/denovo_submission.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.sequences import SCFV_LINKER, TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

COLUMNS = ["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"]
MODEL = "RFD4-Proteina v0.3.0"


def describe(rec: dict) -> str:
    kind = (rec.get("kind") or "").strip()
    if not kind:
        kind = "scfv" if rec.get("vl_sequence", "").strip() else "vhh"
    what = ("de novo VHH; three CDRs designed on a nanobody framework against HER2 domain IV "
            "at the trastuzumab epitope"
            if kind == "vhh" else
            "all six CDRs designed on the trastuzumab Fv framework against HER2 domain IV")
    bits = [
        what,
        f"{MODEL} weight {str(rec.get('weight', '?')).removeprefix('vhh-')}",
        f"seed {rec.get('seed', '?')}",
        f"run {rec.get('run_id', '?')}",
        f"CDRH3 {rec.get('cdrh3', '?')}",
        f"CDR identity to trastuzumab {rec.get('identity_cdr_pct', '?')} %",
        f"nearest provided CDRH3 {rec.get('cdrh3_min_dist_provided', '?')} mutations away",
    ]
    for key, label in (("fold_note", None), ("rank_note", None)):
        if rec.get(key):
            bits.append(rec[key])
    return "; ".join(bits)


def rows_from(path: Path, prefix: str, limit: int | None) -> list[dict]:
    out: list[dict] = []
    with path.open(newline="", encoding="utf-8") as fh:
        records = list(csv.DictReader(fh))
    if "rank" in (records[0] if records else {}):
        records.sort(key=lambda r: float(r["rank"]))
    if limit:
        records = records[:limit]
    for i, rec in enumerate(records, 1):
        kind = (rec.get("kind") or "").strip()
        if not kind:
            # An empty light chain means a VHH. Without this, a VHH gets the
            # scFv linker left dangling after it, and the sequence ends up wrong.
            kind = "scfv" if rec.get("vl_sequence", "").strip() else "vhh"
        if kind == "vhh":
            # one chain, no linker, no light chain
            heavy = rec["vh_sequence"].strip()
        else:
            heavy = rec["vh_sequence"].strip() + SCFV_LINKER + rec["vl_sequence"].strip()
        out.append({
            "sequence_id": f"{prefix}-{i:02d}",
            "kind": kind,
            "heavy_sequence": heavy,
            "light_sequence": "",
            "design_description": describe(rec),
        })
    return out


def reference_row() -> dict:
    return {
        "sequence_id": "tz_reference",
        "kind": "scfv",
        "heavy_sequence": TRASTUZUMAB_VH + SCFV_LINKER + TRASTUZUMAB_VL,
        "light_sequence": "",
        "design_description": "trastuzumab itself: a control that MUST fail the novelty gate",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prefix", default="tz-denovo")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--with-reference", action="store_true",
                    help="put trastuzumab first, as a control that MUST fail")
    args = ap.parse_args(argv)

    rows = rows_from(args.designs, args.prefix, args.limit or None)
    if args.with_reference:
        rows.insert(0, reference_row())

    too_long = [r for r in rows if len(r["heavy_sequence"]) > 250]
    if too_long:
        raise SystemExit(
            f"{len(too_long)} row(s) over 250 residues, e.g. {too_long[0]['sequence_id']} "
            f"({len(too_long[0]['heavy_sequence'])})"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
