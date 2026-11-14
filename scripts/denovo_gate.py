#!/usr/bin/env python3
"""What kind of novelty do the de novo designs really have — and under which reading?

The de novo requirement is "under 70 % identity to a known antibody". Computed
over the whole sequence, nothing built on a real framework passes. Computed
over the CDRs, a design with six new loops passes easily. The script computes
both figures against trastuzumab, so the question can be asked with numbers in
hand.

    python scripts/denovo_gate.py --design runs/<run>/design
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10, cdrh3_from_vh, provided_cdrh3s  # noqa: E402
from abx.gate import ProvidedIndex, count_unpaired_cysteines, identity  # noqa: E402
from abx.sequences import SCFV_LINKER, TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

# CDR boundaries in the 1N8Z numbering, 1-based and inclusive. Same boundaries as
# in configs/her2_all6cdr.yaml; if they change in one place, they must change in both.
CDR_LIGHT = [(24, 34), (50, 56), (89, 97)]
CDR_HEAVY = [(31, 35), (50, 66), (99, 108)]


def cdrs(seq: str, spans) -> str:
    return "".join(seq[a - 1:b] for a, b in spans)


def parse_sequences(field: str) -> dict[str, str]:
    """"A:DIQ...;B:EVQ...;C:..." -> {'A': 'DIQ...', ...}"""
    out = {}
    for part in field.split(";"):
        if ":" in part:
            cid, seq = part.split(":", 1)
            out[cid.strip()] = seq.strip()
    return out


def evaluate(vl: str, vh: str, novelty: ProvidedIndex | None) -> dict:
    scfv = vh + SCFV_LINKER + vl
    row = {
        "vh_len": len(vh),
        "vl_len": len(vl),
        "scfv_len": len(scfv),
        "unpaired_cys": count_unpaired_cysteines(scfv),
    }
    # the whole sequence against trastuzumab
    same_len = len(vh) == len(TRASTUZUMAB_VH) and len(vl) == len(TRASTUZUMAB_VL)
    row["identity_whole_pct"] = (
        round(100 * identity(vh + vl, TRASTUZUMAB_VH + TRASTUZUMAB_VL), 1)
        if same_len else ""
    )
    # only the CDRs
    if same_len:
        mine = cdrs(vl, CDR_LIGHT) + cdrs(vh, CDR_HEAVY)
        theirs = cdrs(TRASTUZUMAB_VL, CDR_LIGHT) + cdrs(TRASTUZUMAB_VH, CDR_HEAVY)
        row["identity_cdr_pct"] = round(100 * identity(mine, theirs), 1)
    else:
        row["identity_cdr_pct"] = ""
    # CDR-H3 against what was provided (the optimization track's requirement B)
    h3 = cdrh3_from_vh(vh)
    row["cdrh3"] = h3 or ""
    if h3 and novelty is not None:
        row["cdrh3_min_dist_provided"] = novelty.min_distance(h3)
        row["cdrh3_identity_tz_pct"] = round(100 * identity(h3, TZ_CDRH3_10), 1)
    else:
        row["cdrh3_min_dist_provided"] = ""
        row["cdrh3_identity_tz_pct"] = ""
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", type=Path, required=True,
                    help="the directory that module-run design wrote to")
    ap.add_argument("--out", type=Path, default=Path("results/denovo"))
    ap.add_argument("--no-novelty", action="store_true",
                    help="skip the novelty index (it reads the whole FACS table)")
    args = ap.parse_args()

    metrics = args.design / "generation" / "inference_metrics" / "per_example_metrics.csv"
    if not metrics.exists():
        raise SystemExit(f"could not find {metrics}")

    novelty = None if args.no_novelty else ProvidedIndex(provided_cdrh3s())

    rows = []
    with metrics.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            chains = parse_sequences(rec.get("sequences", ""))
            vl, vh = chains.get("A"), chains.get("B")
            if not vl or not vh:
                continue
            rows.append({"design_id": rec.get("example_id", ""),
                         **evaluate(vl, vh, novelty),
                         "hotspot_obedience": rec.get("hotspot_obedience.hotspot_obedience", ""),
                         "pre_oracle_pass": rec.get("pre_oracle_pass", "")})

    if not rows:
        raise SystemExit("no designs with both light and heavy chain in the table")

    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / f"{args.design.parent.name or 'design'}_gate.csv"
    with dest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def frac(key, test):
        vals = [r[key] for r in rows if r[key] != ""]
        return f"{sum(1 for v in vals if test(v))}/{len(vals)}" if vals else "-"

    print(f"{len(rows)} designs -> {dest}")
    print(f"  identity to trastuzumab, WHOLE sequence: "
          f"median {sorted(r['identity_whole_pct'] for r in rows)[len(rows)//2]} %")
    print(f"  identity to trastuzumab, ONLY the CDRs:  "
          f"median {sorted(r['identity_cdr_pct'] for r in rows)[len(rows)//2]} %")
    print(f"  below 70 % on the whole sequence: {frac('identity_whole_pct', lambda v: v < 70)}")
    print(f"  below 70 % on the CDRs:           {frac('identity_cdr_pct', lambda v: v < 70)}")
    print(f"  without unpaired cysteines:       {frac('unpaired_cys', lambda v: v == 0)}")
    if novelty is not None:
        print(f"  CDR-H3 at least 3 from all provided data: "
              f"{frac('cdrh3_min_dist_provided', lambda v: v >= 3)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
