"""Compare each design's developability liabilities against trastuzumab's.

WHY THIS EXISTS
  The official checker's liability report gives a submission-wide total: 166 occurrences
  across 21 designs. That number is useless on its own, because trastuzumab
  itself carries most of them. Trastuzumab is a licensed drug with a DG
  isomerisation motif in CDRH3 and a tryptophan sitting right in the
  paratope. Liabilities are not disqualifying; they are a risk budget.

  The question that matters is the delta: did OUR four mutations add
  liabilities that trastuzumab did not have, or remove some it did?

  A design that adds a high-severity motif to CDRH3 is a worse expression
  bet than the baseline, and expression failure means no tier at all --
  not a B, nothing. A design that removes one is a developability argument
  we can actually make to the judges.

WHAT IT DOES NOT DO
  These are sequence motifs, not measurements. A motif flags a residue
  pattern known to be liable to a chemical modification; whether it
  actually degrades depends on local structure and solvent exposure. Treat
  a delta of +1 high-severity motif as a reason to look, not a verdict.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

SEVERITIES = ("high", "medium", "low")
REFERENCE_ID = "tz_reference"


def load(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def summarise(design: dict) -> dict:
    """Counts for one design: by severity, and for heavy CDR3 specifically."""
    liabs = design.get("liabilities") or []
    by_sev = Counter(l["severity"] for l in liabs)
    # The loop we actually redesigned. Everything else is inherited.
    h3 = [
        l for l in liabs
        if l.get("chain") == "H" and l.get("region") == "cdr3"
    ]
    return {
        "sequence_id": design["sequence_id"],
        "n_total": len(liabs),
        "n_high": by_sev["high"],
        "n_medium": by_sev["medium"],
        "n_low": by_sev["low"],
        "n_h3": len(h3),
        "h3_tags": "|".join(sorted(l["tag"] for l in h3)),
        "h3_motifs": "|".join(sorted(f"{l['tag']}:{l['motif']}" for l in h3)),
        # Multiset of CDRH3 tags, so "two tryptophans" differs from "one".
        "_h3_counter": Counter(l["tag"] for l in h3),
        "_h3_high": Counter(l["tag"] for l in h3 if l["severity"] == "high"),
    }


def diff_tags(design: Counter, ref: Counter) -> tuple[str, str]:
    """Which motif tags this design gained and lost relative to trastuzumab."""
    gained = design - ref
    lost = ref - design
    fmt = lambda c: "|".join(
        f"{t}x{n}" if n > 1 else t for t, n in sorted(c.items())
    )
    return fmt(gained), fmt(lost)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--liabilities", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    designs = load(args.liabilities)
    rows = [summarise(d) for d in designs]
    by_id = {r["sequence_id"]: r for r in rows}

    ref = by_id.get(REFERENCE_ID)
    if ref is None:
        raise SystemExit(
            f"no {REFERENCE_ID} row: rebuild the table with --with-reference, "
            "otherwise there is nothing to compare against"
        )

    for r in rows:
        for key in ("n_total", "n_high", "n_medium", "n_low", "n_h3"):
            r[f"d_{key}"] = r[key] - ref[key]
        r["h3_gained"], r["h3_lost"] = diff_tags(r["_h3_counter"], ref["_h3_counter"])
        r["h3_gained_high"], r["h3_lost_high"] = diff_tags(
            r["_h3_high"], ref["_h3_high"]
        )

    for r in rows:
        del r["_h3_counter"]
        del r["_h3_high"]

    fields = list(rows[0].keys())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["sequence_id"] != REFERENCE_ID,
                                                r["d_n_high"], r["d_n_total"])))

    # ---- readable summary -------------------------------------------------
    print(f"trastuzumab baseline: {ref['n_total']} liabilities "
          f"({ref['n_high']} high, {ref['n_medium']} medium, {ref['n_low']} low); "
          f"{ref['n_h3']} in CDRH3 [{ref['h3_motifs']}]")
    print()

    designs_only = [r for r in rows if r["sequence_id"] != REFERENCE_ID]
    better = [r for r in designs_only if r["d_n_high"] < 0]
    same = [r for r in designs_only if r["d_n_high"] == 0]
    worse = [r for r in designs_only if r["d_n_high"] > 0]

    print(f"high-severity vs trastuzumab: "
          f"{len(better)} better, {len(same)} same, {len(worse)} worse")
    print()
    print(f"{'design':<12} {'tot':>4} {'hi':>3} {'d_tot':>6} {'d_hi':>5}  "
          f"{'gained in CDRH3':<22} {'lost in CDRH3'}")
    for r in sorted(designs_only, key=lambda r: (r["d_n_high"], r["d_n_total"])):
        print(f"{r['sequence_id']:<12} {r['n_total']:>4} {r['n_high']:>3} "
              f"{r['d_n_total']:>+6} {r['d_n_high']:>+5}  "
              f"{r['h3_gained'] or '-':<22} {r['h3_lost'] or '-'}")

    # High-severity motifs trastuzumab does not carry are the real warning.
    introduced = [r for r in designs_only if r["h3_gained_high"]]
    if introduced:
        print()
        print("designs introducing a HIGH-severity motif trastuzumab lacks:")
        for r in sorted(introduced, key=lambda r: r["sequence_id"]):
            print(f"  {r['sequence_id']:<12} +{r['h3_gained_high']}")

    print()
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
