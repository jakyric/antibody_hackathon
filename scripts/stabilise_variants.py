"""Apply a framework stabilising mutation and gate-check it.

THE MUTATION
  VL Arg66Gly (Kabat), from Okazaki et al. 2023, ACS Omega 8(25):22945-22954,
  PMID 37396255, DOI 10.1021/acsomega.3c01948.

  Measured by DSF on a trastuzumab scFv:

      wild type   61.8 +/- 0.2 C
      VL-R66G     66.7 +/- 0.2 C      dTm +4.9 C, binding unchanged

  This is the best-evidenced stabilising change available to us, and the
  reason is narrow and specific: it was measured on OUR scaffold (trastuzumab
  scFv), on OUR judged axis (melting temperature), by direct measurement.
  Nothing else we considered is all three.

WHY NOT THE DISULFIDE
  The VH44-VL100 engineered disulfide was the obvious candidate and it is
  dead twice over. The gate rejects it outright (see results/disulfide/), and
  the literature does not support it anyway: Weatherill et al. 2012, PEDS
  25(7):321-329, tested it across six scFvs and found "no or little increase
  in thermal stability". The nearest 4D5-relevant measurement, Worn &
  Pluckthun 1999, bought ~0.4 M GdnHCl of chemical stability for a 29-fold
  collapse in expression yield.

WHY THIS ONE IS STILL A BET
  One paper, one mutation, one antibody. R66 sits in the light-chain
  framework, away from CDRH3, so it does not touch our novelty margin and
  adds no cysteine -- but a single published dTm is not a guarantee, and we
  have not measured it ourselves. Treat a design carrying it as a stabilised
  HYPOTHESIS, paired with its unmutated parent so the lab result can tell
  the two apart.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import immunum  # noqa: E402

from abx.sequences import SCFV_LINKER  # noqa: E402

VL_KABAT = 66
VL_EXPECT = "R"
VL_TO = "G"


def kabat_index(seq: str, chain: str, position: int) -> tuple[int, str]:
    ann = immunum.Annotator(chains=[chain], scheme="kabat")
    res = ann.number(seq)
    if res.error:
        raise SystemExit(f"numbering failed for {chain}: {res.error}")
    for i, (label, aa) in enumerate(res.numbering.items()):
        digits = "".join(c for c in str(label) if c.isdigit())
        suffix = "".join(c for c in str(label) if c.isalpha())
        if digits and int(digits) == position and not suffix:
            idx = res.query_start + i
            if seq[idx] != aa:
                raise SystemExit(
                    f"{chain}: index {idx} holds {seq[idx]}, numbering says {aa}")
            return idx, aa
    raise SystemExit(f"Kabat {position} not present in this {chain} chain")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ids", nargs="*", default=None)
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--kind", choices=["paired", "scfv"], default="paired")
    args = ap.parse_args(argv)

    with args.designs.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if args.ids:
        order = {d: i for i, d in enumerate(args.ids)}
        rows = sorted((r for r in rows if r["design_id"] in order),
                      key=lambda r: order[r["design_id"]])
    else:
        rows = rows[: args.n]
    if not rows:
        raise SystemExit("no designs selected")

    out_rows = []
    for r in rows:
        vh = r["vh_sequence"].strip()
        vl = r["vl_sequence"].strip()
        idx, aa = kabat_index(vl, "L", VL_KABAT)
        if aa != VL_EXPECT:
            raise SystemExit(
                f"{r['design_id']}: VL Kabat{VL_KABAT} is {aa}, expected {VL_EXPECT}")
        vl_mut = vl[:idx] + VL_TO + vl[idx + 1:]
        print(f"{r['design_id']}: VL Kabat{VL_KABAT} {aa}->{VL_TO} at index {idx}")

        for tag, h, l, desc in (
            (r["design_id"], vh, vl, "parent, unmodified framework"),
            (f"{r['design_id']}-stab", vh, vl_mut,
             f"VL R66G (Okazaki 2023, dTm +4.9 C on trastuzumab scFv)"),
        ):
            if args.kind == "paired":
                heavy, light = h, l
            else:
                heavy, light = h + SCFV_LINKER + l, ""
            out_rows.append({
                "sequence_id": tag,
                "kind": args.kind,
                "heavy_sequence": heavy,
                "light_sequence": light,
                "design_description": f"{desc}; CDRH3 {r['cdrh3']}",
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nwrote {len(out_rows)} rows ({len(rows)} pairs) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
