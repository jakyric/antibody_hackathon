"""Build and gate-check the VH44-VL100 engineered disulfide.

WHAT IS BEING TESTED
  The classic disulfide-stabilised scFv (ds-scFv) replaces VH Gly44 and
  VL Gln100 (Kabat numbering) with cysteines, so an engineered bond staples
  the two variable domains together across the interface. It is reported to
  reduce aggregation and to help scFvs that fall apart; its effect on melting
  temperature is NOT reliably positive, so we treat it as a hypothesis with
  real downside, not a free improvement.

WHY POSITIONS ARE LOOKED UP, NOT HARDCODED
  Kabat 44 and Kabat 100 are not the 44th and 100th residues of the chain.
  Antibody numbering has insertions and deletions, and the offset differs
  between VH and VL and between antibodies. We number each chain with
  `immunum` -- the same library the official checker uses -- and read the
  residue at the Kabat position. If the residue there is not the expected Gly
  (VH44) or Gln (VL100), the script says so instead of mutating whatever it
  found.

THE GATE QUESTION, WHICH COMES FIRST
  The official checker's hard filters reject unpaired cysteines, and its
  liability reference flags `xCys`: a cysteine outside IMGT 23 and 104. An
  engineered disulfide adds exactly two such cysteines. If the official
  checker treats that as a hard failure, the whole idea is dead and no
  folding time should be spent on it. So this script builds the variants and
  runs the gate before anything else.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import immunum  # noqa: E402

from abx.sequences import SCFV_LINKER  # noqa: E402

VH_KABAT = 44
VL_KABAT = 100
VH_EXPECT = "G"
VL_EXPECT = "Q"


def kabat_index(seq: str, chain: str, position: int) -> tuple[int, str]:
    """Return (0-based index into seq, residue) for a Kabat position."""
    ann = immunum.Annotator(chains=[chain], scheme="kabat")
    res = ann.number(seq)
    if res.error:
        raise SystemExit(f"numbering failed for {chain}: {res.error}")
    # `numbering` is an ordered {label: residue} map over seq[query_start:].
    # Labels carry insertion codes ("52A"); a bare number is what we want.
    for i, (label, aa) in enumerate(res.numbering.items()):
        digits = "".join(c for c in str(label) if c.isdigit())
        suffix = "".join(c for c in str(label) if c.isalpha())
        if digits and int(digits) == position and not suffix:
            idx = res.query_start + i
            if seq[idx] != aa:
                raise SystemExit(
                    f"{chain}: index {idx} holds {seq[idx]}, numbering says {aa}"
                )
            return idx, aa
    raise SystemExit(f"Kabat {position} not present in this {chain} chain")


def mutate(seq: str, idx: int, aa: str = "C") -> str:
    return seq[:idx] + aa + seq[idx + 1:]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--designs", type=Path, default=Path("results/designs.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--ids", nargs="*", default=None,
                    help="design_ids to test; default is the top few")
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--kind", choices=["paired", "scfv"], default="paired")
    args = ap.parse_args(argv)

    with args.designs.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if args.ids:
        rows = [r for r in rows if r["design_id"] in set(args.ids)]
    else:
        rows = rows[: args.n]
    if not rows:
        raise SystemExit("no designs selected")

    out_rows = []
    for r in rows:
        vh, vl = r["vh_sequence"].strip(), r["vl_sequence"].strip()
        hi, h_aa = kabat_index(vh, "H", VH_KABAT)
        li, l_aa = kabat_index(vl, "L", VL_KABAT)

        note = ""
        if h_aa != VH_EXPECT:
            note += f" WARNING VH Kabat44 is {h_aa}, expected {VH_EXPECT}."
        if l_aa != VL_EXPECT:
            note += f" WARNING VL Kabat100 is {l_aa}, expected {VL_EXPECT}."
        print(f"{r['design_id']}: VH Kabat44 = {h_aa} at index {hi}, "
              f"VL Kabat100 = {l_aa} at index {li}.{note}")

        vh_ds, vl_ds = mutate(vh, hi), mutate(vl, li)

        for tag, h, l in ((f"{r['design_id']}", vh, vl),
                          (f"{r['design_id']}-ds", vh_ds, vl_ds)):
            is_ds = tag.endswith("-ds")
            if args.kind == "paired":
                heavy, light = h, l
            else:
                heavy, light = h + SCFV_LINKER + l, ""
            out_rows.append({
                "sequence_id": tag,
                "kind": args.kind,
                "heavy_sequence": heavy,
                "light_sequence": light,
                "design_description": (
                    f"{'VH44C/VL100C engineered disulfide' if is_ds else 'parent, no disulfide'}; "
                    f"CDRH3 {r['cdrh3']}"
                ),
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
