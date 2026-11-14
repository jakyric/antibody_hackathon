"""Build the +/- R66G folding panel: does the stabilising mutation move the epitope?

THE QUESTION
  Every stabilised design carries one light-chain change, VL R66G, taken from
  the literature on trastuzumab scFv. The published trade is thermostability
  up, measured affinity somewhat down. What the literature does not tell us is
  whether the antibody still binds the SAME PATCH on HER2 -- and epitope
  competition with trastuzumab is one of the three judged axes.

  Kabat L66 sits in the loop just after CDR-L2. It is framework, not CDR, but
  it is not far from the paratope, so "framework, therefore harmless" is an
  assumption, not a finding. This panel tests it.

THE DESIGN
  Matched pairs. Each entry appears twice, identical except for the single
  R66G substitution in the light chain, folded by the same model with the same
  seed and the same template. The difference between a pair is then attributable
  to the mutation and nothing else.

  Trastuzumab itself is in the panel, both ways. It is the control that makes
  the rest readable: the epitope shift seen for trastuzumab +/- R66G is the
  shift the mutation causes on a backbone we KNOW binds. A design whose shift is
  no larger than trastuzumab's has not been shown to move.

  Several seeds, because a single prediction has no error bar. Two structures
  of the SAME sequence at different seeds give the noise floor, and a pair
  difference smaller than that floor means nothing.

USAGE
  python scripts/r66g_panel.py --out runs/r66g/panel.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402

# The stabilising substitution, as a sequence edit rather than a Kabat index:
# ...GVPSRFSGS[R]SGTDFTL... -> ...GVPSRFSGS[G]SGTDFTL...
WT_CONTEXT = "GVPSRFSGSRSGTDFTL"
MUT_CONTEXT = "GVPSRFSGSGSGTDFTL"


def apply_r66g(vl: str) -> str:
    if MUT_CONTEXT in vl:
        return vl  # already stabilised
    if vl.count(WT_CONTEXT) != 1:
        raise ValueError(
            f"expected exactly one {WT_CONTEXT!r} in the light chain, "
            f"found {vl.count(WT_CONTEXT)}"
        )
    return vl.replace(WT_CONTEXT, MUT_CONTEXT)


def strip_r66g(vl: str) -> str:
    if WT_CONTEXT in vl:
        return vl
    if vl.count(MUT_CONTEXT) != 1:
        raise ValueError(f"expected exactly one {MUT_CONTEXT!r} in the light chain")
    return vl.replace(MUT_CONTEXT, WT_CONTEXT)


def r66g_position(vl: str) -> int:
    """1-based index of the substituted residue, for the record."""
    return vl.index(WT_CONTEXT) + WT_CONTEXT.index("RSGTDFTL") + 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--submission", type=Path, default=Path("results/SUBMISSION.csv"))
    ap.add_argument("--out", type=Path, default=Path("runs/r66g/panel.csv"))
    ap.add_argument(
        "--designs",
        default="tz-opt-02,tz-opt-07,tz-opt-18",
        help="design ids to include, each folded both with and without R66G",
    )
    args = ap.parse_args(argv)

    wanted = args.designs.split(",")
    rows = {r["sequence_id"]: r for r in csv.DictReader(args.submission.open())}

    pos = r66g_position(TRASTUZUMAB_VL)
    print(f"R66G is residue {pos} of the light chain (1-based): "
          f"{TRASTUZUMAB_VL[pos-1]} -> G")

    panel = [
        {
            "id": "trastuzumab",
            "group": "control_wt",
            "cdrh3": "WGGDGFYAMD",
            "kd_nM": "0.1345",
            "vh_sequence": TRASTUZUMAB_VH,
            "vl_sequence": strip_r66g(TRASTUZUMAB_VL),
        },
        {
            "id": "trastuzumab-stab",
            "group": "control_stab",
            "cdrh3": "WGGDGFYAMD",
            "kd_nM": "",
            "vh_sequence": TRASTUZUMAB_VH,
            "vl_sequence": apply_r66g(TRASTUZUMAB_VL),
        },
    ]

    for did in wanted:
        if did not in rows:
            raise SystemExit(f"{did} is not in {args.submission}")
        vh = rows[did]["heavy_sequence"]
        vl = rows[did]["light_sequence"]
        cdrh3 = vh[vh.index("YYCSR") + 5 : vh.index("YWGQGT", 10)]
        panel.append({
            "id": did, "group": "design_wt", "cdrh3": cdrh3, "kd_nM": "",
            "vh_sequence": vh, "vl_sequence": strip_r66g(vl),
        })
        panel.append({
            "id": f"{did}-stab", "group": "design_stab", "cdrh3": cdrh3, "kd_nM": "",
            "vh_sequence": vh, "vl_sequence": apply_r66g(vl),
        })

    # The pairs must differ in exactly one residue, or the comparison is not
    # about R66G. Check it rather than trust it.
    for i in range(0, len(panel), 2):
        a, b = panel[i], panel[i + 1]
        diffs = [k for k, (x, y) in enumerate(zip(a["vl_sequence"], b["vl_sequence"])) if x != y]
        assert a["vh_sequence"] == b["vh_sequence"], f"{a['id']}: heavy chains differ"
        assert diffs == [pos - 1], f"{a['id']}: light chains differ at {diffs}, expected [{pos-1}]"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(panel[0]))
        w.writeheader()
        w.writerows(panel)
    print(f"{len(panel)} entries ({len(panel)//2} matched pairs) -> {args.out}")
    for r in panel:
        print(f"  {r['id']:22s} {r['group']:14s} {r['cdrh3']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
