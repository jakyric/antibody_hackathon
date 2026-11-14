"""Build a submission table from the liability-clean pool.

The shape is the one we argued for and are not changing: the stabilising
light-chain mutation on the top loops, because only a stabilised design can
improve all three judged properties, and the unmodified parent of the first
few kept as controls so the lab result can say what the mutation did.

What changes is where the loops come from. The old list was sampled and then
inspected for liabilities; this one is drawn from the complete enumeration of
loops that carry none of trastuzumab's own CDR-H3 motifs and add none.

USAGE
  python scripts/build_final_clean.py --out results/pool/SUBMISSION_clean.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path

import immunum

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.gate import levenshtein  # noqa: E402

VL_KABAT = 66
STAB_NOTE = ("VL R66G, measured dTm +4.9 C on trastuzumab scFv by DSF "
             "(Okazaki et al. 2023, ACS Omega 8:22945, PMID 37396255), "
             "binding unchanged")


def kabat_index(seq: str, chain: str, position: int) -> int:
    ann = immunum.Annotator(chains=[chain], scheme="kabat")
    res = ann.number(seq)
    if res.error:
        raise SystemExit(f"numbering failed: {res.error}")
    for i, (label, aa) in enumerate(res.numbering.items()):
        digits = "".join(c for c in str(label) if c.isdigit())
        suffix = "".join(c for c in str(label) if c.isalpha())
        if digits and int(digits) == position and not suffix:
            return res.query_start + i
    raise SystemExit(f"Kabat {position} not found")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selected", type=Path,
                    default=Path("results/pool/selected.csv"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-stab", type=int, default=6)
    ap.add_argument("--n-control", type=int, default=4)
    ap.add_argument("--kind", choices=["paired", "scfv"], default="paired")
    args = ap.parse_args(argv)

    with args.selected.open(newline="", encoding="utf-8") as fh:
        designs = list(csv.DictReader(fh))
    if len(designs) < args.n_stab:
        raise SystemExit(f"need {args.n_stab} designs, have {len(designs)}")
    designs = designs[: args.n_stab]

    plan = [(d, True) for d in designs]
    plan += [(d, False) for d in designs[: args.n_control]]

    rows = []
    for rank, (d, stab) in enumerate(plan, 1):
        vh, vl = d["vh_sequence"].strip(), d["vl_sequence"].strip()
        if stab:
            i = kabat_index(vl, "L", VL_KABAT)
            if vl[i] != "R":
                raise SystemExit(f"{d['design_id']}: VL Kabat66 is {vl[i]}, not R")
            vl = vl[:i] + "G" + vl[i + 1:]
        sid = f"{d['design_id']}-stab" if stab else d["design_id"]
        if stab:
            reason = ("drops all three of trastuzumab's own CDR-H3 liability "
                      "motifs (DG isomerisation, Met oxidation, Trp oxidation) "
                      "and introduces none; stabilised, so all three judged "
                      "properties are in play")
        else:
            reason = (f"CONTROL: {d['design_id']} without the stabilising "
                      f"mutation, so the lab result isolates what VL R66G did")
        desc = (f"{reason}. CDRH3 {d['cdrh3']}, "
                f"{d['n_mutations_from_trastuzumab']} mutations from "
                f"trastuzumab, {d['min_edit_distance_to_provided']} edits from "
                f"the nearest provided sequence.")
        if stab:
            desc += f" {STAB_NOTE}."
        rows.append({
            "sequence_id": sid,
            "kind": args.kind,
            "heavy_sequence": vh,
            "light_sequence": vl,
            "design_description": desc,
            "_rank": rank,
            "_cdrh3": d["cdrh3"],
            "_stab": "yes" if stab else "no",
        })

    loops = sorted({r["_cdrh3"] for r in rows})
    print(f"{len(rows)} candidates over {len(loops)} distinct CDRH3 loops")
    print(f"stabilised: {sum(r['_stab'] == 'yes' for r in rows)}, "
          f"controls: {sum(r['_stab'] == 'no' for r in rows)}\n")
    for r in rows:
        print(f"{r['_rank']:>2}. {r['sequence_id']:<20} {r['_cdrh3']:<12} "
              f"stab={r['_stab']}")
    d = [levenshtein(a, b) for a, b in combinations(loops, 2)]
    print(f"\npairwise edit distance between the {len(loops)} loops: "
          f"min {min(d)}, median {sorted(d)[len(d) // 2]}, max {max(d)}")

    cols = ["sequence_id", "kind", "heavy_sequence", "light_sequence",
            "design_description"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r[c] for c in cols})
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
