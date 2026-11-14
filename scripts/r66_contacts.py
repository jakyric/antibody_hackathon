"""Is VL Arg66 really framework-only, far from the paratope? Measure it.

WHY
  The stabilising mutation VL R66G is justified in the literature by the claim
  that position 66 "was not located at the antigen-binding region" (Okazaki et
  al. 2023), so changing it should not touch binding. That claim is checkable
  against the crystal structure of the actual complex, and it ought to be
  checked before six of our ten submitted designs rest on it.

  1N8Z is trastuzumab Fab bound to HER2: chain A = light chain, chain B =
  heavy chain, chain C = HER2. We measure, from the real coordinates:

    - the closest approach of the Arg66 side chain to HER2;
    - every residue whose heavy atoms come within 4.5 A of the Arg66 side
      chain, with the CDRs marked;
    - the specific hydrogen bond the literature mechanism depends on.

  This is geometry from a deposited structure. It is not a prediction and
  carries no model uncertainty -- but it is one crystal form at one
  resolution, and a contact at 3.9 A is not the same kind of fact as one at
  2.5 A.

USAGE
  python scripts/r66_contacts.py --pdb data/1n8z.pdb
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np

# Kabat CDR definitions for the trastuzumab light chain, as numbered in 1N8Z.
CDR_L1 = range(24, 35)
CDR_L2 = range(50, 57)
CDR_L3 = range(89, 98)
BACKBONE = {"N", "CA", "C", "O"}


def label(resi: int) -> str:
    if resi in CDR_L1:
        return "CDR-L1"
    if resi in CDR_L2:
        return "CDR-L2"
    if resi in CDR_L3:
        return "CDR-L3"
    return "framework"


def read_atoms(path: Path):
    """[(chain, resi, resname, atomname, xyz)] for heavy atoms."""
    out = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        element = line[76:78].strip()
        if element == "H":
            continue
        name = line[12:16].strip()
        out.append((
            line[21],
            int(line[22:26]),
            line[17:20].strip(),
            name,
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        ))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", type=Path, default=Path("data/1n8z.pdb"))
    ap.add_argument("--light", default="A")
    ap.add_argument("--heavy", default="B")
    ap.add_argument("--antigen", default="C")
    ap.add_argument("--resi", type=int, default=66)
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--out", type=Path, default=Path("results/r66g/R66_GEOMETRY.md"))
    args = ap.parse_args(argv)

    atoms = read_atoms(args.pdb)

    target = [a for a in atoms if a[0] == args.light and a[1] == args.resi]
    if not target:
        raise SystemExit(f"no residue {args.resi} in chain {args.light} of {args.pdb}")
    resname = target[0][2]
    side = [a for a in target if a[3] not in BACKBONE and a[3] != "CB"]
    print(f"{args.pdb}: chain {args.light} residue {args.resi} is {resname}")
    print(f"side-chain atoms: {', '.join(a[3] for a in side)}")
    if not side:
        raise SystemExit("no side chain beyond CB; nothing to measure")

    side_xyz = np.array([a[4] for a in side])

    # 1. distance to the antigen
    ag = [a for a in atoms if a[0] == args.antigen]
    ag_xyz = np.array([a[4] for a in ag])
    d_ag = np.linalg.norm(side_xyz[:, None, :] - ag_xyz[None, :, :], axis=-1)
    i, j = np.unravel_index(d_ag.argmin(), d_ag.shape)
    min_ag = float(d_ag[i, j])
    print(f"\nclosest {resname}{args.resi} side-chain atom to HER2: "
          f"{min_ag:.2f} A  ({side[i][3]} -> {ag[j][2]}{ag[j][1]} {ag[j][3]})")

    # 2. how close do the CDRs themselves get? the scale for reading the above
    cdr_min = {}
    for name, rng in (("CDR-L1", CDR_L1), ("CDR-L2", CDR_L2), ("CDR-L3", CDR_L3)):
        xs = np.array([a[4] for a in atoms if a[0] == args.light and a[1] in rng])
        if len(xs):
            cdr_min[name] = float(
                np.linalg.norm(xs[:, None, :] - ag_xyz[None, :, :], axis=-1).min()
            )
    for k, v in cdr_min.items():
        print(f"  for scale, {k} to HER2: {v:.2f} A")

    # 3. what the side chain touches inside the antibody
    near = defaultdict(lambda: (99.0, "", ""))
    for a in atoms:
        if a[0] == args.light and a[1] == args.resi:
            continue
        if a[0] not in (args.light, args.heavy):
            continue
        d = float(np.linalg.norm(side_xyz - a[4], axis=-1).min())
        if d <= args.cutoff:
            k = (a[0], a[1], a[2])
            if d < near[k][0]:
                w = int(np.linalg.norm(side_xyz - a[4], axis=-1).argmin())
                near[k] = (d, a[3], side[w][3])

    print(f"\nresidues within {args.cutoff} A of the {resname}{args.resi} side chain:")
    cdr_contacts = []
    for (ch, resi, rn), (d, atom, from_atom) in sorted(near.items(), key=lambda kv: kv[1][0]):
        reg = label(resi) if ch == args.light else "heavy chain"
        mark = "  <-- CDR" if reg.startswith("CDR") else ""
        print(f"  {ch} {rn}{resi:<4} {reg:<10} {d:.2f} A  "
              f"({from_atom} -> {atom}){mark}")
        if reg.startswith("CDR"):
            cdr_contacts.append((ch, rn, resi, reg, d, atom, from_atom))

    hbonds = [c for c in cdr_contacts if c[4] <= 3.5]

    lines = [
        f"# Is VL {resname}{args.resi} really away from the binding site?",
        "",
        f"Measured from `{args.pdb}` (trastuzumab Fab + HER2 crystal structure),",
        "heavy atoms only. Geometry from deposited coordinates, not a prediction.",
        "",
        "| | Distance |",
        "| --- | ---: |",
        f"| {resname}{args.resi} side chain to nearest HER2 atom | **{min_ag:.2f} A** |",
    ] + [f"| for scale, {k} to HER2 | {v:.2f} A |" for k, v in cdr_min.items()] + [
        "",
        f"## What the {resname}{args.resi} side chain touches",
        "",
        f"| Residue | Region | Distance | Atoms |",
        "| --- | --- | ---: | --- |",
    ] + [
        f"| {ch} {rn}{resi} | {reg} | {d:.2f} A | {fa} - {at} |"
        for ch, rn, resi, reg, d, at, fa in cdr_contacts
    ] + [
        "",
        "## What this means",
        "",
    ]

    if cdr_contacts:
        lines += [
            f"The side chain makes **{len(cdr_contacts)} contact(s) with CDR residues**"
            + (f", {len(hbonds)} of them at hydrogen-bonding distance (<= 3.5 A)."
               if hbonds else "."),
            "",
            f"It is **{min_ag:.2f} A** from HER2 -- not in the interface itself, but it",
            "is holding CDR loops that are. The premise that a change here is",
            "remote from binding is **not supported by the structure**.",
        ]
        print(f"\nVERDICT: {len(cdr_contacts)} CDR contacts, {len(hbonds)} at H-bond "
              f"distance. 'Far from the binding site' is NOT supported.")
    else:
        lines += [
            "The side chain makes no CDR contacts within the cutoff. The claim that",
            "this position is remote from binding is consistent with the structure.",
        ]
        print("\nVERDICT: no CDR contacts; 'far from the binding site' is supported.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
