#!/usr/bin/env python3
"""Create a clean framework file from the organizers' pre-annotated VHH file.

Data/03-vhh/pmvk_caplacizumab.cif.zst is built to be run as it is: it
has its own target (PMVK, chain A), placeholders (UNK) where the CDRs go, and
the conditions baked in as separate columns. To use the framework against a
DIFFERENT target, all of this has to go, otherwise parsing crashes on the
placeholders
("Unresolved link valence at H/23/13-13/UNK").

This writes a minimal mmCIF with only chain H, only real residues, and only the
common _atom_site columns.

    python scripts/vhh_framework_cif.py --input vhh.cif --chain H --out framework.cif
"""

from __future__ import annotations

import argparse
from pathlib import Path

HEADER = """data_{name}
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""

SKIP_COMP = {"UNK", "X"}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--chain", default="H")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="vhh_framework")
    args = ap.parse_args(argv)

    cols: list[str] = []
    rows = []
    for line in args.input.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
            continue
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        if len(f) < len(cols):
            continue
        rec = dict(zip(cols, f))
        chain = rec.get("auth_asym_id") or rec.get("label_asym_id")
        if chain != args.chain:
            continue
        comp = (rec.get("label_comp_id") or "").upper()
        if comp in SKIP_COMP or rec.get("type_symbol") == "X":
            continue          # the placeholders where the CDRs are to be designed
        rows.append(rec)

    if not rows:
        raise SystemExit(f"found no residues in chain {args.chain}")

    out = [HEADER.format(name=args.name)]
    for i, rec in enumerate(rows, 1):
        seq = rec.get("auth_seq_id") or rec.get("label_seq_id")
        out.append(" ".join([
            "ATOM", str(i), rec["type_symbol"], rec["label_atom_id"], ".",
            rec["label_comp_id"], args.chain, "1", str(seq), "?",
            rec["Cartn_x"], rec["Cartn_y"], rec["Cartn_z"],
            rec.get("occupancy", "1.0"), rec.get("B_iso_or_equiv", "0.0"), "?",
            str(seq), rec["label_comp_id"], args.chain, rec["label_atom_id"], "1",
        ]))
    out.append("#")
    args.out.write_text("\n".join(out) + "\n")

    resnums = sorted({int(r.get("auth_seq_id") or r["label_seq_id"]) for r in rows})
    runs = []
    start = prev = resnums[0]
    for n in resnums[1:]:
        if n != prev + 1:
            runs.append((start, prev))
            start = n
        prev = n
    runs.append((start, prev))
    print(f"{len(rows)} atoms, {len(resnums)} residues -> {args.out}")
    print("  framework regions: " + ", ".join(f"{a}-{b}" for a, b in runs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
