"""Stage 2: diverse CDR-H3 shortlist and Vernier-zone panel.

Shortlist: greedy by BLOSUM joint probability (Moritz's recipe), members differ from each other at >= N positions.
Vernier: for every shortlisted CDR-H3, each Vernier residue that touches a position changed from the base loop is
varied one at a time over its BLOSUM-45 >= +1 substitutions (plus the all-WT-Vernier genotype).
Output: panel_in.csv (abguard / folding input) and panel_meta.csv.
"""
from __future__ import annotations

import csv

import pandas as pd

from .common import VERNIER, blosum, build_fv, hamming, out, save_summary


def vernier_options(min_score: int, residues: list[str]) -> dict[str, list[str]]:
    B = blosum("BLOSUM45")
    opts = {}
    for name in residues:
        wt = VERNIER[name][2]
        opts[name] = sorted((a for a in B[wt] if a not in (wt, "C", "P") and B[wt][a] >= min_score), key=lambda a: -B[wt][a])
    return opts


def run(cfg: dict) -> None:
    d = cfg["design"]
    enum = pd.read_csv(out(cfg, "enum.csv"))
    short = []
    for r in enum.itertuples():
        if all(hamming(r.cdrh3, s.cdrh3) >= d["shortlist_min_diff"] for s in short):
            short.append(r)
        if len(short) >= d["shortlist"]:
            break
    opts = vernier_options(cfg["vernier"]["min_blosum45"], cfg["vernier"]["residues"])
    print("  Vernier options:", {k: "".join(v) for k, v in opts.items()})

    rows, meta = [], []
    for r in short:
        changed = {int(m[1:-1]) for m in r.mutations_vs_base.split(";")}
        genotypes = [("WT", {})]
        for name, aas in opts.items():
            if VERNIER[name][3] & changed:                   # touches a position we changed
                genotypes += [(f"{name}{VERNIER[name][2]}>{a}", {name: a}) for a in aas]
        for tag, sub in genotypes:
            vh, vl = build_fv(r.cdrh3, sub)
            sid = f"{r.cdrh3}_{tag.replace('>', '')}"
            rows.append([sid, "paired", vh, vl, f"{r.mutations_vs_wt}|{tag}"])
            meta.append((sid, r.cdrh3, r.blosum_logp, r.ridge_sum, r.mutations_vs_wt, r.mutations_vs_base, tag))
    with open(out(cfg, "panel_in.csv"), "w", newline="") as f:      # column layout of `abguard team-validate`
        w = csv.writer(f)
        w.writerow(["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"])
        w.writerows(rows)
    pd.DataFrame(meta, columns=["id", "cdrh3", "blosum_logp", "ridge_sum", "mutations_vs_wt", "mutations_vs_base",
                                "vernier_changes"]).to_csv(out(cfg, "panel_meta.csv"), index=False)
    print(f"  {len(short)} CDR-H3 variants, {len(rows)} sequences")
    save_summary(cfg, n_shortlist=len(short), n_panel=len(rows))
