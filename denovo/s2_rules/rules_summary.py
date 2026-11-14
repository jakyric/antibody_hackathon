"""Per-generator summary of results/denovo/rules.csv -> stdout (markdown). Used for decisions.md / the final report."""
import csv, collections, statistics as st, sys
from pathlib import Path
rows = list(csv.DictReader(open(Path(__file__).resolve().parents[2] / "results/denovo/rules.csv")))
by = collections.defaultdict(list)
for r in rows:
    by[r["generator"]].append(r)
b = lambda rs, k: sum(r[k] == "True" for r in rs)
print("| generator | n | basic_pass | canonical Cys | abguard_pass | de_novo class (VM) | identity_nearest_known < 0.70 | strict (both) | trast VH id. median [min-max] | nearest known id. median [min-max] |")
print("|---|---|---|---|---|---|---|---|---|---|")
for g, rs in sorted(by.items()):
    tr = [float(r["identity_to_trastuzumab"]) for r in rs if r["identity_to_trastuzumab"]]
    nk = [float(r["identity_nearest_known"]) for r in rs]
    print(f"| {g} | {len(rs)} | {b(rs,'basic_pass')} | {b(rs,'canon_cys_ok')} | {b(rs,'abguard_pass')} | {b(rs,'denovo_vm_pass')} | {b(rs,'ident_lt70')} | {b(rs,'denovo_strict_pass')} | "
          f"{st.median(tr):.3f} [{min(tr):.2f}-{max(tr):.2f}] | {st.median(nk):.3f} [{min(nk):.2f}-{max(nk):.2f}] |")
