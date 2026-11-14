#!/usr/bin/env python3
"""Markdown summary of results/denovo/props.csv: reference nanobodies vs designs (medians per generator).

    python3 denovo/s4_props/report.py > denovo/s4_props/report.md
"""
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dn_common as dn  # noqa: E402

rows = dn.read_csv(dn.OUT / "props.csv")
num = lambda r, k: float(r[k]) if r.get(k, "") not in ("", None) else None  # noqa: E731


def med(rs, k):
    v = [num(r, k) for r in rs if num(r, k) is not None]
    return f"{st.median(v):.2f} (n={len(v)})" if v else "-"


print(f"# S4 report: {len(rows)} sequences in props.csv\n")
print("## Reference nanobodies and controls\n")
print("| id | nm_tm | tnp red/amber | PSH | PPC | my_yield | my_mdratio | liab flags |\n|---|---|---|---|---|---|---|---|")
for r in rows:
    if r["is_control"] in ("1", 1):
        print(f"| {r['design_id']} | {r.get('nm_tm') or '-'} | {r.get('tnp_n_red') or '-'}/{r.get('tnp_n_amber') or '-'} | "
              f"{r.get('tnp_psh') or '-'} | {r.get('tnp_ppc') or '-'} | {r.get('my_yield') or '-'} | {r.get('my_mdratio') or '-'} | "
              f"{r.get('liab_flags') or '-'} |")
print("\n## Designs by generator / framework (medians)\n")
print("| generator | framework | n | nm_tm | tnp_n_red | my_yield | liab_n_flags | unpaired Cys | no canonical SS |\n|---|---|---|---|---|---|---|---|---|")
groups = {}
for r in rows:
    if r["is_control"] in ("1", 1):
        continue
    groups.setdefault((r["generator"], r["framework"]), []).append(r)
for (g, f), rs in sorted(groups.items()):
    print(f"| {g} | {f} | {len(rs)} | {med(rs, 'nm_tm')} | {med(rs, 'tnp_n_red')} | {med(rs, 'my_yield')} | {med(rs, 'liab_n_flags')} | "
          f"{sum(r.get('liab_unpaired_cys') == '1' for r in rs)} | {sum(r.get('liab_no_canon_ss') == '1' for r in rs)} |")
