"""Markdown summary of the rfd4 campaign from designs_rfd4jf.csv + geom_rfd4jf.csv (+ S3 scores when present)."""
import sys, statistics as st
from collections import defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import dn_common as dn

rows = dn.read_csv(dn.OUT / "designs_rfd4jf.csv") if (dn.OUT / "designs_rfd4jf.csv").exists() else []
geom = {r["design_id"]: r for r in dn.read_csv(dn.OUT / "geom_rfd4jf.csv")} if (dn.OUT / "geom_rfd4jf.csv").exists() else {}
f = lambda x: float(x) if x not in ("", None) else None
print(f"designs: {len(rows)}; geometry rows: {len(geom)}; unique vhh_seq: {len({r['vhh_seq'] for r in rows})}; unique CDR3: {len({r['cdr3'] for r in rows})}")
cells = defaultdict(list)
for r in rows:
    g = geom.get(r["design_id"])
    if g:
        cells[(g["framework"], g["hotspot_set"], g["c3_range"], g["pose"])].append(g)
print("\n| framework | hotspots | CDR3 | pose | n | contact pairs | CDR3 res in contact | epitope precision | hotspots touched | clash>2 pairs |")
print("|---|---|---|---|---|---|---|---|---|---|")
m = lambda xs: f"{st.mean(xs):.1f}" if xs else "-"
for k in sorted(cells):
    g = cells[k]
    ht = [f(x["hot_touched"]) / max(1, f(x["hot_total"])) for x in g]
    print(f"| {' | '.join(k)} | {len(g)} | {m([f(x['contact_pairs']) for x in g])} | {m([f(x['cdr3_res_contact']) for x in g])} | "
          f"{st.mean([f(x['epi_precision']) for x in g]):.2f} | {st.mean(ht):.2f} | {sum(1 for x in g if f(x['clash_lt2p2']) > 2)} |")
# identity of each design to its parent framework sequence (bcii10: full parent from 3EAK; capla: only the kept framework is known,
# so the fraction of the design that is inherited framework is an upper bound on identity apart from chance matches in the CDRs)
import json


def nw_identity(a, b):
    """global alignment (match 1, mismatch 0, gap -1); identity = matches / max(len)"""
    n, m = len(a), len(b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1): S[i][0] = -i
    for j in range(1, m + 1): S[0][j] = -j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            S[i][j] = max(S[i - 1][j - 1] + (a[i - 1] == b[j - 1]), S[i - 1][j] - 1, S[i][j - 1] - 1)
    i, j, mt = n, m, 0
    while i > 0 and j > 0:
        if S[i][j] == S[i - 1][j - 1] + (a[i - 1] == b[j - 1]):
            mt += a[i - 1] == b[j - 1]; i -= 1; j -= 1
        elif S[i][j] == S[i - 1][j] - 1: i -= 1
        else: j -= 1
    return mt / max(n, m)


fws = {x: json.load(open(HERE / "frameworks" / f"{x}.json")) for x in ("bcii10", "capla")}
parent = {"bcii10": "".join(open(HERE / "frameworks" / "bcii10_clean.cif").read() and [])}
import prep_inputs as PI
res = PI.read_atoms(HERE / "frameworks" / "3EAK.cif", "A")
parent["bcii10"] = "".join(PI.AA3[r["aa3"]] for r in res)
ids = defaultdict(list); kept = defaultdict(list)
for r in rows:
    fr = fws[r["framework"]]["fr"]
    kept[r["framework"]].append(sum(map(len, fr)) / len(r["vhh_seq"]))
    if r["framework"] == "bcii10":
        ids[r["framework"]].append(nw_identity(r["vhh_seq"], parent["bcii10"]))
for k in kept:
    msg = f"{k}: inherited-framework fraction {st.mean(kept[k]):.2f} (min {min(kept[k]):.2f}, max {max(kept[k]):.2f})"
    if ids[k]:
        msg += f"; identity to parent 3EAK {st.mean(ids[k]):.2f} (min {min(ids[k]):.2f}, max {max(ids[k]):.2f}, n>=0.70: {sum(1 for x in ids[k] if x >= 0.70)})"
    print(msg)
lens = defaultdict(int)
for r in rows:
    lens[len(r["cdr3"])] += 1
print("\nCDR3 length histogram:", dict(sorted(lens.items())))
