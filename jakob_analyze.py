#!/usr/bin/env python3
"""Analyse trastuzumab CDRH3 FACS + SPR data. Pure stdlib."""
import csv, zipfile, io, sys, math
from collections import Counter, defaultdict

WT = "WGGDGFYAMD"
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"

# ---------- FACS ----------
rows = []
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    name = [n for n in z.namelist() if n.endswith(".csv") and not n.startswith("__MACOSX")][0]
    with z.open(name) as fh:
        r = csv.DictReader(io.TextIOWrapper(fh, "utf-8"))
        for d in r:
            rows.append((d["cdrh3"], d["binding_class"], int(d["edit_distance"])))

print(f"FACS rows: {len(rows)}  unique cdrh3: {len(set(x[0] for x in rows))}")
print("class counts:", Counter(x[1] for x in rows))

# edit distance vs class
tab = defaultdict(Counter)
for s, c, ed in rows:
    tab[ed][c] += 1
print("\n== binding class by edit distance ==")
print(f"{'ed':>3} {'n':>7} {'high%':>7} {'med%':>7} {'low%':>7}")
for ed in sorted(tab):
    c = tab[ed]; n = sum(c.values())
    print(f"{ed:>3} {n:>7} {100*c['high']/n:>7.1f} {100*c['medium']/n:>7.1f} {100*c['low']/n:>7.1f}")

# position-wise enrichment: log2( P(aa|high) / P(aa|low) )
AAs = "ACDEFGHIKLMNPQRSTVWY"
hi = [Counter() for _ in range(10)]
lo = [Counter() for _ in range(10)]
for s, c, ed in rows:
    if c == "high":
        for i, a in enumerate(s): hi[i][a] += 1
    elif c == "low":
        for i, a in enumerate(s): lo[i][a] += 1

nhi = sum(hi[0].values()); nlo = sum(lo[0].values())
print(f"\n== position-wise log2 enrichment high vs low (n_high={nhi}, n_low={nlo}) ==")
print("WT residue marked *. Only aa with >=200 total counts at that position shown.")
enrich = []
for i in range(10):
    scores = {}
    for a in AAs:
        h, l = hi[i][a], lo[i][a]
        if h + l < 200: continue
        ph = (h + 1) / (nhi + 20); pl = (l + 1) / (nlo + 20)
        scores[a] = math.log2(ph / pl)
    enrich.append(scores)
    top = sorted(scores.items(), key=lambda kv: -kv[1])
    wt = WT[i]
    wtmark = f"WT={wt}({scores.get(wt, float('nan')):+.2f})"
    best = "  ".join(f"{a}{'*' if a==wt else ''}:{v:+.2f}" for a, v in top[:8])
    worst = "  ".join(f"{a}{'*' if a==wt else ''}:{v:+.2f}" for a, v in top[-4:])
    print(f"\npos {i+1} (H3 res {99+i}) {wtmark}\n   best : {best}\n   worst: {worst}")

# ---------- SPR ----------
print("\n\n================ SPR ================")
spr = []
with open(f"{BASE}/trastuzumab_her2_spr_all_768.csv") as fh:
    for d in csv.DictReader(fh):
        spr.append(d)
print("SPR rows:", len(spr), Counter(d["binding_result"] for d in spr))

TZ_VH_WT = None
for d in spr:
    if d["cdrh3"] == WT: TZ_VH_WT = d["vh_sequence"]

def ed(a, b):
    if len(a) != len(b): return None
    return sum(x != y for x, y in zip(a, b))

quant = [d for d in spr if d["binding_result"] == "quantified_binding" and d["kd_nM"]]
print("quantified:", len(quant))
# drop the other antibody
quant = [d for d in quant if d["cdrh3"] != "ARDRIVDDYYYYGMDV"]
quant.sort(key=lambda d: float(d["kd_nM"]))
print(f"\n== all quantified variants sorted by KD ==")
print(f"{'cdrh3':<18}{'len':>4}{'ed_wt':>6}{'KD_nM':>12}  vh_mut_outside_cdrh3")
for d in quant:
    s = d["cdrh3"]; e = ed(s, WT)
    # check VH differences outside the CDRH3 region
    vh = d["vh_sequence"]
    extra = ""
    if TZ_VH_WT and len(vh) == len(TZ_VH_WT):
        diffs = [(i, TZ_VH_WT[i], vh[i]) for i in range(len(vh)) if vh[i] != TZ_VH_WT[i]]
        # cdrh3 sits at index 98..107 in the 120-mer
        outside = [f"{o}{i+1}{n}" for i, o, n in diffs if not (98 <= i <= 107)]
        extra = ",".join(outside)
    print(f"{s:<18}{len(s):>4}{str(e):>6}{float(d['kd_nM']):>12.4f}  {extra}")
