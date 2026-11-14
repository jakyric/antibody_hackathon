#!/usr/bin/env python
"""Trastuzumab CDRH3 optimisation: data landscape.

Part 1: what is actually in the two tables, and the single most important
question for the task - does FACS binding class carry affinity information?
"""
import csv, zipfile, io, math, json
from collections import Counter, defaultdict

import numpy as np

WT = "WGGDGFYAMD"
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"
OTHER_AB = "ARDRIVDDYYYYGMDV"

# ---------------------------------------------------------------- load FACS
facs = []
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    name = [n for n in z.namelist() if n.endswith(".csv") and "MACOS" not in n][0]
    with z.open(name) as fh:
        for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
            facs.append((d["cdrh3"], d["binding_class"], int(d["edit_distance"])))

facs_cls = {s: c for s, c, _ in facs}
print(f"FACS: {len(facs)} rows, {len(facs_cls)} unique CDRH3")
print("  classes:", dict(Counter(c for _, c, _ in facs)))
print("  lengths:", dict(Counter(len(s) for s, _, _ in facs)))

# which positions actually vary in the library?
print("\n== FACS library design: per-position amino-acid diversity ==")
print(f"{'pos':>4} {'WT':>3} {'n_aa':>5}  {'WT_freq':>8}  top residues")
for i in range(10):
    c = Counter(s[i] for s, _, _ in facs)
    tot = sum(c.values())
    top = " ".join(f"{a}:{100*n/tot:.0f}%" for a, n in c.most_common(6))
    print(f"{i+1:>4} {WT[i]:>3} {len(c):>5}  {100*c[WT[i]]/tot:>7.1f}%  {top}")

print("\n== binding class by edit distance from WT ==")
tab = defaultdict(Counter)
for s, c, ed in facs:
    tab[ed][c] += 1
print(f"{'ed':>3} {'n':>8} {'high%':>7} {'med%':>7} {'low%':>7}")
for ed in sorted(tab):
    c = tab[ed]; n = sum(c.values())
    print(f"{ed:>3} {n:>8} {100*c['high']/n:>7.1f} {100*c['medium']/n:>7.1f} {100*c['low']/n:>7.1f}")

# ---------------------------------------------------------------- load SPR
spr = list(csv.DictReader(open(f"{BASE}/trastuzumab_her2_spr_all_768.csv")))
print(f"\n\nSPR: {len(spr)} rows")
print("  binding_result:", dict(Counter(d["binding_result"] for d in spr)))
print("  cdrh3 lengths:", dict(Counter(len(d["cdrh3"]) for d in spr)))
vh_set = Counter(d["vh_sequence"][:20] + "..." for d in spr)
print("  n distinct VH:", len(set(d["vh_sequence"] for d in spr)),
      " n distinct VL:", len(set(d["vl_sequence"] for d in spr)))

# the trastuzumab reference rows
ref = [d for d in spr if d["cdrh3"] == WT]
print(f"\n  WT CDRH3 rows: {len(ref)}; KDs:",
      [d["kd_nM"] for d in ref])
print("  other-antibody rows:", sum(1 for d in spr if d["cdrh3"] == OTHER_AB))

# trastuzumab-variant subset: drop the other antibody
tz = [d for d in spr if d["cdrh3"] != OTHER_AB]
# is the VH outside CDRH3 identical across tz rows?
WT_VH = ref[0]["vh_sequence"]
print(f"\n  WT VH len {len(WT_VH)}; CDRH3 found at index {WT_VH.find(WT)}")
odd = []
for d in tz:
    vh, s = d["vh_sequence"], d["cdrh3"]
    j = vh.find(s)
    if j < 0:
        odd.append(("cdrh3 not in vh", s))
        continue
    scaffold = vh[:j] + vh[j + len(s):]
    wt_j = WT_VH.find(WT)
    wt_scaffold = WT_VH[:wt_j] + WT_VH[wt_j + len(WT):]
    if scaffold != wt_scaffold:
        odd.append(("framework differs", s))
print(f"  rows whose VH framework differs from WT outside CDRH3: {len(odd)}")
for o in odd[:10]:
    print("   ", o)

quant = [d for d in tz if d["binding_result"] == "quantified_binding" and d["kd_nM"]]
print(f"\n  quantified trastuzumab variants: {len(quant)}")
kds = np.array([float(d["kd_nM"]) for d in quant])
print(f"  KD nM: min {kds.min():.3f}  median {np.median(kds):.1f}  max {kds.max():.0f}")

def ed(a, b):
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))

print("\n== quantified variants by CDRH3 length and edit distance ==")
t2 = defaultdict(list)
for d in quant:
    s = d["cdrh3"]
    t2[(len(s), ed(s, WT))].append(float(d["kd_nM"]))
for k in sorted(t2, key=lambda k: (k[0], k[1] if k[1] is not None else -1)):
    v = np.array(t2[k])
    print(f"  len {k[0]}  ed {str(k[1]):>4}  n={len(v):>3}  "
          f"KD median {np.median(v):>9.2f}  best {v.min():>8.3f}")

print("\n== top 30 tightest binders (10-mers, excluding WT sequence) ==")
q10 = sorted([d for d in quant if len(d["cdrh3"]) == 10],
             key=lambda d: float(d["kd_nM"]))
print(f"{'cdrh3':<12}{'ed':>3}{'KD_nM':>11}  {'FACS':>7}  mutations")
for d in q10[:30]:
    s = d["cdrh3"]
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i]) or "-WT-"
    print(f"{s:<12}{ed(s,WT):>3}{float(d['kd_nM']):>11.3f}  {facs_cls.get(s,'-'):>7}  {muts}")

print("\n== worst 10 (weakest) 10-mers ==")
for d in q10[-10:]:
    s = d["cdrh3"]
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i]) or "-WT-"
    print(f"{s:<12}{ed(s,WT):>3}{float(d['kd_nM']):>11.1f}  {facs_cls.get(s,'-'):>7}  {muts}")

# -------------------------------------------- THE key question: FACS vs KD
print("\n\n================ FACS class vs measured KD (overlap) ================")
ov = [(d["cdrh3"], float(d["kd_nM"]), facs_cls[d["cdrh3"]])
      for d in quant if d["cdrh3"] in facs_cls]
print(f"quantified SPR variants also in FACS table: {len(ov)}")
byc = defaultdict(list)
for s, k, c in ov:
    byc[c].append(k)
for c in ("high", "medium", "low"):
    v = np.array(byc.get(c, []))
    if len(v):
        print(f"  FACS {c:<6} n={len(v):>3}  KD nM: median {np.median(v):>9.2f} "
              f" geo-mean {np.exp(np.mean(np.log(v))):>9.2f}  range {v.min():.2f}-{v.max():.0f}")
# all SPR rows (not just quantified) vs FACS class
print("\n  binding_result of ALL SPR trastuzumab rows present in FACS, by FACS class:")
t3 = defaultdict(Counter)
for d in tz:
    if d["cdrh3"] in facs_cls:
        t3[facs_cls[d["cdrh3"]]][d["binding_result"]] += 1
for c in ("high", "medium", "low"):
    print(f"   {c:<7}", dict(t3[c]))

json.dump({"n_facs": len(facs_cls), "n_quant": len(quant)},
          open("/home/hackathon/jakob/landscape.json", "w"))
