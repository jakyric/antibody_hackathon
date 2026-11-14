#!/usr/bin/env python
"""Part 11: the tiered panel implied by the retrospective test.

With the 10 sub-nM winners hidden, WT was still the best residue at positions
1,2,5,6,7,8,9 and was beaten only at positions 3, 4 and 10. The challenge
filter needs >=3 changes. So the minimum-risk design mutates exactly 3, 4, 10.
"""
import csv, itertools, zipfile, io
from collections import defaultdict
import numpy as np
from sklearn.linear_model import Ridge

WT = "WGGDGFYAMD"
AAs = "ACDEFGHIKLMNPQRSTVWY"; AI = {a: i for i, a in enumerate(AAs)}
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"
OTHER = "ARDRIVDDYYYYGMDV"

facs = {}
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    n = [x for x in z.namelist() if x.endswith(".csv") and "MACOS" not in x][0]
    with z.open(n) as fh:
        for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
            facs[d["cdrh3"]] = d["binding_class"]
spr = list(csv.DictReader(open(f"{BASE}/trastuzumab_her2_spr_all_768.csv")))
provided10 = set(s for s in (set(facs) | {d["cdrh3"] for d in spr}) if len(s) == 10)
VH = [d["vh_sequence"] for d in spr if d["cdrh3"] == WT][0]
VL = [d["vl_sequence"] for d in spr if d["cdrh3"] == WT][0]

def onehot(ss):
    X = np.zeros((len(ss), 200))
    for r, s in enumerate(ss):
        for i, a in enumerate(s): X[r, i*20+AI[a]] = 1
    return X
def ed(a, b=WT): return sum(x != y for x, y in zip(a, b))
def passes(s):
    if ed(s) < 3 or s in provided10: return False
    for i in range(10):
        for a in AAs:
            if a != s[i] and s[:i]+a+s[i+1:] in provided10: return False
    return True

agg = defaultdict(list)
for d in spr:
    if (d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
            and d["binding_result"] == "quantified_binding" and d["kd_nM"]):
        agg[d["cdrh3"]].append(float(d["kd_nM"]))
kd = {s: float(np.mean(v)) for s, v in agg.items()}
seqs = sorted(kd)
X = onehot(seqs); y = np.array([np.log10(kd[s]) for s in seqs])
m_all = Ridge(alpha=10).fit(X, y)
m3 = np.array([ed(s) >= 3 for s in seqs])
m_ed3 = Ridge(alpha=10).fit(X[m3], y[m3])
# leave-winners-out model (the unbiased one from part 10)
hid = [s for s in seqs if ed(s) >= 3 and kd[s] < 1.0]
tr = [s for s in seqs if s not in hid]
m_lwo = Ridge(alpha=10).fit(onehot(tr), np.array([np.log10(kd[s]) for s in tr]))

def score3(s):
    X1 = onehot([s])
    return (10**m_all.predict(X1)[0], 10**m_ed3.predict(X1)[0],
            10**m_lwo.predict(X1)[0])

LIAB = ("NG", "NS", "DG", "DS")
def liab(s):
    out = [x for x in LIAB if x in s and x not in WT]
    if "C" in s: out.append("freeCys")
    for i in range(8):
        if s[i] == "N" and s[i+1] != "P" and s[i+2] in "ST": out.append("Nglyc")
    return out

print("="*78)
print("TIER A - minimum risk: change ONLY positions 3, 4, 10 (ed = 3 exactly)")
print("  the three positions where WT lost even with all winners hidden,")
print("  and all three make zero contact with HER2 in 1N8Z")
print("="*78)
A3 = "SVDAYTNG"      # pos3, observed n>=4, best coefficients
A4 = "QESYTCNV"      # pos4
A10 = "EITGVAN"      # pos10
tierA = []
for a, b, c in itertools.product(A3, A4, A10):
    s = WT[0:2] + a + b + WT[4:9] + c
    if passes(s) and not liab(s):
        tierA.append(s)
rows = [(s,) + score3(s) for s in tierA]
rows.sort(key=lambda r: (r[1]*r[2]*r[3])**(1/3))
print(f"  {len(tierA)} designs pass the filter and carry no liability motif\n")
print(f"{'#':>3} {'cdrh3':<12}{'pred_all':>10}{'pred_ed3':>10}{'pred_lwo':>10}   mutations")
for i, (s, a, b, c) in enumerate(rows[:15], 1):
    muts = ",".join(f"{WT[j]}{j+1}{s[j]}" for j in range(10) if s[j] != WT[j])
    print(f"{i:>3} {s:<12}{a:>10.2f}{b:>10.2f}{c:>10.2f}   {muts}")

print("\n" + "="*78)
print("TIER B - add position 9 and/or 8 (8 of the 10 measured sub-nM")
print("  ed>=3 variants mutate position 9; M9Y/F/L/W/H all appear there)")
print("="*78)
B = []
for a, b, e8, e9, c in itertools.product("SVDAY", "QESY", "AELSMV", "MYFLWH", "DEITGA"):
    s = WT[0:2] + a + b + WT[4:7] + e8 + e9 + c
    if ed(s) >= 3 and passes(s) and not liab(s):
        B.append(s)
rowsB = [(s,) + score3(s) for s in B]
rowsB.sort(key=lambda r: (r[1]*r[2]*r[3])**(1/3))
picked = []
for r in rowsB:
    if all(ed(r[0], t[0]) >= 2 for t in picked):
        picked.append(r)
    if len(picked) == 15: break
print(f"  {len(B)} designs pass; 15 shown, >=2 apart from each other\n")
print(f"{'#':>3} {'cdrh3':<12}{'ed':>3}{'pred_all':>10}{'pred_ed3':>10}{'pred_lwo':>10}   mutations")
for i, (s, a, b, c) in enumerate(picked, 1):
    muts = ",".join(f"{WT[j]}{j+1}{s[j]}" for j in range(10) if s[j] != WT[j])
    print(f"{i:>3} {s:<12}{ed(s):>3}{a:>10.2f}{b:>10.2f}{c:>10.2f}   {muts}")

print("\n== reference: the same three scorers on the known reference points ==")
print(f"{'cdrh3':<12}{'measured':>10}{'pred_all':>10}{'pred_ed3':>10}{'pred_lwo':>10}")
for s in ("WGGDGFYAMD", "WGSQGFYAYD", "WGGSGPYLFD", "WTYNSFYAMD", "WKGYGSYVHD"):
    a, b, c = score3(s)
    print(f"{s:<12}{kd.get(s, float('nan')):>10.3f}{a:>10.2f}{b:>10.2f}{c:>10.2f}")
print("\n  all three scorers under-predict tight binders by ~1-1.5 log units")
print("  (ridge shrinkage). Use them to RANK, never as a KD estimate.")

# write submission-ready file
out = "/home/hackathon/jakob/panel_tiered.csv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["tier", "rank", "cdrh3", "ed_wt", "pred_kd_nM_allfit",
                "pred_kd_nM_ed3fit", "pred_kd_nM_winnersheldout",
                "mutations", "vh_sequence", "vl_sequence"])
    for i, (s, a, b, c) in enumerate(rows[:15], 1):
        w.writerow(["A", i, s, ed(s), f"{a:.3f}", f"{b:.3f}", f"{c:.3f}",
                    ";".join(f"{WT[j]}{j+1}{s[j]}" for j in range(10) if s[j] != WT[j]),
                    VH[:98] + s + VH[108:], VL])
    for i, (s, a, b, c) in enumerate(picked, 1):
        w.writerow(["B", i, s, ed(s), f"{a:.3f}", f"{b:.3f}", f"{c:.3f}",
                    ";".join(f"{WT[j]}{j+1}{s[j]}" for j in range(10) if s[j] != WT[j]),
                    VH[:98] + s + VH[108:], VL])
print(f"\nwrote {out} (30 designs, tiers A and B, with full VH/VL)")
