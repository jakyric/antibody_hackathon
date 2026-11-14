#!/usr/bin/env python
"""Part 2: the challenge filter, exactly; and per-position / per-mutation effects.

Filter (12-her2.md): a design passes if its CDRH3
  (a) differs at >=2 of its 10 positions from EVERY sequence in the provided data, and
  (b) is at most 70% identical to WGGDGFYAMD  (i.e. >=3 of 10 positions changed).

(a) is equivalent to: neither the candidate nor any of its 190 single-mutant
neighbours is in the provided set -> O(1) per candidate via a hash set.
"""
import csv, zipfile, io, math, itertools, random
from collections import Counter, defaultdict
import numpy as np

WT = "WGGDGFYAMD"
AAs = "ACDEFGHIKLMNPQRSTVWY"
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"
OTHER_AB = "ARDRIVDDYYYYGMDV"

facs = {}
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    n = [x for x in z.namelist() if x.endswith(".csv") and "MACOS" not in x][0]
    with z.open(n) as fh:
        for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
            facs[d["cdrh3"]] = d["binding_class"]

spr = list(csv.DictReader(open(f"{BASE}/trastuzumab_her2_spr_all_768.csv")))
spr_seqs = set(d["cdrh3"] for d in spr)

print("== what makes up the 525,063 'provided data' sequences? ==")
print(f"  unique FACS                 : {len(facs)}")
print(f"  unique SPR                  : {len(spr_seqs)}")
print(f"  union                       : {len(set(facs) | spr_seqs)}")
print(f"  union, 10-mers only         : {len([s for s in set(facs)|spr_seqs if len(s)==10])}")
provided = set(facs) | spr_seqs
provided10 = set(s for s in provided if len(s) == 10)
print(f"  -> doc says 525,063; union = {len(set(facs)|spr_seqs)}")

def passes_novelty(s, prov=provided10):
    """(a) min Hamming to provided 10-mers >= 2."""
    if s in prov:
        return False
    for i in range(10):
        for a in AAs:
            if a == s[i]:
                continue
            if s[:i] + a + s[i+1:] in prov:
                return False
    return True

def ed_wt(s):
    return sum(x != y for x, y in zip(s, WT))

def passes_filter(s):
    return ed_wt(s) >= 3 and passes_novelty(s)

# sanity: every provided sequence must fail (a)
assert not passes_novelty(WT)
print("\n  self-check: a random provided sequence passes novelty?",
      passes_novelty(next(iter(provided10))))

# ------------------------------------------------- how covered is the space?
print("\n== how much of 10-mer space does the provided data block? ==")
blocked = set()
for s in provided10:
    blocked.add(s)
print(f"  provided 10-mers: {len(provided10)}")
print(f"  total 10-mer space: 20^10 = {20**10:.3e} -- blocking is negligible globally,")
print("  but the library is concentrated, so blocking matters near the good region.")

# ------------------------------------- do the best measured binders pass?
tz = [d for d in spr if d["cdrh3"] != OTHER_AB and len(d["cdrh3"]) == 10]
quant = [d for d in tz if d["binding_result"] == "quantified_binding" and d["kd_nM"]]
quant.sort(key=lambda d: float(d["kd_nM"]))
print("\n== would the best measured variants themselves pass the filter? ==")
print("   (they can't pass (a) - they ARE provided data - so we report ed_wt only,")
print("    plus how many of their 1-neighbours are free, i.e. local novelty room)")
print(f"{'cdrh3':<12}{'KD_nM':>10}{'ed_wt':>6}{'ok_edwt':>9}  free_1nb/190  free_2nb_pass/sample")
for d in quant[:15]:
    s = d["cdrh3"]
    nb1 = [s[:i]+a+s[i+1:] for i in range(10) for a in AAs if a != s[i]]
    free1 = sum(1 for x in nb1 if x not in provided10)
    # 2-mutant neighbours of s that pass the full filter
    nb2 = set()
    for i, j in itertools.combinations(range(10), 2):
        for a in AAs:
            if a == s[i]: continue
            for b in AAs:
                if b == s[j]: continue
                t = list(s); t[i] = a; t[j] = b
                nb2.add("".join(t))
    nb2 = list(nb2)
    samp = random.Random(0).sample(nb2, 400)
    npass = sum(1 for x in samp if passes_filter(x))
    print(f"{s:<12}{float(d['kd_nM']):>10.3f}{ed_wt(s):>6}"
          f"{'yes' if ed_wt(s)>=3 else 'NO':>9}  {free1:>8}/190  {npass}/400 "
          f"(~{100*npass/400:.0f}% of {len(nb2)} 2-mutants)")

# ------------------------------------------------- per-position tolerance, FACS
print("\n\n== FACS: per-position, per-residue log2 enrichment (high vs low) ==")
hi = [Counter() for _ in range(10)]
lo = [Counter() for _ in range(10)]
for s, c in facs.items():
    if c == "high":
        for i, a in enumerate(s): hi[i][a] += 1
    elif c == "low":
        for i, a in enumerate(s): lo[i][a] += 1
nhi, nlo = sum(hi[0].values()), sum(lo[0].values())

enrich = {}
for i in range(10):
    for a in AAs:
        h, l = hi[i][a], lo[i][a]
        if h + l < 100:
            enrich[(i, a)] = None
            continue
        enrich[(i, a)] = math.log2(((h+1)/(nhi+20)) / ((l+1)/(nlo+20)))

print("  (blank = <100 observations at that position)")
hdr = "pos WT  " + "".join(f"{a:>6}" for a in AAs)
print(hdr)
for i in range(10):
    row = f"{i+1:>3} {WT[i]:>2}  "
    for a in AAs:
        v = enrich[(i, a)]
        row += "      " if v is None else f"{v:>+6.2f}"
    print(row)

print("\n  position tolerance = spread of enrichment (max-min) and WT's rank:")
print(f"{'pos':>4}{'WT':>4}{'WT_enrich':>11}{'best_aa':>9}{'best':>7}{'spread':>8}  residues better than WT")
for i in range(10):
    vals = {a: v for a, v in ((a, enrich[(i, a)]) for a in AAs) if v is not None}
    best = max(vals, key=vals.get)
    wtv = vals.get(WT[i])
    better = [a for a in sorted(vals, key=lambda a: -vals[a])
              if wtv is not None and vals[a] > wtv]
    print(f"{i+1:>4}{WT[i]:>4}{(f'{wtv:+.2f}' if wtv is not None else 'n/a'):>11}"
          f"{best:>9}{vals[best]:>+7.2f}{max(vals.values())-min(vals.values()):>8.2f}  "
          f"{''.join(better)}")

# ------------------------------------------------- per-mutation effects, SPR
print("\n\n== SPR: effect of each single mutation on log10 KD ==")
print("   (from all quantified 10-mers: mean log10KD of variants carrying the")
print("    mutation minus mean log10KD of those that do not, at matched ed)")
rows = [(d["cdrh3"], math.log10(float(d["kd_nM"]))) for d in quant]
glob = np.mean([k for _, k in rows])
eff = defaultdict(list)
for s, k in rows:
    for i in range(10):
        if s[i] != WT[i]:
            eff[(i, s[i])].append(k)
print(f"{'mut':>8}{'n':>5}{'mean_log10KD':>14}{'median_KD_nM':>14}")
items = [(k, v) for k, v in eff.items() if len(v) >= 4]
for (i, a), v in sorted(items, key=lambda kv: np.mean(kv[1])):
    print(f"{WT[i]}{i+1}{a:<5}{len(v):>5}{np.mean(v):>14.2f}{10**np.median(v):>14.2f}")

print("\n== SPR: per-position, does mutating it hurt? ==")
print(f"{'pos':>4}{'WT':>4}{'n_mut':>7}{'mean_log10KD_if_mut':>21}{'if_WT':>8}{'delta':>8}")
for i in range(10):
    m = [k for s, k in rows if s[i] != WT[i]]
    w = [k for s, k in rows if s[i] == WT[i]]
    if m and w:
        print(f"{i+1:>4}{WT[i]:>4}{len(m):>7}{np.mean(m):>21.2f}{np.mean(w):>8.2f}"
              f"{np.mean(m)-np.mean(w):>+8.2f}")
