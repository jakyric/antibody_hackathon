#!/usr/bin/env python
"""Part 10: retrospective test of the whole pipeline.

Hide all sub-nanomolar ed>=3 variants from training. Re-derive the position
constraints and the scorer from what is left. Then ask: does the pipeline rank
the hidden winners near the top of the measured ed>=3 pool?

This is the only honest test of the recipe, because it mimics the real task:
score sequences whose affinity you have never seen.
"""
import csv, zipfile, io
from collections import defaultdict
import numpy as np
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr

WT = "WGGDGFYAMD"
AAs = "ACDEFGHIKLMNPQRSTVWY"; AI = {a: i for i, a in enumerate(AAs)}
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"
OTHER = "ARDRIVDDYYYYGMDV"
spr = list(csv.DictReader(open(f"{BASE}/trastuzumab_her2_spr_all_768.csv")))

def onehot(ss):
    X = np.zeros((len(ss), 200))
    for r, s in enumerate(ss):
        for i, a in enumerate(s): X[r, i*20+AI[a]] = 1
    return X
def ed(a, b=WT): return sum(x != y for x, y in zip(a, b))

agg = defaultdict(list)
for d in spr:
    if (d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
            and d["binding_result"] == "quantified_binding" and d["kd_nM"]):
        agg[d["cdrh3"]].append(float(d["kd_nM"]))
kd = {s: float(np.mean(v)) for s, v in agg.items()}
pool = sorted([s for s in kd if ed(s) >= 3])          # the filter-relevant regime
hidden = [s for s in pool if kd[s] < 1.0]
train = [s for s in kd if s not in hidden]
print(f"measured ed>=3 pool: {len(pool)}   hidden winners (sub-nM): {len(hidden)}")
print(f"training set after hiding them: {len(train)}")

Xt = onehot(train); yt = np.array([np.log10(kd[s]) for s in train])

# --- does the constraint 'freeze positions 1,5,6,7' survive without them?
m = Ridge(alpha=10).fit(Xt, yt); C = m.coef_.reshape(10, 20)
nobs = np.zeros((10, 20))
for s in train:
    for i, a in enumerate(s): nobs[i, AI[a]] += 1
print("\n== is WT still the best residue at each position, winners hidden? ==")
print(f"{'pos':>4}{'WT':>4}{'WT_rank_among_observed':>24}{'best':>6}")
for i in range(10):
    obs = sorted([(a, C[i, AI[a]]) for a in AAs if nobs[i, AI[a]] >= 4],
                 key=lambda t: t[1])
    rk = [a for a, _ in obs].index(WT[i]) + 1 if any(a == WT[i] for a, _ in obs) else None
    print(f"{i+1:>4}{WT[i]:>4}{f'{rk} of {len(obs)}':>24}{obs[0][0]:>6}")

# --- rank the whole measured ed>=3 pool with this scorer
Xp = onehot(pool); pp = m.predict(Xp)
yp = np.array([np.log10(kd[s]) for s in pool])
o = np.argsort(pp)
print(f"\n== ranking of the {len(pool)} measured ed>=3 variants by the scorer ==")
print(f"   Spearman over the whole pool: {spearmanr(pp, yp)[0]:+.3f}")
ranks = {pool[j]: r for r, j in enumerate(o, 1)}
print(f"\n   where the {len(hidden)} hidden winners landed (out of {len(pool)}):")
for s in sorted(hidden, key=lambda s: kd[s]):
    print(f"     {s}  measured {kd[s]:>7.3f} nM   rank {ranks[s]:>3}/{len(pool)}"
          f"   (top {100*ranks[s]/len(pool):>4.1f}%)")
med = np.median([ranks[s] for s in hidden])
print(f"   median rank of hidden winners: {med:.0f}/{len(pool)} "
      f"(random would be {len(pool)/2:.0f})")

for N in (10, 20, 30, 50):
    got = sum(1 for s in hidden if ranks[s] <= N)
    base = N * len(hidden) / len(pool)
    print(f"   top {N:>3}: captured {got}/{len(hidden)} winners "
          f"(expected {base:.1f} by chance, enrichment {got/base if base else 0:.1f}x)")

# --- how much of this is just 'edit distance' within ed>=3?
edp = np.array([ed(s) for s in pool], float)
print(f"\n   for reference, Spearman(ed, log10KD) within the ed>=3 pool: "
      f"{spearmanr(edp, yp)[0]:+.3f}")
oe = np.argsort(edp)
re_ = {pool[j]: r for r, j in enumerate(oe, 1)}
print(f"   median rank of winners if you rank by edit distance alone: "
      f"{np.median([re_[s] for s in hidden]):.0f}/{len(pool)}")

# --- and the FACS-only alternative, for completeness
facs = {}
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    n = [x for x in z.namelist() if x.endswith(".csv") and "MACOS" not in x][0]
    with z.open(n) as fh:
        for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
            facs[d["cdrh3"]] = d["binding_class"]
from sklearn.linear_model import LogisticRegression
hl = [(s, 1 if c == "high" else 0) for s, c in facs.items() if c in ("high", "low")]
Xf = onehot([s for s, _ in hl]); yf = np.array([v for _, v in hl])
cf = LogisticRegression(max_iter=3000).fit(Xf, yf)
CF = cf.coef_.reshape(10, 20)
fs = np.array([sum(CF[i, AI[s[i]]] for i in range(10)) for s in pool])
of = np.argsort(-fs)
rf = {pool[j]: r for r, j in enumerate(of, 1)}
print(f"\n   median rank of winners if you rank by the 524k-sequence FACS "
      f"classifier: {np.median([rf[s] for s in hidden]):.0f}/{len(pool)}")
print(f"   Spearman(FACS score, log10KD) in this pool: {spearmanr(-fs, yp)[0]:+.3f}")
print("\n   summary of rankers on the ed>=3 pool (lower median rank = better):")
print(f"     SPR additive model (winners hidden) : {med:.0f}")
print(f"     edit distance from WT               : {np.median([re_[s] for s in hidden]):.0f}")
print(f"     FACS classifier (524k sequences)    : {np.median([rf[s] for s in hidden]):.0f}")
print(f"     random                              : {len(pool)/2:.0f}")
