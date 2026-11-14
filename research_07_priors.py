#!/usr/bin/env python
"""Part 7: base rates, and is the FACS data useful as a BINDING gate
(rather than as an affinity ranker, which part 3 showed it is not)?"""
import csv, zipfile, io, math
from collections import defaultdict, Counter
import numpy as np
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score

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

def onehot(ss):
    X = np.zeros((len(ss), 200))
    for r, s in enumerate(ss):
        for i, a in enumerate(s): X[r, i*20+AI[a]] = 1
    return X
def ed(a, b=WT): return sum(x != y for x, y in zip(a, b))

# ---------------------------------------------- base rates at ed>=3
q = defaultdict(list)
for d in spr:
    if (d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
            and d["binding_result"] == "quantified_binding" and d["kd_nM"]):
        q[d["cdrh3"]].append(float(d["kd_nM"]))
kd = {s: np.mean(v) for s, v in q.items()}
e3 = {s: k for s, k in kd.items() if ed(s) >= 3}
print("== base rates among MEASURED variants with >=3 changes from WT ==")
print(f"   (this is the regime the challenge filter forces you into)")
print(f"   n measured with ed>=3: {len(e3)}")
for thr in (0.2, 1, 10, 100):
    n = sum(1 for k in e3.values() if k < thr)
    print(f"     KD < {thr:>5} nM : {n:>3} / {len(e3)} = {100*n/len(e3):>5.1f}%")
print(f"   WT itself: {np.mean(kd[WT]):.3f} nM (8 replicates, "
      f"{min(q[WT]):.3f}-{max(q[WT]):.3f})")
print(f"\n   sub-nanomolar variants with ed>=3 (all of them):")
for s, k in sorted(e3.items(), key=lambda t: t[1]):
    if k < 1:
        muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i])
        print(f"     {s}  ed {ed(s)}  {k:>7.3f} nM   {muts}")

# ---------------------------------------------- FACS as a binding gate
print("\n== is FACS useful as a BINDING gate (not an affinity ranker)? ==")
hl = [(s, 1 if c == "high" else 0) for s, c in facs.items() if c in ("high", "low")]
Xf = onehot([s for s, _ in hl]); yf = np.array([v for _, v in hl])
clf = LogisticRegression(max_iter=3000).fit(Xf, yf)
CO = clf.coef_.reshape(10, 20)
def fscore(s): return sum(CO[i, AI[s[i]]] for i in range(10))

# task: among SPR trastuzumab 10-mers, does FACS score predict whether SPR
# returned a KD at all?
rows = [(d["cdrh3"], d["binding_result"]) for d in spr
        if d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
        and d["binding_result"] in ("quantified_binding", "no_reportable_kd")]
seen = {}
for s, r in rows: seen.setdefault(s, r)
ss = list(seen); lab = np.array([1 if seen[s] == "quantified_binding" else 0 for s in ss])
sc = np.array([fscore(s) for s in ss])
edv = np.array([ed(s) for s in ss], float)
print(f"   n={len(ss)}  ({lab.sum()} gave a KD, {len(ss)-lab.sum()} did not)")
print(f"   AUC, FACS classifier score -> 'SPR returned a KD' : {roc_auc_score(lab, sc):.3f}")
print(f"   AUC, edit distance from WT (inverted)             : {roc_auc_score(lab, -edv):.3f}")
# combined
from sklearn.linear_model import LogisticRegression as LR
Xc = np.c_[sc, edv]
pr = np.zeros(len(ss))
from sklearn.model_selection import KFold
for a, b in KFold(5, shuffle=True, random_state=0).split(Xc):
    pr[b] = LR(max_iter=1000).fit(Xc[a], lab[a]).predict_proba(Xc[b])[:, 1]
print(f"   AUC, FACS score + ed (5-fold CV)                  : {roc_auc_score(lab, pr):.3f}")
print("   -> FACS DOES carry binder/non-binder signal that transfers to SPR;")
print("      part 3 showed it carries almost no AFFINITY signal.")
print("      Correct use: a gate on 'will it bind at all', not a ranker.")

# ---------------------------------------------- pool sizes
print("\n== pool sizes for the obvious alternative strategies ==")
hi = [s for s, c in facs.items() if c == "high"]
print(f"   FACS high binders: {len(hi)}  -- but every one of them is provided")
print(f"     data, so NONE can pass the novelty rule. ed>=3 among them: "
      f"{sum(1 for s in hi if ed(s)>=3)}")
# 2-mutant neighbours of FACS high sequences that pass the filter: sample
import random
rnd = random.Random(0)
def passes(s):
    if ed(s) < 3 or s in provided10: return False
    for i in range(10):
        for a in AAs:
            if a != s[i] and s[:i]+a+s[i+1:] in provided10: return False
    return True
samp = rnd.sample(hi, 200)
tot = ok = 0
for s in samp:
    for _ in range(20):
        t = list(s)
        for i in rnd.sample(range(10), 2):
            t[i] = rnd.choice(AAs.replace(s[i], ""))
        tot += 1
        ok += passes("".join(t))
print(f"   random 2-substitution neighbours of FACS-high sequences that pass "
      f"the filter: {100*ok/tot:.0f}% ({ok}/{tot})")
print("   -> novelty is easy to satisfy; the hard part is affinity.")
