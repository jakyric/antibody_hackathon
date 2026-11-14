#!/usr/bin/env python
"""Part 5: where the model fails, and the champion-anchored alternative.

A. Calibration: the ridge shrinks hard - how badly, and does that matter?
B. Can it rank WITHIN the tight binders? (the actual objective)
C. Non-linear models and epistasis evidence.
D. Leave-combination-out: can it predict a mutation combination it never saw?
E. Champion-anchored candidate pool: >=2 substitutions from WGSQGFYAYD (23 pM).
"""
import csv, zipfile, io, math, itertools
from collections import Counter, defaultdict
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import KFold
from scipy.stats import spearmanr

WT = "WGGDGFYAMD"
CHAMP = "WGSQGFYAYD"
AAs = "ACDEFGHIKLMNPQRSTVWY"
AI = {a: i for i, a in enumerate(AAs)}
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
        for i, a in enumerate(s):
            X[r, i*20 + AI[a]] = 1
    return X
def ed(a, b): return sum(x != y for x, y in zip(a, b))
def passes(s):
    if ed(s, WT) < 3 or s in provided10: return False
    for i in range(10):
        for a in AAs:
            if a != s[i] and s[:i]+a+s[i+1:] in provided10: return False
    return True

agg = defaultdict(list)
for d in spr:
    if (d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
            and d["binding_result"] == "quantified_binding" and d["kd_nM"]):
        agg[d["cdrh3"]].append(float(d["kd_nM"]))
seqs = sorted(agg); y = np.array([np.mean(np.log10(agg[s])) for s in seqs])
X = onehot(seqs)

# ---------------------------------------------------- A. calibration
print("== A. calibration of the additive ridge (5-fold CV predictions) ==")
kf = KFold(5, shuffle=True, random_state=1)
pcv = np.zeros_like(y)
for a, b in kf.split(X):
    pcv[b] = Ridge(alpha=10).fit(X[a], y[a]).predict(X[b])
sl, ic = np.polyfit(y, pcv, 1)
print(f"   regression of predicted on measured: slope {sl:.2f} (1.0 = calibrated),"
      f" intercept {ic:.2f}")
print(f"   -> predictions are compressed {1/sl:.1f}x toward the mean "
      f"({10**np.mean(y):.0f} nM)")
print("\n   measured decile vs mean predicted:")
q = np.quantile(y, np.linspace(0, 1, 6))
for i in range(5):
    m = (y >= q[i]) & (y <= q[i+1])
    print(f"     measured {10**q[i]:>9.3f}-{10**q[i+1]:>9.1f} nM (n={m.sum():>3}): "
          f"mean predicted {10**pcv[m].mean():>9.2f} nM")
print("\n   consequence: absolute predicted KD is meaningless for tight binders;")
print("   only the RANK is usable, and even the rank must be checked at the top end.")

# ---------------------------------------------------- B. ranking at the top
print("\n== B. can it rank WITHIN the tight binders? ==")
for thr in (1.0, 10.0, 100.0):
    m = y < np.log10(thr)
    if m.sum() > 8:
        r = spearmanr(pcv[m], y[m])[0]
        print(f"   restricted to measured KD < {thr:>5.0f} nM (n={m.sum():>3}): "
              f"Spearman {r:+.3f}")
print("   (full set Spearman was +0.71)")
m = y < 0
print(f"\n   the {m.sum()} sub-nanomolar variants, measured vs CV-predicted:")
for j in np.argsort(y)[:m.sum()]:
    print(f"     {seqs[j]}  ed {ed(seqs[j],WT)}  measured {10**y[j]:>7.3f} nM   "
          f"predicted {10**pcv[j]:>8.2f} nM")

# ---------------------------------------------------- C. non-linear / epistasis
print("\n== C. non-linear models (same CV) ==")
def cvscore(est, X, n_rep=5):
    sp = []
    for rep in range(n_rep):
        kf = KFold(5, shuffle=True, random_state=rep)
        pr = np.zeros_like(y)
        for a, b in kf.split(X):
            pr[b] = est().fit(X[a], y[a]).predict(X[b])
        sp.append(spearmanr(pr, y)[0])
    return np.mean(sp), np.std(sp)
for nm, est in [("ridge additive", lambda: Ridge(alpha=10)),
                ("random forest", lambda: RandomForestRegressor(
                    400, min_samples_leaf=2, random_state=0, n_jobs=-1)),
                ("gradient boosting", lambda: GradientBoostingRegressor(
                    random_state=0, max_depth=2, n_estimators=300,
                    learning_rate=0.05))]:
    s, sd = cvscore(est, X)
    print(f"   {nm:<20} Spearman {s:+.3f} +- {sd:.3f}")
# pairwise terms
Xp = [X]
cols = []
for i, j in itertools.combinations(range(10), 2):
    pass
print("   (pairwise one-hot = 19,000 features for 180 points: not fittable;")
print("    instead test whether residuals are structured)")
m = Ridge(alpha=10).fit(X, y)
res = y - m.predict(X)
print(f"   residual sd {res.std():.2f} log10 units; "
      f"Spearman(|residual|, ed) = {spearmanr(np.abs(res), [ed(s,WT) for s in seqs])[0]:+.3f}")
worst = np.argsort(-np.abs(res))[:8]
print("   8 largest residuals (model badly wrong -> candidate epistasis):")
for j in worst:
    print(f"     {seqs[j]}  ed {ed(seqs[j],WT):>2}  measured {10**y[j]:>9.3f} nM  "
          f"fitted {10**m.predict(X[j:j+1])[0]:>9.2f} nM  resid {res[j]:+.2f}")

# ---------------------------------------------------- D. leave-combination-out
print("\n== D. the extrapolation we actually rely on: predict a NEW combination ==")
print("   leave-one-out, but also removing every training variant that shares")
print("   >=2 mutated positions with the held-out one (so the combination is unseen)")
mutpos = [frozenset(i for i in range(10) if s[i] != WT[i]) for s in seqs]
preds, trues, ntr = [], [], []
for j in range(len(seqs)):
    keep = [k for k in range(len(seqs))
            if k != j and len(mutpos[k] & mutpos[j]) < 2]
    if len(keep) < 30: continue
    mm = Ridge(alpha=10).fit(X[keep], y[keep])
    preds.append(mm.predict(X[j:j+1])[0]); trues.append(y[j]); ntr.append(len(keep))
preds, trues = np.array(preds), np.array(trues)
print(f"   n evaluated {len(preds)} (mean train size {np.mean(ntr):.0f})")
print(f"   Spearman {spearmanr(preds, trues)[0]:+.3f}   "
      f"MAE {np.mean(np.abs(preds-trues)):.3f} log10")
print("   -> this is the honest number for ranking unseen mutation combinations.")

# ---------------------------------------------------- E. champion-anchored pool
print("\n== E. champion-anchored pool: >=2 substitutions from WGSQGFYAYD ==")
print(f"   anchor {CHAMP}: measured {np.mean(agg[CHAMP]):.3f} nM, "
      f"ed_wt {ed(CHAMP,WT)}, best in the whole dataset")
model = Ridge(alpha=10).fit(X, y)
# conservative alphabet: residues observed >=4x and not catastrophic
nobs = np.zeros((10, 20))
for s in seqs:
    for i, a in enumerate(s): nobs[i, AI[a]] += 1
C = model.coef_.reshape(10, 20)
alpha = {}
for i in range(10):
    wtc = C[i, AI[WT[i]]]
    keep = [a for a in AAs if nobs[i, AI[a]] >= 4 and C[i, AI[a]] <= wtc + 0.45]
    for extra in (WT[i], CHAMP[i]):
        if extra not in keep: keep.append(extra)
    alpha[i] = keep
pool = []
for s in ("".join(t) for t in itertools.product(*[alpha[i] for i in range(10)])):
    if ed(s, CHAMP) >= 2 and ed(s, CHAMP) <= 3 and passes(s):
        pool.append(s)
print(f"   candidates 2-3 substitutions from the champion that pass the filter: {len(pool)}")
if pool:
    pp = model.predict(onehot(pool))
    o = np.argsort(pp)
    print(f"{'cdrh3':<12}{'ed_wt':>6}{'ed_champ':>9}{'pred_nM':>10}  diff-from-champion")
    for j in o[:20]:
        s = pool[j]
        dm = ",".join(f"{CHAMP[i]}{i+1}{s[i]}" for i in range(10) if s[i] != CHAMP[i])
        print(f"{s:<12}{ed(s,WT):>6}{ed(s,CHAMP):>9}{10**pp[j]:>10.3f}  {dm}")
    with open("/home/hackathon/jakob/candidates_champion_anchored.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["cdrh3","ed_wt","ed_champion","pred_kd_nM","diff_from_champion"])
        for j in o:
            s = pool[j]
            w.writerow([s, ed(s,WT), ed(s,CHAMP), f"{10**pp[j]:.4f}",
                        ";".join(f"{CHAMP[i]}{i+1}{s[i]}" for i in range(10) if s[i]!=CHAMP[i])])
    print(f"\n   wrote ~/jakob/candidates_champion_anchored.csv ({len(pool)} rows)")
