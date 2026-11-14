#!/usr/bin/env python
"""Part 3: can we rank novel CDRH3 designs for AFFINITY?

Tests, in order of increasing honesty:
  M1  edit distance from WT alone               (the trivial baseline)
  M2  FACS high-vs-low position enrichment sum  (free, 524k sequences)
  M3  FACS-trained classifier P(high)           (ridge logistic on 524k)
  M4  ridge on one-hot, fitted to 182 SPR log10 KD
  M5  M4 + ed + FACS features
evaluated by random CV and by a harder extrapolation split.
"""
import csv, zipfile, io, math, random, itertools
from collections import Counter, defaultdict
import numpy as np
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.model_selection import KFold
from scipy.stats import spearmanr, pearsonr

WT = "WGGDGFYAMD"
AAs = "ACDEFGHIKLMNPQRSTVWY"
AI = {a: i for i, a in enumerate(AAs)}
BASE = "/home/hackathon/data/hackathon-data/her2-trastuzumab"
OTHER = "ARDRIVDDYYYYGMDV"
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- data
facs = {}
with zipfile.ZipFile(f"{BASE}/tz_her2_facs_3classes_524346.csv.zip") as z:
    n = [x for x in z.namelist() if x.endswith(".csv") and "MACOS" not in x][0]
    with z.open(n) as fh:
        for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")):
            facs[d["cdrh3"]] = d["binding_class"]

spr = list(csv.DictReader(open(f"{BASE}/trastuzumab_her2_spr_all_768.csv")))
quant = [(d["cdrh3"], float(d["kd_nM"])) for d in spr
         if d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
         and d["binding_result"] == "quantified_binding" and d["kd_nM"]]
# average replicates of the same sequence (WT appears 8x)
agg = defaultdict(list)
for s, k in quant:
    agg[s].append(k)
seqs = sorted(agg)
y = np.array([np.mean(np.log10(agg[s])) for s in seqs])
print(f"SPR regression set: {len(seqs)} unique 10-mers "
      f"({len(quant)} measurements), log10KD {y.min():.2f}..{y.max():.2f}")

def onehot(ss):
    X = np.zeros((len(ss), 200))
    for r, s in enumerate(ss):
        for i, a in enumerate(s):
            X[r, i*20 + AI[a]] = 1
    return X

def edwt(s):
    return sum(x != y_ for x, y_ in zip(s, WT))

# ------------------------------------------------- FACS-derived features
hi = [Counter() for _ in range(10)]; lo = [Counter() for _ in range(10)]
med = [Counter() for _ in range(10)]
for s, c in facs.items():
    tgt = hi if c == "high" else (lo if c == "low" else med)
    for i, a in enumerate(s):
        tgt[i][a] += 1
nhi, nlo = sum(hi[0].values()), sum(lo[0].values())
ENR = np.zeros((10, 20))
for i in range(10):
    for a in AAs:
        ENR[i, AI[a]] = math.log2(((hi[i][a]+1)/(nhi+20)) / ((lo[i][a]+1)/(nlo+20)))

def enrich_sum(s):
    return sum(ENR[i, AI[a]] for i, a in enumerate(s))

# ------------------------------------------- M3: classifier on all 524k FACS
Xf_s = [s for s, c in facs.items() if c in ("high", "low")]
yf = np.array([1 if facs[s] == "high" else 0 for s in Xf_s])
print(f"\nTraining FACS high/low classifier on {len(Xf_s)} sequences ...")
Xf = onehot(Xf_s)
clf = LogisticRegression(max_iter=2000, C=1.0)
# held-out check of the classifier itself
idx = rng.permutation(len(Xf)); cut = int(.8*len(Xf))
clf.fit(Xf[idx[:cut]], yf[idx[:cut]])
from sklearn.metrics import roc_auc_score
auc_in = roc_auc_score(yf[idx[cut:]], clf.predict_proba(Xf[idx[cut:]])[:, 1])
print(f"  FACS classifier held-out AUC (high vs low, random split): {auc_in:.3f}")
clf.fit(Xf, yf)
FACS_COEF = clf.coef_.reshape(10, 20)

def facs_score(s):
    return sum(FACS_COEF[i, AI[s[i]]] for i in range(10))

# ------------------------------------------------- build feature matrices
X_oh = onehot(seqs)
f_ed = np.array([edwt(s) for s in seqs], float)
f_enr = np.array([enrich_sum(s) for s in seqs])
f_cls = np.array([facs_score(s) for s in seqs])

print("\n== univariate agreement with measured log10 KD (lower KD = better) ==")
for nm, v in [("edit distance from WT", f_ed),
              ("FACS log2 enrichment sum", -f_enr),
              ("FACS classifier score", -f_cls)]:
    r, p = spearmanr(v, y)
    print(f"  {nm:<28} Spearman r = {r:+.3f}  (p={p:.1e})")
print("  (sign flipped so that positive r = the feature ranks affinity correctly)")

# how much of the FACS features' signal is just edit distance?
print(f"\n  Spearman(FACS enrichment sum, ed)   = {spearmanr(f_enr, f_ed)[0]:+.3f}")
print(f"  Spearman(FACS classifier score, ed) = {spearmanr(f_cls, f_ed)[0]:+.3f}")
# partial: residualise the FACS score on ed, then correlate with y
for nm, v in [("FACS enrichment", f_enr), ("FACS classifier", f_cls)]:
    b = np.polyfit(f_ed, v, 1); res = v - np.polyval(b, f_ed)
    yb = np.polyfit(f_ed, y, 1); yres = y - np.polyval(yb, f_ed)
    print(f"  {nm} vs KD, both residualised on ed: r = {spearmanr(-res, yres)[0]:+.3f}")

# ------------------------------------------------- CV of ridge models
def cv(X, y, alpha=10.0, n_rep=20, k=5):
    sp, pr, mae = [], [], []
    for rep in range(n_rep):
        kf = KFold(k, shuffle=True, random_state=rep)
        pred = np.zeros_like(y)
        for tr, te in kf.split(X):
            m = Ridge(alpha=alpha).fit(X[tr], y[tr])
            pred[te] = m.predict(X[te])
        sp.append(spearmanr(pred, y)[0]); pr.append(pearsonr(pred, y)[0])
        mae.append(np.mean(np.abs(pred - y)))
    return np.mean(sp), np.std(sp), np.mean(pr), np.mean(mae)

print("\n== 5-fold CV (20 repeats) of models fitted to the 182 SPR values ==")
print(f"{'model':<42}{'Spearman':>10}{'+-':>7}{'Pearson':>9}{'MAE(log10)':>12}")
models = {
    "M1 ed only": f_ed.reshape(-1, 1),
    "M2 FACS enrichment sum only": f_enr.reshape(-1, 1),
    "M3 FACS classifier score only": f_cls.reshape(-1, 1),
    "M1+M3 ed + FACS classifier": np.c_[f_ed, f_cls],
    "M4 one-hot (additive, 200 feats)": X_oh,
    "M5 one-hot + ed + FACS feats": np.c_[X_oh, f_ed, f_enr, f_cls],
}
for nm, X in models.items():
    s, sd, p, m = cv(X, y)
    print(f"{nm:<42}{s:>10.3f}{sd:>7.3f}{p:>9.3f}{m:>12.3f}")

for a in (1, 3, 10, 30, 100):
    s, sd, p, m = cv(X_oh, y, alpha=a, n_rep=10)
    print(f"   one-hot alpha={a:<5} Spearman {s:.3f}  MAE {m:.3f}")

# ------------------------------------------------- harder: extrapolation split
print("\n== harder split: train on ed<=4, test on ed>=5 (extrapolating outward) ==")
tr = f_ed <= 4; te = f_ed >= 5
print(f"   n_train={tr.sum()} n_test={te.sum()}")
for nm, X in models.items():
    m = Ridge(alpha=10).fit(X[tr], y[tr])
    pr = m.predict(X[te])
    print(f"   {nm:<42} Spearman {spearmanr(pr, y[te])[0]:+.3f}  "
          f"MAE {np.mean(np.abs(pr-y[te])):.3f}")

print("\n== the question that matters: can a model find the tight binders? ==")
print("   enrichment of true sub-10 nM variants in the model's top-20 by CV prediction")
kf = KFold(5, shuffle=True, random_state=1)
for nm, X in models.items():
    pred = np.zeros_like(y)
    for a, b in kf.split(X):
        pred[b] = Ridge(alpha=10).fit(X[a], y[a]).predict(X[b])
    top = np.argsort(pred)[:20]
    hit = (y[top] < 1.0).sum()          # <10 nM
    base = (y < 1.0).sum()
    print(f"   {nm:<42} {hit}/20 sub-10nM in top-20 "
          f"(base rate {base}/{len(y)} = {100*base/len(y):.0f}%)")

# ------------------------------------------------- additive coefficients of M4
print("\n== M4 additive one-hot coefficients (fitted to all 182; "
      "negative = tighter) ==")
m = Ridge(alpha=10).fit(X_oh, y)
C = m.coef_.reshape(10, 20)
print("pos WT  " + "".join(f"{a:>6}" for a in AAs))
for i in range(10):
    row = f"{i+1:>3} {WT[i]:>2}  "
    for a in AAs:
        seen = sum(1 for s in seqs if s[i] == a)
        row += f"{C[i,AI[a]]:>+6.2f}" if seen >= 3 else "     ."
    print(row)
print("  '.' = fewer than 3 observations of that residue at that position")
print("\n  best (most affinity-improving) residue per position, among observed:")
for i in range(10):
    obs = [(a, C[i, AI[a]]) for a in AAs if sum(1 for s in seqs if s[i] == a) >= 3]
    obs.sort(key=lambda t: t[1])
    wtc = C[i, AI[WT[i]]]
    print(f"   pos {i+1} (WT {WT[i]} {wtc:+.2f}): " +
          "  ".join(f"{a}{'*' if a==WT[i] else ''}{c:+.2f}" for a, c in obs[:5]))

np.save("/home/hackathon/jakob/ridge_coef.npy", C)
np.save("/home/hackathon/jakob/facs_coef.npy", FACS_COEF)
np.save("/home/hackathon/jakob/facs_enrich.npy", ENR)
print("\nsaved coefficient matrices to ~/jakob/*.npy")
