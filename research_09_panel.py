#!/usr/bin/env python
"""Part 9: the confound check that decides the design rule, then the panel.

Confound: in the FACS library WT M appears at position 9 in only 1.2% of
sequences, and 148/180 quantified SPR variants mutate position 9. So
"position 9 == M" is nearly a proxy for "close to wild type", and an additive
ridge can load the whole near-WT affinity advantage onto M9M.
Yet 8 of the 10 sub-nanomolar ed>=3 variants DO mutate position 9.
Refit on ed>=3 only to break the confound.
"""
import csv, zipfile, io, itertools, math
from collections import defaultdict
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from scipy.stats import spearmanr

WT = "WGGDGFYAMD"; CHAMP = "WGSQGFYAYD"
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
seqs = sorted(agg); y = np.array([np.mean(np.log10(agg[s])) for s in seqs])
X = onehot(seqs); edv = np.array([ed(s) for s in seqs])

print("== confound check: residue frequencies at each position in the training set ==")
print(f"{'pos':>4}{'WT':>4}{'n_WT':>7}{'frac_WT':>9}   corr(is_WT_at_pos, ed)")
for i in range(10):
    isw = np.array([1.0 if s[i] == WT[i] else 0.0 for s in seqs])
    print(f"{i+1:>4}{WT[i]:>4}{int(isw.sum()):>7}{isw.mean():>9.2f}"
          f"{spearmanr(isw, edv)[0]:>24.3f}")

print("\n== refit using only variants with ed>=3 (the filter-relevant regime) ==")
m3 = edv >= 3
print(f"   n = {m3.sum()}")
for lbl, mask in (("all 180", np.ones(len(seqs), bool)), ("ed>=3 only", m3)):
    mod = Ridge(alpha=10).fit(X[mask], y[mask])
    C = mod.coef_.reshape(10, 20)
    nobs = np.zeros((10, 20))
    for s in np.array(seqs)[mask]:
        for i, a in enumerate(s): nobs[i, AI[a]] += 1
    print(f"\n   --- fitted on {lbl} ---")
    print(f"   {'pos':>4}{'WT':>4}{'WT_coef':>10}  best 4 observed (n>=4), coef (lower=tighter)")
    for i in range(10):
        obs = sorted([(a, C[i, AI[a]], int(nobs[i, AI[a]])) for a in AAs
                      if nobs[i, AI[a]] >= 4], key=lambda t: t[1])
        wtc = C[i, AI[WT[i]]] if nobs[i, AI[WT[i]]] >= 1 else float("nan")
        print(f"   {i+1:>4}{WT[i]:>4}{wtc:>+10.2f}  " +
              "  ".join(f"{a}{'*' if a==WT[i] else ''}{c:+.2f}(n{n})" for a, c, n in obs[:4]))
    # CV on ed>=3
    if lbl == "ed>=3 only":
        pr = np.zeros(mask.sum())
        Xs, ys = X[mask], y[mask]
        for a, b in KFold(5, shuffle=True, random_state=0).split(Xs):
            pr[b] = Ridge(alpha=10).fit(Xs[a], ys[a]).predict(Xs[b])
        print(f"\n   5-fold CV Spearman within ed>=3 only: {spearmanr(pr, ys)[0]:+.3f}")
        print("   (vs +0.71 on the full set, which is inflated by the easy")
        print("    near-WT-vs-far-from-WT contrast the filter forbids us to exploit)")

print("\n== position 9 specifically ==")
p9 = defaultdict(list)
for s in seqs:
    p9[s[8]].append(10**y[seqs.index(s)] if False else np.mean(agg[s]))
for a in sorted(p9, key=lambda a: np.median(p9[a])):
    v = np.array(p9[a])
    if len(v) >= 4:
        print(f"   pos9={a}{'*' if a=='M' else ''}  n={len(v):>3}  median KD "
              f"{np.median(v):>8.1f} nM  best {v.min():>7.3f} nM")
sub = [s for s in seqs if np.mean(agg[s]) < 1 and ed(s) >= 3]
print(f"   among the {len(sub)} sub-nM variants with ed>=3: "
      f"{sum(1 for s in sub if s[8] != 'M')} mutate position 9")
print("   -> position 9 is NOT off limits. The additive penalty on mutating it")
print("      is largely the edit-distance confound.")

# --------------------------------------------------------------- final panel
print("\n\n" + "="*72)
print("FINAL PANEL: structure-constrained, evidence-weighted, filter-passing")
print("="*72)
# design rule from part 8 (structure) + part 9 (confound):
#   FREEZE positions 1 (W), 5 (G), 6 (F), 7 (Y)  -- the only HER2 contacts
#   VARY  positions 2, 3, 4, 8, 9, 10            -- zero HER2 contacts
FREEZE = {0: "W", 4: "G", 5: "F", 6: "Y"}
VARY = {
    1: "GYWTA",       # pos2: WT + measured-good
    2: "SGVDAY",      # pos3: S strongly good
    3: "QESYDN",      # pos4: Q/E best measured
    7: "ELASVMY",     # pos8: E/L best
    8: "MYFLWH",      # pos9: WT M + the residues in sub-nM variants
    9: "DIEGTAV",     # pos10: I/D/E
}
alpha = {}
for i in range(10):
    alpha[i] = FREEZE[i] if i in FREEZE else VARY[i]
space = 1
for i in range(10): space *= len(alpha[i])
print(f"\nconstrained space: {space:,} sequences "
      f"(positions 1,5,6,7 frozen to WT = the HER2 contacts)")

cand = ["".join(t) for t in itertools.product(*[alpha[i] for i in range(10)])]
cand = [s for s in cand if ed(s) >= 3]
print(f"  with ed>=3: {len(cand):,}")
cand = [s for s in cand if passes(s)]
print(f"  passing the full challenge filter: {len(cand):,}")

# two scorers: fitted on all data, and fitted on ed>=3 only; average the ranks
mA = Ridge(alpha=10).fit(X, y)
mB = Ridge(alpha=10).fit(X[m3], y[m3])
Xc = onehot(cand)
pA, pB = mA.predict(Xc), mB.predict(Xc)
rA, rB = pA.argsort().argsort(), pB.argsort().argsort()
score = (rA + rB) / 2.0
# similarity to the measured sub-nM set (an empirical prior)
subnm = [s for s in seqs if np.mean(agg[s]) < 1.0]
def minsim(s): return min(ed(s, t) for t in subnm)
LIAB = ("NG", "NS", "DG", "DS", "NN")
def liab(s):
    out = [m for m in LIAB if m in s and m not in WT]
    if "C" in s: out.append("freeCys")
    for i in range(8):
        if s[i] == "N" and s[i+1] != "P" and s[i+2] in "ST": out.append("Nglyc")
    return out
order = np.argsort(score)
print("\n== top 30 by consensus rank of both scorers ==")
print(f"{'#':>3} {'cdrh3':<12}{'ed':>3}{'predA_nM':>10}{'predB_nM':>10}"
      f"{'d_to_subnM':>12}  mutations  [liabilities]")
rows = []
for r, j in enumerate(order[:30], 1):
    s = cand[j]
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i])
    lb = ";".join(liab(s)) or "-"
    print(f"{r:>3} {s:<12}{ed(s):>3}{10**pA[j]:>10.2f}{10**pB[j]:>10.2f}"
          f"{minsim(s):>12}  {muts}  [{lb}]")

# diverse final 24: >=3 apart, no liabilities, spread over ed 3-6
print("\n== recommended 24-design panel (>=3 apart from each other, "
      "no liability motifs) ==")
picked = []
for j in order:
    s = cand[j]
    if liab(s): continue
    if all(ed(s, t) >= 3 for t in picked):
        picked.append(s)
    if len(picked) == 24: break
print(f"{'#':>3} {'cdrh3':<12}{'ed_wt':>6}{'predA_nM':>10}{'predB_nM':>10}"
      f"{'d_subnM':>9}  mutations")
for r, s in enumerate(picked, 1):
    j = cand.index(s)
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i])
    print(f"{r:>3} {s:<12}{ed(s):>6}{10**pA[j]:>10.2f}{10**pB[j]:>10.2f}"
          f"{minsim(s):>9}  {muts}")

with open("/home/hackathon/jakob/panel_final.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["rank", "cdrh3", "ed_wt", "pred_kd_nM_allfit", "pred_kd_nM_ed3fit",
                "min_dist_to_measured_subnM", "mutations", "liabilities",
                "vh_sequence"])
    VH = [d["vh_sequence"] for d in spr if d["cdrh3"] == WT][0]
    for r, s in enumerate(picked, 1):
        j = cand.index(s)
        vh = VH[:98] + s + VH[108:]
        w.writerow([r, s, ed(s), f"{10**pA[j]:.3f}", f"{10**pB[j]:.3f}",
                    minsim(s),
                    ";".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i]),
                    ";".join(liab(s)), vh])
print("\nwrote ~/jakob/panel_final.csv (24 designs, with full VH sequences)")
with open("/home/hackathon/jakob/candidates_constrained_all.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["cdrh3","ed_wt","pred_kd_nM_allfit","pred_kd_nM_ed3fit","consensus_rank","mutations"])
    for j in order:
        s = cand[j]
        w.writerow([s, ed(s), f"{10**pA[j]:.3f}", f"{10**pB[j]:.3f}", f"{score[j]:.1f}",
                    ";".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i]!=WT[i])])
print(f"wrote ~/jakob/candidates_constrained_all.csv ({len(cand)} rows)")
