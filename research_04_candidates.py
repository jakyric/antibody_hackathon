#!/usr/bin/env python
"""Part 4: enumerate the filter-passing design space, rank it, quantify the risk.

1. Validate the additive SPR model on data it never saw: the 460
   `no_reportable_kd` rows (weak/unmeasurable) vs the 190 quantified ones.
2. Enumerate every 10-mer within a restricted, evidence-backed alphabet,
   apply the exact challenge filter, rank by predicted log10 KD.
3. Bootstrap the model to get per-candidate uncertainty and report how many
   SPR observations back each substitution used.
4. Flag developability liabilities.
"""
import csv, zipfile, io, math, itertools, random
from collections import Counter, defaultdict
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

WT = "WGGDGFYAMD"
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

def edwt(s):
    return sum(x != y for x, y in zip(s, WT))

def passes(s):
    if edwt(s) < 3 or s in provided10:
        return False
    for i in range(10):
        for a in AAs:
            if a != s[i] and s[:i] + a + s[i+1:] in provided10:
                return False
    return True

# ----------------------------------------------------- fit the model
agg = defaultdict(list)
for d in spr:
    if (d["cdrh3"] != OTHER and len(d["cdrh3"]) == 10
            and d["binding_result"] == "quantified_binding" and d["kd_nM"]):
        agg[d["cdrh3"]].append(float(d["kd_nM"]))
seqs = sorted(agg)
y = np.array([np.mean(np.log10(agg[s])) for s in seqs])
X = onehot(seqs)
model = Ridge(alpha=10).fit(X, y)
nobs = np.zeros((10, 20))
for s in seqs:
    for i, a in enumerate(s):
        nobs[i, AI[a]] += 1

def pred(ss):
    return model.predict(onehot(ss))

# ----------------------------------------------- 1. out-of-distribution check
nork = sorted({d["cdrh3"] for d in spr
               if d["binding_result"] == "no_reportable_kd"
               and len(d["cdrh3"]) == 10 and d["cdrh3"] != OTHER})
print("== validation on data the model never saw ==")
print(f"   quantified (trained on): n={len(seqs)}   "
      f"no_reportable_kd (unseen): n={len(nork)}")
p_q, p_n = pred(seqs), pred(nork)
print(f"   predicted log10 KD: quantified mean {p_q.mean():.2f}, "
      f"no_reportable mean {p_n.mean():.2f}")
lab = np.r_[np.zeros(len(seqs)), np.ones(len(nork))]
sc = np.r_[p_q, p_n]
print(f"   AUC separating no_reportable from quantified: {roc_auc_score(lab, sc):.3f}")
print("   (the model is trained only on KD values, yet it ranks the")
print("    unmeasurable variants as weaker -- the direction generalises.")
print("    Caveat from the data README: no_reportable is NOT confirmed non-binding.)")

# ------------------------------------- 2. evidence-backed alphabet per position
print("\n== alphabet for enumeration ==")
C = model.coef_.reshape(10, 20)
alpha = {}
for i in range(10):
    cand = [(a, C[i, AI[a]], int(nobs[i, AI[a]])) for a in AAs
            if nobs[i, AI[a]] >= 4]
    # keep WT plus any residue no worse than WT + 0.25 log10 units
    wtc = C[i, AI[WT[i]]]
    keep = [a for a, c, n_ in cand if c <= wtc + 0.25]
    if WT[i] not in keep:
        keep.append(WT[i])
    alpha[i] = sorted(keep, key=lambda a: C[i, AI[a]])
    print(f"   pos {i+1} (WT {WT[i]}): " +
          " ".join(f"{a}{'*' if a==WT[i] else ''}[{C[i,AI[a]]:+.2f},n={int(nobs[i,AI[a]])}]"
                   for a in alpha[i]))
space = 1
for i in range(10):
    space *= len(alpha[i])
print(f"   enumeration space: {space:,} sequences")

# ----------------------------------------------- 3. enumerate, filter, rank
all_seqs = ["".join(t) for t in itertools.product(*[alpha[i] for i in range(10)])]
print(f"\n   enumerated {len(all_seqs):,}")
ed_ok = [s for s in all_seqs if edwt(s) >= 3]
print(f"   with >=3 changes from WT (filter b): {len(ed_ok):,}")
novel = [s for s in ed_ok if passes(s)]
print(f"   also >=2 from every provided sequence (filter a): {len(novel):,}")
blocked = len(ed_ok) - len(novel)
print(f"   -> the screen blocks {100*blocked/len(ed_ok):.1f}% of this space")

p = pred(novel)
order = np.argsort(p)

# bootstrap uncertainty
B = 300
rng = np.random.default_rng(0)
boot = np.zeros((B, len(novel)))
Xn = onehot(novel)
for b in range(B):
    idx = rng.integers(0, len(seqs), len(seqs))
    mb = Ridge(alpha=10).fit(X[idx], y[idx])
    boot[b] = mb.predict(Xn)
lo, hi_ = np.percentile(boot, [10, 90], axis=0)

# support: min number of SPR observations behind the substitutions used
def support(s):
    return min(int(nobs[i, AI[s[i]]]) for i in range(10) if s[i] != WT[i])

LIAB = [("NG deamidation", "NG"), ("NS deamidation", "NS"), ("DG isomerisation", "DG"),
        ("DS isomerisation", "DS"), ("DD", "DD"), ("free Cys", "C"),
        ("Met oxidation", "M"), ("Trp oxidation", "W")]

def liabilities(s):
    out = []
    for nm, m in LIAB:
        if m in s and not (m in WT and s.count(m) <= WT.count(m)):
            out.append(nm)
    if "N" in s:
        for i in range(len(s)-2):
            if s[i] == "N" and s[i+1] != "P" and s[i+2] in "ST":
                out.append("N-glyc NxS/T")
    return out

print("\n== top 25 filter-passing candidates by predicted affinity ==")
print(f"{'rank':>4} {'cdrh3':<12}{'ed':>3}{'pred_KD_nM':>12}{'80% CI (nM)':>22}"
      f"{'min_n':>7}  mutations / liabilities")
for r, j in enumerate(order[:25], 1):
    s = novel[j]
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i])
    lb = ";".join(liabilities(s)) or "-"
    print(f"{r:>4} {s:<12}{edwt(s):>3}{10**p[j]:>12.3f}"
          f"{f'{10**lo[j]:.3f} - {10**hi_[j]:.2f}':>22}{support(s):>7}  {muts}  [{lb}]")

print("\n== where the measured reference points sit on this scale ==")
refs = ["WGGDGFYAMD", "WGSQGFYAYD", "WGSDGWYAMD", "WGGQGFYAMD", "WGGSGPYLFD"]
for s in refs:
    meas = np.mean(agg[s]) if s in agg else None
    print(f"   {s}  ed {edwt(s)}  measured {meas if meas is None else f'{meas:.3f}'} nM"
          f"   predicted {10**pred([s])[0]:.3f} nM   passes filter: {passes(s)}")

# ----------------------------- 4. how novel is the top region, really?
print("\n== sanity: is the top-ranked region unexplored, or just barely dodged? ==")
for j in order[:8]:
    s = novel[j]
    d2 = []
    for t in provided10:
        dd = sum(a != b for a, b in zip(s, t))
        if dd <= 3:
            d2.append((dd, t, facs.get(t, "SPR-only")))
    d2.sort()
    near = "; ".join(f"d{d}:{t}({c})" for d, t, c in d2[:4])
    print(f"   {s}: {len(d2)} provided seqs within 3 subs | {near}")

# ----------------------------- 5. diversity-aware pick
print("\n== a diverse shortlist (greedy, >=3 apart from each other) ==")
picked = []
for j in order:
    s = novel[j]
    if all(sum(a != b for a, b in zip(s, t)) >= 3 for t in picked):
        picked.append(s)
    if len(picked) == 12:
        break
for r, s in enumerate(picked, 1):
    jj = novel.index(s)
    muts = ",".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i])
    print(f"   {r:>2} {s}  ed {edwt(s)}  pred {10**p[jj]:>8.3f} nM  "
          f"min_n {support(s):>3}  {muts}")

with open("/home/hackathon/jakob/candidates_ranked.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["cdrh3", "ed_wt", "pred_kd_nM", "ci10_nM", "ci90_nM",
                "min_spr_support", "mutations", "liabilities"])
    for j in order:
        s = novel[j]
        w.writerow([s, edwt(s), f"{10**p[j]:.4f}", f"{10**lo[j]:.4f}",
                    f"{10**hi_[j]:.4f}", support(s),
                    ";".join(f"{WT[i]}{i+1}{s[i]}" for i in range(10) if s[i] != WT[i]),
                    ";".join(liabilities(s))])
print(f"\nwrote {len(novel)} ranked candidates to ~/jakob/candidates_ranked.csv")
