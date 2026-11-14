"""Paper-faithful BLOSUM library design (Chinery et al., bioRxiv 10.1101/2024.03.26.586756), phase 1.

What the paper did, and what is done here
  1. Generate 1,000,000 CDRH3 variants by sampling every one of the ten CDRH3 positions independently from
     P(b | a) proportional to f_b * exp(0.25 * s_ab) (BLOSUM-45, SAbDab CDRH3 background, 'stay' included).  [same]
  2. Remove redundant sequences.                                                                              [same]
  3. Keep sequences with CNN P(binder) above 90% (CNN trained on the FACS screen).         [authors' CNN weights]
  4. Remove Cys, 'GPR' motifs, sequences with 4 or more liability motifs and sequences with 5 or more of one
     residue. The paper counted motifs over the Fv; the rest of the Fv is trastuzumab's own, so we count in CDRH3. [adapted]
  5. Model with an antibody structure predictor and drop designs whose CDRH3 RMSD to the wild-type model exceeds
     3.5 A or whose total CDR RMSD exceeds 1.5 A. The paper used ABodyBuilder2; we use ABodyBuilder3.  [phase 2]
  6. Sample test designs equally across edit distances (the paper used 1-10).                                  [same]
Added for the hackathon: designs must be more than 2 mutations from every provided FACS/SPR sequence (abguard's
novelty gate), and we use edit distances 4-10 (4 satisfies both readings of the CDRH3-identity rule).

Run: baseline/.venv/bin/python baseline/paper_pipeline.py [--n 1000000] [--per-bucket 4] [--seed 0]
"""
import argparse
import json
import math
import re

import blosum as bl
import joblib
import numpy as np
import pandas as pd

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=1_000_000)
ap.add_argument("--per-bucket", type=int, default=4, help="designs per edit-distance bucket sent to structure modelling")
ap.add_argument("--min-ed", type=int, default=4)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
wt = C.encode_idx([C.WT])[0]
funnel = {}

# ------------------------------------------------------------------ 1. BLOSUM sampling (paper recipe)
H3_COUNTS = {"A": 6710, "C": 555, "D": 6442, "E": 1636, "F": 3404, "G": 6582, "H": 1081, "I": 1342, "K": 970,
             "L": 2509, "M": 1218, "N": 1658, "P": 2090, "Q": 767, "R": 5196, "S": 4350, "T": 2835, "V": 2975,
             "W": 1686, "Y": 9028}
tot = sum(H3_COUNTS.values())
bg = np.array([H3_COUNTS[x] / tot for x in C.AA])
B45 = bl.BLOSUM(45)
LAMBDA = 0.25
S = np.array([[B45[x][y] for y in C.AA] for x in C.AA], dtype=float)
P = bg[None, :] * np.exp(LAMBDA * S)  # f_b * exp(lambda s_ab); rows renormalised below
P = P / P.sum(1, keepdims=True)  # P[a, b], 'stay' (b == a) included, as in the paper's sampling
pos_cdf = np.cumsum(P[wt], axis=1)  # (10, 20): distribution at each position given the wild-type residue
print("P(stay) per position:", {f"{C.WT[i]}{i + 1}": round(float(P[wt[i], wt[i]]), 3) for i in range(10)})

draws = []
CH = 250_000
for s in range(0, a.n, CH):
    m = min(CH, a.n - s)
    u = rng.random((m, 10))
    draws.append((u[:, :, None] > pos_cdf[None, :, :]).sum(-1).astype(np.int8))
idx = np.concatenate(draws)
funnel["generated"] = int(len(idx))

# ------------------------------------------------------------------ 2. dedupe, edit distance
codes = (idx.astype(np.int64) * C._PW).sum(1)
_, first = np.unique(codes, return_index=True)
idx = idx[np.sort(first)]
ed = (idx != wt[None, :]).sum(1)
funnel["unique"] = int(len(idx))
facs = C.load_facs()
facs_ed = np.bincount(facs.edit_distance, minlength=11)[1:11] / len(facs)
gen_ed = np.bincount(ed, minlength=11)[:11]
print("\nedit distance      " + " ".join(f"{k:>7d}" for k in range(0, 11)))
print("unique samples     " + " ".join(f"{v:>7d}" for v in gen_ed))
print("FACS library (%)   " + " ".join(f"{100 * x:7.1f}" for x in [0] + list(facs_ed)))

# ------------------------------------------------------------------ 3. CNN binder probability above 90%
cnn = C.AuthorsCNN(C.MEHTA / "saved_weights.h5")
p_cnn = cnn.predict_proba(idx)
keep = p_cnn > 0.9
funnel["cnn_over_0.9"] = int(keep.sum())
print(f"\nCNN P(binder) > 0.9: {int(keep.sum()):,} of {len(idx):,}")

# ------------------------------------------------------------------ 4. liabilities and composition
seqs = np.array(["".join(C.AA[i] for i in r) for r in idx[keep]])
LIAB = ["NG", "NS", "NT", "NN", "DG", "DS", "DD", "DT", "DH", "DN", "DP", "DQ", "KE", "KD", "EK", "ED", "M", "W"]


def n_liab(s):
    n = sum(len(re.findall(f"(?={r})", s)) for r in LIAB)
    return n + len(re.findall(r"N[^P][ST]", s)) + sum(s.count(x) for x in ("RGD", "RYD", "LDV"))


liab = np.array([n_liab(s) for s in seqs])
no_cys = ~np.char.find(seqs, "C").astype(bool) | (np.char.find(seqs, "C") < 0)
no_gpr = np.char.find(seqs, "GPR") < 0
few_liab = liab < 4
max_rep = np.array([max(s.count(x) for x in set(s)) for s in seqs])
rep_ok = max_rep < 5
sel = no_cys & no_gpr & few_liab & rep_ok
print(f"after Cys/GPR/liabilities/repeats filters: {int(sel.sum()):,}  "
      f"(Cys {int((~no_cys).sum()):,}, GPR {int((~no_gpr).sum()):,}, 4+ liabilities {int((~few_liab).sum()):,}, 5+ repeats {int((~rep_ok).sum()):,} removed)")
funnel["after_sequence_filters"] = int(sel.sum())
idx2, ed2, pc2, liab2, seqs2 = idx[keep][sel], ed[keep][sel], p_cnn[keep][sel], liab[sel], seqs[sel]

# ------------------------------------------------------------------ hackathon rules: novelty, edit distance
prov = C.ProvidedData(C.provided_sequences())
near2 = prov.has_neighbour_within_2(idx2)
funnel["not_within_2_of_provided"] = int((~near2).sum())
ok = (~near2) & (ed2 >= a.min_ed)
funnel[f"and_edit_distance_ge_{a.min_ed}"] = int(ok.sum())
print(f"hackathon novelty (more than 2 mutations from every FACS/SPR sequence): {int((~near2).sum()):,}; "
      f"and edit distance >= {a.min_ed}: {int(ok.sum()):,}")
pool = pd.DataFrame({"cdrh3": seqs2[ok], "edit_distance": ed2[ok], "p_cnn": pc2[ok], "liabilities": liab2[ok]})
pool_idx = idx2[ok]
pool["p_hgb"] = joblib.load(C.OUT / "hgb_full.joblib")["model"].predict_proba(pool_idx.astype(np.float32))[:, 1]
pool["p_mean"] = (pool.p_cnn + pool.p_hgb) / 2
pool["blosum_logp"] = np.log(P[wt[None, :], pool_idx]).sum(1)  # paper generative log-probability (stay included)
pool["changed_positions"] = [",".join(str(i + 1) for i in np.where(r != wt)[0]) for r in pool_idx]
print("pool by edit distance:", pool.edit_distance.value_counts().sort_index().to_dict())

# ------------------------------------------------------------------ 6. equal sampling across edit distances
rows = []
for k in sorted(pool.edit_distance.unique()):
    g = pool[pool.edit_distance == k]
    rows.append(g.sample(n=min(a.per_bucket, len(g)), random_state=a.seed))
mod = pd.concat(rows).reset_index(drop=True)
mod.insert(0, "id", [f"p{i + 1:02d}" for i in range(len(mod))])
mod.to_csv(C.OUT / "paper_modelling_set.csv", index=False)
funnel["sent_to_structure_modelling"] = int(len(mod))
(C.OUT / "paper_funnel.json").write_text(json.dumps(funnel, indent=2))
print(f"\nsampled {len(mod)} designs ({a.per_bucket} per edit-distance bucket) for structure modelling:")
print(mod[["id", "cdrh3", "edit_distance", "p_cnn", "liabilities", "changed_positions"]].round(3).to_string(index=False))
print("\nfunnel:", json.dumps(funnel))
