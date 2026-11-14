"""Train simple CDRH3 binding baselines on the trastuzumab FACS data and test them.

Training/eval protocol
  * Train on the authors' clonotype-split train set, choose settings on their val set, report on
    their test set (clusters of similar CDRH3s never straddle the split).
  * Compare with (a) the authors' CNN (their weights, numpy re-implementation) and (b) a trivial
    "closer to trastuzumab is better" score.
  * External checks: the SPR table (sequences not present in the FACS table) and the Absci
    trastuzumab designs (Shanehsazzadeh et al.).
  * Finally refit the best simple model on all 524,346 FACS variants and save it for the gate.

Run:  baseline/.venv/bin/python baseline/train_eval.py
"""
import json
import time

import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, log_loss, roc_auc_score

import common as C

SEED = 0
C.OUT.mkdir(exist_ok=True)
t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:6.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------ data
tr = pd.read_csv(C.MEHTA / "train.csv")
va = pd.read_csv(C.MEHTA / "val.csv")
te = pd.read_csv(C.MEHTA / "test.csv")
Xtr, Xva, Xte = (C.encode_idx(d.seq) for d in (tr, va, te))
log(f"authors' clonotype split: train {len(tr):,} / val {len(va):,} / test {len(te):,}; "
    f"positives (high) {100 * tr.label.mean():.1f}%")

# how novel is the test set relative to train (distance-<=1 neighbours)?
train_index = C.ProvidedData(tr.seq)
nb = train_index.has_neighbour_within_1(Xte)
log(f"test sequences with a train sequence within 1 mutation: {100 * nb.mean():.1f}% "
    f"(the hackathon novelty rule requires 0%)")


# ------------------------------------------------------------------ features / models
def pair_features(idx):
    """Sparse one-hot of positions (200) plus all pairwise position interactions (45 x 400)."""
    n = len(idx)
    cols = [idx[:, p].astype(np.int64) + 20 * p for p in range(C.L)]
    k = 200
    for p in range(C.L):
        for q in range(p + 1, C.L):
            cols.append(k + idx[:, p].astype(np.int64) * 20 + idx[:, q])
            k += 400
    cols = np.stack(cols, 1)
    rows = np.repeat(np.arange(n), cols.shape[1])
    return sp.csr_matrix((np.ones(cols.size, np.float32), (rows, cols.ravel())), shape=(n, k))


def onehot_flat(idx):
    return C.onehot(idx).reshape(len(idx), -1)


def metrics(y, p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return {
        "auc": float(roc_auc_score(y, p)),
        "ap": float(average_precision_score(y, p)),
        "acc@0.5": float(accuracy_score(y, p > 0.5)),
        "logloss": float(log_loss(y, p)),
    }


results = {}
scorers = {}  # name -> function(idx) -> P(high)

log("model A: logistic regression on one-hot positions (additive)")
lr = LogisticRegression(C=1.0, max_iter=500).fit(onehot_flat(Xtr), tr.label)
scorers["LR one-hot"] = lambda idx, m=lr: m.predict_proba(onehot_flat(idx))[:, 1]

log("model B: logistic regression with pairwise position interactions")
lrp = LogisticRegression(C=0.5, max_iter=300).fit(pair_features(Xtr), tr.label)
scorers["LR pairwise"] = lambda idx, m=lrp: m.predict_proba(pair_features(idx))[:, 1]

log("model C: gradient-boosted trees, positions as categorical features")
hgb_kw = dict(learning_rate=0.1, max_leaf_nodes=63, min_samples_leaf=40, early_stopping=False,
              categorical_features=[True] * C.L, random_state=SEED)
probe = HistGradientBoostingClassifier(max_iter=400, **hgb_kw).fit(Xtr.astype(np.float32), tr.label)
best_auc, best_iter = -1, 0
for i, pv in enumerate(probe.staged_predict_proba(Xva.astype(np.float32)), 1):
    if i % 20 == 0:
        a = roc_auc_score(va.label, pv[:, 1])
        if a > best_auc:
            best_auc, best_iter = a, i
log(f"   best iteration on val: {best_iter} (val AUC {best_auc:.4f})")
hgb = HistGradientBoostingClassifier(max_iter=best_iter, **hgb_kw).fit(Xtr.astype(np.float32), tr.label)
scorers["HGB trees"] = lambda idx, m=hgb: m.predict_proba(idx.astype(np.float32))[:, 1]

cnn = C.AuthorsCNN(C.MEHTA / "saved_weights.h5")
scorers["Authors' CNN"] = cnn.predict_proba
scorers["trivial: -mutations vs WT"] = lambda idx: -C.hamming_to_wt(idx).astype(float)

# ------------------------------------------------------------------ test on the authors' split
log("evaluating on the authors' clonotype test set")
bins = [(1, 3), (4, 5), (6, 7), (8, 10)]
results["test_split"] = {}
strict = ~nb  # test sequences with no train sequence within 1 mutation (what the novelty rule demands)
log(f"   strict-novel test subset: {int(strict.sum()):,} sequences, {100 * te.label[strict].mean():.1f}% positives")
print(f"\n{'model':28s} {'AUC':>7s} {'AP':>7s} {'acc':>7s} {'logloss':>8s} {'AUC strict':>10s} | AUC by mutations vs WT: " +
      "  ".join(f"{a}-{b}" for a, b in bins))
for name, f in scorers.items():
    p = f(Xte)
    if name.startswith("trivial"):
        m = {"auc": float(roc_auc_score(te.label, p)), "ap": float(average_precision_score(te.label, p))}
    else:
        m = metrics(te.label, p)
    per_bin = {}
    for a, b in bins:
        sel = (te.edit_distance >= a) & (te.edit_distance <= b)
        if te.label[sel].nunique() == 2:
            per_bin[f"{a}-{b}"] = float(roc_auc_score(te.label[sel], p[sel]))
    m["auc_by_distance"] = per_bin
    m["auc_strict_novel"] = float(roc_auc_score(te.label[strict], p[strict]))
    results["test_split"][name] = m
    print(f"{name:28s} {m['auc']:7.4f} {m['ap']:7.4f} {m.get('acc@0.5', float('nan')):7.4f} "
          f"{m.get('logloss', float('nan')):8.4f} {m['auc_strict_novel']:10.4f} |                        " +
          "  ".join(f"{per_bin.get(f'{a}-{b}', float('nan')):.3f}" for a, b in bins))

# ------------------------------------------------------------------ external check 1: SPR
log("\nexternal check: SPR table (only sequences NOT in the FACS table)")
facs = C.load_facs()
facs_set = set(facs.cdrh3)
spr_all = C.load_spr()
wt_kd = spr_all.loc[spr_all.cdrh3 == C.WT, "kd_nM"].dropna()
wt_worst = float(wt_kd.max())
spr = spr_all[(spr_all.cdrh3.str.len() == C.L) & (spr_all.cdrh3 != C.WT) & spr_all.cdrh3.str.isalpha()].copy()
n_overlap = int(spr.cdrh3.isin(facs_set).sum())
spr = spr[~spr.cdrh3.isin(facs_set)]
g = spr.groupby("cdrh3")
agg = pd.DataFrame({
    "kd": g["kd_nM"].median(),
    "status": g["binding_result"].agg(
        lambda s: "binder" if (s == "quantified_binding").any()
        else ("no_kd" if (s == "no_reportable_kd").any() else "other")),
}).reset_index()
Xs = C.encode_idx(agg.cdrh3)
agg["hamming"] = C.hamming_to_wt(Xs)
log(f"   {len(agg)} unique SPR CDRH3s (dropped {n_overlap} rows overlapping FACS): "
    f"{(agg.status == 'binder').sum()} with a KD, {(agg.status == 'no_kd').sum()} with no reportable KD; "
    f"trastuzumab reference KD up to {wt_worst:.3f} nM")
binders = agg[agg.status == "binder"].copy()
top = binders.kd <= wt_worst
log(f"   variants with KD <= trastuzumab's worst reference: {int(top.sum())} of {len(binders)}")
results["spr"] = {"n_unique": int(len(agg)), "n_binders": int((agg.status == "binder").sum()),
                  "n_no_kd": int((agg.status == "no_kd").sum()), "n_beating_wt": int(top.sum())}
b_idx = C.encode_idx(binders.cdrh3)
sel_cls = agg.status.isin(["binder", "no_kd"])
y_cls = (agg.status[sel_cls] == "binder").astype(int)
print(f"\n{'model':28s} {'AUC binder vs no-KD':>20s} {'Spearman vs -KD':>16s} {'Spearman, 3-6 mut.':>19s} {'best-variants in top 20%':>25s}")
for name, f in scorers.items():
    p_all = f(Xs)
    auc = roc_auc_score(y_cls, p_all[sel_cls.values])
    pb = f(b_idx)
    rho = spearmanr(pb, -binders.kd.values).statistic
    mid = ((binders.hamming >= 3) & (binders.hamming <= 6)).values
    rho_mid = spearmanr(pb[mid], -binders.kd.values[mid]).statistic if mid.sum() > 5 else float("nan")
    # rank percentile of the variants that beat trastuzumab (1.0 = best of all binders)
    order = (-pb).argsort().argsort() / max(len(pb) - 1, 1)
    in_top = int((order[top.values] <= 0.2).sum())
    results["spr"][name] = {"auc_binder_vs_nokd": float(auc), "spearman_vs_neg_kd": float(rho),
                            "spearman_3to6_mutations": float(rho_mid), "beating_wt_in_top20pct": in_top}
    print(f"{name:28s} {auc:20.3f} {rho:16.3f} {rho_mid:19.3f} {in_top:>14d} of {int(top.sum())}")

# ------------------------------------------------------------------ external check 2: Absci designs
log("\nexternal check: Absci trastuzumab designs (Shanehsazzadeh et al.), all measured binders")
ab = pd.read_csv(C.MEHTA / "shanehsazzadeh_test.csv")
ab = ab[ab.seq.str.len() == C.L]
ab = ab[~ab.seq.isin(facs_set)]
Xa = C.encode_idx(ab.seq)
print(f"   {len(ab)} sequences not in the FACS table; affinity column range {ab.affinity.min():.1f}-{ab.affinity.max():.1f}")
results["absci"] = {"n": int(len(ab))}
for name, f in scorers.items():
    if name.startswith("trivial") or len(ab) < 5:
        continue
    p = f(Xa)
    results["absci"][name] = {"frac_predicted_binder": float((p > 0.5).mean()),
                              "spearman_vs_affinity": float(spearmanr(p, ab.affinity).statistic)}
    print(f"   {name:28s} predicted binder (p>0.5): {100 * (p > 0.5).mean():5.1f}%   Spearman vs affinity: "
          f"{results['absci'][name]['spearman_vs_affinity']:.3f}")

# ------------------------------------------------------------------ final model on all FACS data
log("\nrefitting the HGB model on all 524,346 FACS variants for the gate")
Xall = C.encode_idx(facs.cdrh3)
final = HistGradientBoostingClassifier(max_iter=best_iter, **hgb_kw).fit(Xall.astype(np.float32), facs.label)
joblib.dump({"model": final, "alphabet": C.AA, "trained_on": "all 524,346 FACS variants, label=high",
             "best_iter": best_iter}, C.OUT / "hgb_full.joblib")
(C.OUT / "metrics.json").write_text(json.dumps(results, indent=2))
log(f"saved {C.OUT / 'hgb_full.joblib'} and metrics.json")
