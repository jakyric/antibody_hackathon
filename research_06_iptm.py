#!/usr/bin/env python
"""Part 6: how much should we trust folding-model confidence to rank designs?

The docs (18-vm.md) say: rank HER2 designs on min_ipae / iptm / binder_plddt.
The Smorodina et al. data shipped with the hackathon measures exactly how well
interface confidence separates a REAL antibody-antigen pair from a shuffled one.
Diagonal = real pair, off-diagonal = shuffled.
"""
import zipfile, io, csv, os
import numpy as np
from sklearn.metrics import roc_auc_score

Z = "/home/hackathon/data/hackathon-data/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip"

for model in ("af3", "boltz2", "chai1"):
    with zipfile.ZipFile(Z) as z:
        with z.open(f"iptm_confidence_scores/{model}_matrix_clean.csv") as fh:
            txt = io.TextIOWrapper(fh, "utf-8")
            rdr = csv.reader(txt)
            hdr = next(rdr)
            rows, labels = [], []
            for r in rdr:
                labels.append(r[0])
                rows.append([float(x) if x not in ("", "NA", "nan") else np.nan
                             for x in r[1:]])
    cols = hdr[1:]
    M = np.array(rows)
    print(f"\n{'='*70}\n{model.upper()}: matrix {M.shape}  "
          f"(rows=nanobody, cols=antigen), {np.isnan(M).sum()} missing")

    # real pairs: the diagonal, matched by name
    real, shuf = [], []
    for i, rl in enumerate(labels):
        for j, cl in enumerate(cols):
            v = M[i, j]
            if np.isnan(v):
                continue
            (real if (i == j) else shuf).append(v)
    real, shuf = np.array(real), np.array(shuf)
    print(f"  real (diagonal) n={len(real)}: ipTM mean {real.mean():.3f} "
          f"median {np.median(real):.3f}  IQR {np.percentile(real,25):.3f}-{np.percentile(real,75):.3f}")
    print(f"  shuffled       n={len(shuf)}: ipTM mean {shuf.mean():.3f} "
          f"median {np.median(shuf):.3f}  IQR {np.percentile(shuf,25):.3f}-{np.percentile(shuf,75):.3f}")
    lab = np.r_[np.ones(len(real)), np.zeros(len(shuf))]
    sc = np.r_[real, shuf]
    auc = roc_auc_score(lab, sc)
    print(f"  AUC(real vs shuffled by ipTM) = {auc:.3f}")

    # the operationally relevant question: if I keep the top-N by ipTM, how
    # many are real binders?
    print("  if you rank ALL pairs by ipTM and keep the top N, precision is:")
    o = np.argsort(-sc)
    for N in (10, 25, 50, 100, 200):
        if N <= len(sc):
            prec = lab[o[:N]].mean()
            # base rate
            print(f"    top {N:>4}: {100*prec:>5.1f}% real   "
                  f"(base rate {100*lab.mean():.2f}%, enrichment {prec/lab.mean():.0f}x)")
    # how many shuffled pairs look "confident"?
    for thr in (0.6, 0.7, 0.8, 0.9):
        print(f"    ipTM > {thr}: {100*(shuf>thr).mean():>5.1f}% of shuffled pairs, "
              f"{100*(real>thr).mean():>5.1f}% of real pairs")
    # what threshold would you need for 90% precision?
    for want in (0.5, 0.9):
        best = None
        for thr in np.arange(0.3, 0.99, 0.01):
            sel = sc > thr
            if sel.sum() >= 5 and lab[sel].mean() >= want:
                best = (thr, sel.sum(), lab[sel].mean()); break
        print(f"    threshold for >={100*want:.0f}% precision: "
              f"{'none exists' if best is None else f'ipTM>{best[0]:.2f} keeps {best[1]} pairs, {100*best[2]:.0f}% real'}")
