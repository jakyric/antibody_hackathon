"""Stage 6: ipSAE and iLIS from the PAE matrices of the designs that pass the contact pre-filter.

PAE matrices are large (~1.4 MB), so only designs that are docked (ipTM) and keep the crystal contacts are downloaded.
  ipSAE (Dunbrack 2025, d0res form, PAE cut-off 10 A): max over residues of the mean TM-like score over partner
         residues with PAE < cut-off; asymmetric values are combined by max.
  iLIS  (Kim et al.): sqrt(LIS * cLIS), LIS = mean of (12 - PAE)/12 over PAE < 12 A, cLIS the same restricted to
         residue pairs in contact (C-beta within 8 A). This is a re-implementation of the published definition.
Both are computed for the Fv (VH + VL) against HER2.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import out, save_summary
from .japanfold import JapanFold


def d0(n: int) -> float:
    return max(1.0, 1.24 * (max(n, 27) - 15) ** (1 / 3) - 1.8)


def ipsae(pae: np.ndarray, X: np.ndarray, Y: np.ndarray, cut: float = 10.0) -> float:
    best = 0.0
    for i in X:
        row = pae[i, Y]
        ok = row < cut
        n = int(ok.sum())
        if n:
            best = max(best, float((1.0 / (1.0 + (row[ok] / d0(n)) ** 2)).mean()))
    return best


def ilis(pae: np.ndarray, X: np.ndarray, Y: np.ndarray, xyz: np.ndarray, cut: float = 12.0, contact: float = 8.0):
    t = lambda blk: np.where(blk < cut, (cut - blk) / cut, 0.0)
    A, B = t(pae[np.ix_(X, Y)]), t(pae[np.ix_(Y, X)])
    vals = np.concatenate([A[A > 0], B[B > 0]])
    lis = float(vals.mean()) if len(vals) else 0.0
    c = np.linalg.norm(xyz[X][:, None, :] - xyz[Y][None, :, :], axis=2) <= contact
    cv = np.concatenate([A[c], B[c.T]])
    clis = float(cv.mean()) if len(cv) else 0.0
    return lis, clis, float(np.sqrt(lis * clis))


def cbeta_and_lengths(pdb) -> tuple[np.ndarray, dict]:
    res, order = {}, []
    for line in open(pdb):
        if line.startswith("ATOM"):
            key = (line[21], int(line[22:26]))
            if key not in res:
                res[key] = {}
                order.append(key)
            res[key][line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz = np.array([res[k].get("CB", res[k].get("CA")) for k in order])
    lengths = {}
    for chain, _ in order:
        lengths[chain] = lengths.get(chain, 0) + 1
    return xyz, lengths


def score_one(npz, pdb) -> dict:
    pae = np.load(npz)["pae"].astype(float)
    xyz, L = cbeta_and_lengths(pdb)
    a, b, c = L["A"], L["B"], L["C"]
    assert pae.shape[0] == a + b + c, (pae.shape, L)
    fv, her2 = np.arange(0, a + b), np.arange(a + b, a + b + c)
    lis, clis, il = ilis(pae, fv, her2, xyz)
    return dict(ipsae=max(ipsae(pae, fv, her2), ipsae(pae, her2, fv)), lis=lis, clis=clis, ilis=il)


def prefilter(cfg: dict) -> pd.DataFrame:
    r = cfg["rank"]
    d = pd.read_csv(out(cfg, "fold_scores.csv")).merge(pd.read_csv(out(cfg, "contacts.csv")), on="id")
    design = ~d.id.str.startswith("ctrl_")
    ok = ((d.iptm >= r["min_iptm"]) & (d.if_recall >= r["min_interface_recall"]) & (d.if_h3_missing <= r["max_h3_pairs_missing"])
          & (d.in_recall >= r["min_h3_rest_recall"]) & (d.in_polar_recall >= r["min_h3_rest_polar_recall"]))
    return d[(design & ok) | ~design]


def run(cfg: dict) -> None:
    sel = prefilter(cfg)
    print(f"  {int((~sel.id.str.startswith('ctrl_')).sum())} designs pass the contact pre-filter (+ controls)")
    job = out(cfg, "job_id.txt").read_text().strip()
    pae_dir = out(cfg, "pae")
    JapanFold(cfg["fold"]["base_url"]).fetch_many(job, [f"{i}_pae.npz" for i in sel.id], pae_dir)
    rows = [dict(id=i, **score_one(pae_dir / f"{i}_pae.npz", out(cfg, "structures") / f"{i}.pdb")) for i in sel.id]
    pd.DataFrame(rows).to_csv(out(cfg, "scores.csv"), index=False)
    save_summary(cfg, n_pass_contacts=int((~sel.id.str.startswith("ctrl_")).sum()),
                 n_pass_contacts_h3=int(sel[~sel.id.str.startswith("ctrl_")].id.str.split("_").str[0].nunique()))
