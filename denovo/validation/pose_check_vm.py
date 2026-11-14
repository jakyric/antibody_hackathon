"""Does the templated OpenFold3 (VM, HER2 template, no MSA) reproduce the DESIGNED binding pose better than the other scorers?

For every VM prediction (OpenFold3 `of3_*`, ESMFold-2 `esm2_*`; structures rsynced to results/denovo/vm_models/) and, for comparison, the OpenDDE
prediction of the same design (results/denovo/dde_pose_check.csv, same definitions):
  pose RMSD  = VHH CA RMSD to the designed VHH after superposing the predicted HER2 on the designed HER2 (CA atoms matched by HER2 number)
  scRMSD     = VHH CA RMSD after a VHH-only superposition (self-consistency of the fold)
  epi frac   = fraction of the HER2 residues that the predicted VHH touches (heavy atoms <= 4.5 A) that lie in trastuzumab's epitope
  pDockQ2    = same implementation as denovo/validation/dde_pdockq2.py (pLDDT from the CIF B-factors, PAE from the detailed scores)
VM chains: A = VHH, C = HER2 489-580 (resi 1-92), D = HER2 591-607 (resi 1-17). Design CIF: A = HER2 489-580, B = VHH, C = HER2 591-607.
Run: denovo/.venv/bin/python denovo/validation/pose_check_vm.py  -> results/denovo/pose_check_vm.csv (+ summary printed)
"""
import glob
import gzip
import json
import os
import sys

import cbor2
import numpy as np
import pandas as pd

ROOT = os.getcwd()
sys.path.insert(0, f"{ROOT}/denovo")
EPI = {int(r["auth_seq"]) for r in json.load(open(f"{ROOT}/denovo/targets/epitope.json"))["residues"]}
D = pd.read_csv(f"{ROOT}/results/denovo/designs_rfd4jf.csv").set_index("design_id")


def read(path, seq="auth_seq_id"):
    hdr, atoms = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            atoms.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", seq, "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "B_iso_or_equiv", "label_comp_id")}
    return [(a[ix["label_asym_id"]], int(a[ix[seq]]), a[ix["label_atom_id"]], np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])]),
             float(a[ix["B_iso_or_equiv"]])) for a in atoms if a[ix["type_symbol"]] != "H"]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731


def arr(tag):
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 40:
        shape, data = tag.value
        return arr(data).reshape(shape)
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 85:
        return np.frombuffer(tag.value, dtype="<f4")
    return np.array(tag)


def pdockq2(P, pae):
    """P: atoms of the predicted structure in token order VHH(A) then HER2 (C, D); pae: token x token."""
    toks = []  # (chain, resid) in order of first appearance; the CIF lists chains A, C, D in token order
    res = {}
    for c, r, n, x, b in P:
        k = (c, r)
        if k not in res:
            res[k] = {}
            toks.append(k)
        if n == "CA":
            res[k]["ca"], res[k]["plddt"] = x, b
        if n == "CB":
            res[k]["cb"] = x
    if pae.shape[0] != len(toks) or any("ca" not in res[k] for k in toks):
        return np.nan
    A = [i for i, k in enumerate(toks) if k[0] == "A"]
    B = [i for i, k in enumerate(toks) if k[0] != "A"]
    pos = lambda i: res[toks[i]].get("cb", res[toks[i]]["ca"])  # noqa: E731
    XA, XB = np.array([pos(i) for i in A]), np.array([pos(i) for i in B])
    d = np.linalg.norm(XA[:, None] - XB[None], axis=2)
    ii, jj = np.where(d <= 8.0)
    if not len(ii):
        return 0.005
    a, b = np.array(A)[ii], np.array(B)[jj]
    pa = (pae[a, b] + pae[b, a]) / 2.0
    plddt = np.mean([res[toks[i]]["plddt"] for i in set(a) | set(b)])
    x = plddt * np.mean(1.0 / (1.0 + (pa / 8.0) ** 2))
    return float(1.31 / (1 + np.exp(-0.075 * (x - 84.733))) + 0.005)


rows = []
for model, tag in (("of3", "vmof3"), ("esm2", "vmesm2")):
    for cif in glob.glob(f"{ROOT}/results/denovo/vm_models/{model}_*/*/seed_*/*_model.cif"):
        did = cif.split("/")[-3]
        if did not in D.index:  # controls
            continue
        dcif = f"{ROOT}/{D.loc[did, 'design_cif']}" if isinstance(D.loc[did, "design_cif"], str) and D.loc[did, "design_cif"] else f"{ROOT}/results/denovo/designs/{did}.cif"
        if not os.path.exists(dcif):
            continue
        P, Dd = read(cif, "label_seq_id"), read(dcif)  # VM CIFs keep HER2's original numbering in auth_seq_id, 1..N in label_seq_id
        pv = {r: x for c, r, n, x, b in P if c == "A" and n == "CA"}
        ph = {(488 + r if c == "C" else 590 + r): x for c, r, n, x, b in P if c in ("C", "D") and n == "CA"}
        dv = {r: x for c, r, n, x, b in Dd if c == "B" and n == "CA"}
        dh = {(488 + r if c == "A" else 590 + r): x for c, r, n, x, b in Dd if c in ("A", "C") and n == "CA"}
        common = sorted(set(ph) & set(dh))
        if len(common) < 80 or len(pv) != len(dv):
            continue
        R, t = kabsch(np.array([ph[k] for k in common]), np.array([dh[k] for k in common]))
        her_fit = rmsd(np.array([ph[k] for k in common]) @ R.T + t, np.array([dh[k] for k in common]))
        vp, vd = np.array([pv[k] for k in sorted(pv)]), np.array([dv[k] for k in sorted(dv)])
        pose = rmsd(vp @ R.T + t, vd)
        R2, t2 = kabsch(vp, vd)
        fold = rmsd(vp @ R2.T + t2, vd)
        va = np.array([x for c, r, n, x, b in P if c == "A"])
        ha = [((488 + r if c == "C" else 590 + r), x) for c, r, n, x, b in P if c in ("C", "D")]
        touched = {r for r, x in ha if (np.linalg.norm(va - x, axis=1) <= 4.5).any()}
        cb = cif.replace("_model.cif", "_detailed_scores.cbor.gz").replace("vm_models", "vm_raw")
        pq = np.nan
        if os.path.exists(cb):
            try:
                pq = pdockq2(P, arr(cbor2.loads(gzip.open(cb).read())["token_token_predicted_aligned_error"]))
            except Exception as e:  # noqa: BLE001
                print("pdockq2 skip", did, e)
        rows.append({"design_id": did, "scorer": tag, "her2_fit": round(her_fit, 2), "pose_rmsd": round(pose, 1), "scrmsd": round(fold, 2),
                     "epi_frac": round(len(touched & EPI) / max(1, len(touched)), 2), "n_her2_contact": len(touched), "pdockq2": round(pq, 3) if pq == pq else np.nan})
V = pd.DataFrame(rows)
dde = pd.read_csv(f"{ROOT}/results/denovo/dde_pose_check.csv")
x = pd.read_csv(f"{ROOT}/results/denovo/score_dde_extra.csv")[["design_id", "ddex_pdockq2"]]
dde = dde.merge(x, on="design_id", how="left")
O = pd.DataFrame({"design_id": dde.design_id, "scorer": "opendde", "her2_fit": dde.chk_her2_fit_rmsd, "pose_rmsd": dde.chk_pose_rmsd, "scrmsd": dde.chk_fold_rmsd,
                  "epi_frac": dde.chk_pred_epi_frac, "n_her2_contact": dde.chk_pred_n_her2_contact, "pdockq2": dde.ddex_pdockq2})
A = pd.concat([V, O], ignore_index=True).drop_duplicates(["design_id", "scorer"], keep="last")
A.to_csv(f"{ROOT}/results/denovo/pose_check_vm.csv", index=False)

# same designs for every scorer
common = set.intersection(*[set(A[A.scorer == s].design_id) for s in ("vmof3", "vmesm2", "opendde")])
print(f"designs with all three predictions: {len(common)} (of3 {A[A.scorer == 'vmof3'].design_id.nunique()}, esm2 {A[A.scorer == 'vmesm2'].design_id.nunique()}, opendde {A[A.scorer == 'opendde'].design_id.nunique()})")
for label, sub in (("all designs with the scorer", A), (f"the {len(common)} designs scored by all three", A[A.design_id.isin(common)])):
    print(f"\n== {label}")
    g = sub.groupby("scorer")
    out = pd.DataFrame({"n": g.size(), "HER2 fit A": g.her2_fit.median(), "pose RMSD med": g.pose_rmsd.median(), "pose<=5": g.pose_rmsd.apply(lambda s: int((s <= 5).sum())),
                        "pose<=10": g.pose_rmsd.apply(lambda s: int((s <= 10).sum())), "pose<=20": g.pose_rmsd.apply(lambda s: int((s <= 20).sum())),
                        "scRMSD med": g.scrmsd.median(), "epi frac med": g.epi_frac.median(), "pDockQ2 med": g.pdockq2.median(), "pDockQ2 p90": g.pdockq2.quantile(.9)})
    print(out.round(2).to_string())
W = A[A.design_id.isin(common)].pivot(index="design_id", columns="scorer", values="pose_rmsd")
print("\nSpearman of pose RMSD between scorers:\n", W.corr(method="spearman").round(2).to_string())
S = pd.read_csv(f"{ROOT}/results/denovo/submission_top20_pareto.csv").head(10)
print("\nCurrent top 10, pose RMSD A (designed pose vs prediction) / pDockQ2:")
T = A[A.design_id.isin(S.design_id)].pivot(index="design_id", columns="scorer", values=["pose_rmsd", "pdockq2"]).reindex(S.design_id)
print(T.round(2).to_string())
