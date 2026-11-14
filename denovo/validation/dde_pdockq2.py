"""pDockQ2-style interface quality for the OpenDDE-abag predictions (VHH chain A vs HER2 chain B), from the per-residue pLDDT (CIF B-factors)
and the PAE matrix that S3 saved (results/denovo/s3_raw/dde/*/<design>__her2.cif and __her2_pae.npz).

    pDockQ2 = L / (1 + exp(-k (x - x0))) + b,   L=1.31, x0=84.733, k=0.075, b=0.005   (Zhu et al., Bioinformatics 2023)
    x = <pLDDT of the interface residues> * <1 / (1 + (PAE_ij / d0)^2)> over the interface residue pairs, d0 = 8 A
Interface pair = VHH residue i and HER2 residue j with CB (CA for Gly) within 8 A; PAE_ij symmetrised as (PAE[i,j] + PAE[j,i]) / 2.
This is my implementation of the published sigmoid; the constants were checked against documentation of the ipSAE module (not the
original paper code), the interface definition and d0 are the part to double-check. Output: results/denovo/score_dde_extra.csv (ddex_*).
"""
import glob
import os

import numpy as np
import pandas as pd

ROOT = os.getcwd()
rows = {}
for cif in sorted(glob.glob(f"{ROOT}/results/denovo/s3_raw/dde/*/*__her2.cif"), key=os.path.getmtime):  # later jobs overwrite earlier ones
    did = os.path.basename(cif)[: -len("__her2.cif")]
    npz = cif.replace("__her2.cif", "__her2_pae.npz")
    if not os.path.exists(npz):
        continue
    hdr, atoms = [], []
    for l in open(cif):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            atoms.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", "B_iso_or_equiv", "Cartn_x", "Cartn_y", "Cartn_z", "label_comp_id", "auth_seq_id")}
    res = {}
    for a in atoms:
        ch, sid, name, rn = a[ix["label_asym_id"]], int(a[ix["auth_seq_id"]]), a[ix["label_atom_id"]], a[ix["label_comp_id"]]
        r = res.setdefault((ch, sid), {"rn": rn})
        xyz = np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])])
        if name == "CA":
            r["ca"], r["plddt"] = xyz, float(a[ix["B_iso_or_equiv"]])
        if name == "CB":
            r["cb"] = xyz
    keys = sorted(res, key=lambda k: (k[0], k[1]))
    A = [k for k in keys if k[0] == "A"]
    B = [k for k in keys if k[0] == "B"]
    pae = np.load(npz)["pae"]
    if pae.shape[0] != len(A) + len(B):
        continue
    pos = lambda k: res[k].get("cb", res[k]["ca"])  # noqa: E731
    XA, XB = np.array([pos(k) for k in A]), np.array([pos(k) for k in B])
    d = np.linalg.norm(XA[:, None] - XB[None], axis=2)
    ii, jj = np.where(d <= 8.0)
    if len(ii) == 0:
        rows[did] = {"design_id": did, "ddex_n_if_pairs": 0, "ddex_pdockq2": 0.005}
        continue
    pa = (pae[ii, len(A) + jj] + pae[len(A) + jj, ii]) / 2.0
    int_res = {A[i] for i in ii} | {B[j] for j in jj}
    plddt = np.mean([res[k]["plddt"] for k in int_res])
    x = plddt * np.mean(1.0 / (1.0 + (pa / 8.0) ** 2))
    rows[did] = {"design_id": did, "ddex_n_if_pairs": int(len(ii)), "ddex_if_plddt": round(float(plddt), 1), "ddex_if_pae": round(float(pa.mean()), 2),
                 "ddex_pdockq2": round(float(1.31 / (1 + np.exp(-0.075 * (x - 84.733))) + 0.005), 3)}
out = pd.DataFrame(list(rows.values()))
out.to_csv(f"{ROOT}/results/denovo/score_dde_extra.csv", index=False)
print(len(out), "OpenDDE predictions scored; pDockQ2 median %.3f, p90 %.3f, max %.3f" % (out.ddex_pdockq2.median(), out.ddex_pdockq2.quantile(.9), out.ddex_pdockq2.max()))
