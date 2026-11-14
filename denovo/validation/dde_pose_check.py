"""Independent re-check of the OpenDDE 'pose RMSD vs the designed pose' (S3's dde_own_vhh_rmsd_her2frame), to rule out a bug from the HER2
chain break / numbering. Predicted complex: chain A = VHH, chain B = HER2 domain IV (119 aa, HER2 489-607, contiguous). Design CIF: chain A = HER2
489-580 (resi 1-92), chain B = VHH, chain C = HER2 591-607 (resi 1-17). HER2 residues are matched BY HER2 NUMBER (A: 488+resi, C: 590+resi, predicted
B: 488+resi), Kabsch fit on the common HER2 CA atoms, the same transform applied to the predicted VHH CA atoms. Also: VHH-only fit (scRMSD), and where
the predicted VHH touches HER2 (fraction of contacted HER2 residues inside trastuzumab's epitope).
Run: baseline/.venv/bin/python denovo/validation/dde_pose_check.py   -> results/denovo/dde_pose_check.csv"""
import glob
import json
import os

import numpy as np
import pandas as pd

ROOT = os.getcwd()
EPI = {int(r["auth_seq"]) for r in json.load(open(f"{ROOT}/denovo/targets/epitope.json"))["residues"]}


def read(path):
    hdr, rows = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            rows.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol")}
    out = []
    for a in rows:
        if a[ix["type_symbol"]] == "H":
            continue
        out.append((a[ix["label_asym_id"]], int(a[ix["auth_seq_id"]]), a[ix["label_atom_id"]], np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])])))
    return out


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731
S3 = pd.read_csv(f"{ROOT}/results/denovo/score_orth.csv").set_index("design_id")
res = []
for cif in sorted(glob.glob(f"{ROOT}/results/denovo/s3_raw/dde/*/*__her2.cif"), key=os.path.getmtime):
    did = os.path.basename(cif)[: -len("__her2.cif")]
    dcif = f"{ROOT}/results/denovo/designs/{did}.cif"
    if not os.path.exists(dcif):
        continue
    P, D = read(cif), read(dcif)
    pv = {r: x for c, r, n, x in P if c == "A" and n == "CA"}
    ph = {488 + r: x for c, r, n, x in P if c == "B" and n == "CA"}
    dv = {r: x for c, r, n, x in D if c == "B" and n == "CA"}
    dh = {(488 + r if c == "A" else 590 + r): x for c, r, n, x in D if c in ("A", "C") and n == "CA"}
    common = sorted(set(ph) & set(dh))
    if len(common) < 80 or len(pv) != len(dv):
        continue
    R, t = kabsch(np.array([ph[k] for k in common]), np.array([dh[k] for k in common]))
    her_fit = rmsd(np.array([ph[k] for k in common]) @ R.T + t, np.array([dh[k] for k in common]))
    vp, vd = np.array([pv[k] for k in sorted(pv)]), np.array([dv[k] for k in sorted(dv)])
    pose = rmsd(vp @ R.T + t, vd)
    R2, t2 = kabsch(vp, vd)
    fold = rmsd(vp @ R2.T + t2, vd)
    # where does the predicted VHH touch HER2? (heavy atoms within 4.5 A) -> fraction of contacted HER2 residues inside the trastuzumab epitope
    va = np.array([x for c, r, n, x in P if c == "A"]); ha = [(r + 488, x) for c, r, n, x in P if c == "B"]
    touched = {r for r, x in ha if (np.linalg.norm(va - x, axis=1) <= 4.5).any()}
    dva = np.array([x for c, r, n, x in D if c == "B"]); dha = [((r + 488) if c == "A" else (r + 590), x) for c, r, n, x in D if c in ("A", "C")]
    dtouched = {r for r, x in dha if (np.linalg.norm(dva - x, axis=1) <= 4.5).any()}
    res.append({"design_id": did, "chk_her2_fit_rmsd": round(her_fit, 2), "chk_pose_rmsd": round(pose, 1), "chk_fold_rmsd": round(fold, 2),
                "chk_pred_epi_frac": round(len(touched & EPI) / max(1, len(touched)), 2), "chk_design_epi_frac": round(len(dtouched & EPI) / max(1, len(dtouched)), 2),
                "chk_pred_n_her2_contact": len(touched), "s3_pose": S3.dde_own_vhh_rmsd_her2frame.get(did, np.nan), "s3_fold": S3.dde_own_vhh_rmsd_fold.get(did, np.nan)})
R = pd.DataFrame(res).drop_duplicates("design_id", keep="last")
R.to_csv(f"{ROOT}/results/denovo/dde_pose_check.csv", index=False)
d = R.dropna(subset=["s3_pose"])
print(f"{len(R)} OpenDDE predictions re-checked")
print("HER2 fit RMSD (mapping sanity, should be ~1-3 A): median %.2f, p90 %.2f" % (R.chk_her2_fit_rmsd.median(), R.chk_her2_fit_rmsd.quantile(.9)))
print("pose RMSD  mine vs S3: median %.1f vs %.1f | mean abs difference %.2f A | Spearman %.3f" % (d.chk_pose_rmsd.median(), d.s3_pose.median(), (d.chk_pose_rmsd - d.s3_pose).abs().mean(), d.chk_pose_rmsd.corr(d.s3_pose, method="spearman")))
print("fold RMSD  mine vs S3: median %.2f vs %.2f | mean abs difference %.2f A" % (d.chk_fold_rmsd.median(), d.s3_fold.median(), (d.chk_fold_rmsd - d.s3_fold).abs().mean()))
print("pose <= 2 A: %d | <= 10 A: %d | <= 20 A: %d of %d" % ((R.chk_pose_rmsd <= 2).sum(), (R.chk_pose_rmsd <= 10).sum(), (R.chk_pose_rmsd <= 20).sum(), len(R)))
print("designed pose touches the epitope: median fraction %.2f | OpenDDE pose touches it: median fraction %.2f" % (R.chk_design_epi_frac.median(), R.chk_pred_epi_frac.median()))
