"""Compare OpenDDE-abag on the cropped HER2 (results/denovo/crop/dde_crop.csv) with the full-domain OpenDDE runs on the same designs, and
print the decoy-based gate thresholds on the crop.   baseline/.venv/bin/python denovo/s3_score/crop_compare.py"""
import numpy as np
import pandas as pd

R = "results/denovo/"
C = pd.read_csv(R + "crop/dde_crop.csv")
S = pd.read_csv(R + "crop/set.csv").set_index("design_id")
F = pd.read_csv(R + "dde_pose_check.csv").set_index("design_id")
O = pd.read_csv(R + "score_orth.csv").set_index("design_id")
X = pd.read_csv(R + "score_dde_extra.csv").set_index("design_id")
C["why"] = C.design_id.map(S.why)
des = C[~C.design_id.str.startswith("ctl-")].copy()
des["full_pose"] = des.design_id.map(F.chk_pose_rmsd)
des["full_scrmsd"] = des.design_id.map(F.chk_fold_rmsd)
des["full_iptm"] = des.design_id.map(O.dde_iptm)
des["full_ipsae"] = des.design_id.map(O.dde_ipsae)
des["full_pdockq2"] = des.design_id.map(X.ddex_pdockq2)
des["full_epi"] = des.design_id.map(F.chk_pred_epi_frac)
d = des.dropna(subset=["full_pose", "dde_c_pose_rmsd"])
print(f"crop rows {len(C)} (designs {len(des)}, controls {len(C) - len(des)}); designs with both full and crop pose: {len(d)}")
t = pd.DataFrame({
    "median": [d.full_pose.median(), d.dde_c_pose_rmsd.median(), d.full_scrmsd.median(), d.dde_c_scrmsd.median(), d.full_iptm.median(), d.dde_c_iptm.median(),
               d.full_pdockq2.median(), d.dde_c_pdockq2.median(), d.full_epi.median(), d.dde_c_epi_frac.median()]},
    index=["pose RMSD full", "pose RMSD crop", "scRMSD full", "scRMSD crop", "iPTM full", "iPTM crop", "pDockQ2 full", "pDockQ2 crop", "epi frac full", "epi frac crop"])
print(t.round(3).to_string())
for lab, col in (("full", "full_pose"), ("crop", "dde_c_pose_rmsd")):
    print(f"{lab}: pose <= 5 A: {(d[col] <= 5).sum()} | <= 10 A: {(d[col] <= 10).sum()} | <= 20 A: {(d[col] <= 20).sum()} of {len(d)}")
print("Spearman full vs crop: pose %.2f | iPTM %.2f | pDockQ2 %.2f" % (d.full_pose.corr(d.dde_c_pose_rmsd, method="spearman"), d.full_iptm.corr(d.dde_c_iptm, method="spearman"),
                                                                      d.full_pdockq2.corr(d.dde_c_pdockq2, method="spearman")))
print("\ndesigns with the pose recovered on the crop (pose <= 10 A):")
print(des[des.dde_c_pose_rmsd <= 10][["design_id", "dde_c_pose_rmsd", "full_pose", "dde_c_scrmsd", "dde_c_iptm", "dde_c_ipsae", "dde_c_pdockq2", "dde_c_epi_frac"]].round(2).to_string(index=False))
print("\nclosest 8 poses on the crop:")
print(des.nsmallest(8, "dde_c_pose_rmsd")[["design_id", "dde_c_pose_rmsd", "full_pose", "dde_c_scrmsd", "dde_c_iptm", "dde_c_ipsae", "dde_c_pdockq2", "dde_c_epi_frac"]].round(2).to_string(index=False))
ctl = C[C.design_id.str.startswith("ctl-")]
dec = ctl[ctl.design_id.str.startswith("ctl-fwrand")]
print(f"\ncontrols on the crop: {len(ctl)} ({len(dec)} framework-only decoys)")
for lab, sub in (("decoys ctl-fwrand", dec), ("ctl-bcii10-wt", ctl[ctl.design_id == "ctl-bcii10-wt"]), ("ctl-kri-*", ctl[ctl.design_id.str.startswith("ctl-kri")])):
    if len(sub):
        print(f"  {lab:18s} n={len(sub):2d} | iPTM median {sub.dde_c_iptm.median():.2f} max {sub.dde_c_iptm.max():.2f} | ipSAE median {sub.dde_c_ipsae.median():.2f} max {sub.dde_c_ipsae.max():.2f} | "
              f"ipSAE chn15 median {sub.dde_c_ipsae_chn15.median():.2f} | pDockQ2 median {sub.dde_c_pdockq2.median():.2f} | epi frac median {sub.dde_c_epi_frac.median():.2f}")
if len(dec):
    q = lambda c: float(np.nanpercentile(dec[c], 95))  # noqa: E731
    print(f"\nNEW gate thresholds from the decoy 95th percentile on the crop: iPTM {q('dde_c_iptm'):.3f}, ipSAE {q('dde_c_ipsae'):.3f}, ipSAE chn15 {q('dde_c_ipsae_chn15'):.3f}, pDockQ2 {q('dde_c_pdockq2'):.3f}")
    gate = des[(des.dde_c_iptm >= q("dde_c_iptm")) & (des.dde_c_ipsae >= q("dde_c_ipsae"))]
    print(f"designs above both decoy thresholds on the crop: {len(gate)} of {len(des)}; gate_strict designs among them: {(gate.why == 'gate_strict').sum()} of {(des.why == 'gate_strict').sum()}")
