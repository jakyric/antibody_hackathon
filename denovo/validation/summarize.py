"""VAL: summary numbers for REPORT.md (reads results/denovo/validity.csv + designs; prints a compact text report)."""
import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import validate as V  # noqa: E402

OUT = V.OUT
D = pd.read_csv(OUT / "designs_rfd4jf.csv")
Vd = pd.read_csv(OUT / "validity.csv")
M = D.merge(Vd, on="design_id", suffixes=("", "_v"))
M["pose"] = M.params_json.map(lambda s: json.loads(s)["pose"])
M["caa"] = M.params_json.map(lambda s: json.loads(s).get("caa_fixed", 0))
n = len(M)
print(f"## N designs in designs_rfd4jf.csv: {n}  (bcii10 {int((M.framework=='bcii10').sum())}, capla {int((M.framework=='capla').sum())}; p0 {int((M.pose=='p0').sum())}, p1 {int((M.pose=='p1').sum())}; caa_fixed {int(M.caa.sum())})")
print("val_ok:", int(M.val_ok.sum()))
flags = [c for c in Vd.columns if c.startswith("val_") and c != "val_rama_cdr3_outlier_frac"]
print(M.groupby("framework")[flags].sum().T.to_string())

print("\n## 2 fixed parts (median / p95 / max)")
for c in ("kept_group_rmsd", "her2_ca_rmsd", "fw_pose_rmsd_on_her2"):
    g = M[M.val_motif_bad == 0][c]
    print(f"{c:24s} valid designs: median {g.median():.2f} p95 {g.quantile(.95):.2f} max {g.max():.2f}")
print("by framework (valid only) kept_group median:", M[M.val_motif_bad == 0].groupby("framework").kept_group_rmsd.median().round(2).to_dict())
print("wave1 bcii10 designs with changed A106 (AA after Cys104):", int(M[(M.framework == "bcii10") & (M.caa == 0)].wave1_AA_changed.sum()), "of", int(((M.framework == "bcii10") & (M.caa == 0)).sum()))

print("\n## 3 covalent geometry")
for c in ("cn_sd_A", "cn_sd_B", "cn_sd_C", "cn_min_B", "cn_max_B", "cn_outlier_frac_A", "cn_outlier_frac_B", "caca_min_A", "caca_max_A", "caca_min_B", "caca_max_B"):
    print(f"{c:20s} median {M[c].median():.3f}  min {M[c].min():.2f}  max {M[c].max():.2f}")
print("breaks (CA-CA>4.3):", {k: int(M[f'breaks_{k}'].sum()) for k in "ABC"}, "| cis non-Pro:", int((M.cis_nonpro > 0).sum()), "| D-amino:", int((M.d_amino > 0).sum()),
      "| incomplete residues:", int((M.incomplete_res > 0).sum()), "| H atoms:", int(M.H_atoms.sum()), "| zero/nan xyz:", int((M.zero_or_nan_xyz > 0).sum()))
tot = {r: (int(M[f"rama_out_{r}"].sum()), int(M[f"rama_n_{r}"].sum())) for r in ("fw", "cdr1", "cdr2", "cdr3")}
print("Ramachandran outliers (coarse):", {k: f"{a}/{b} = {a/max(b,1):.1%}" for k, (a, b) in tot.items()})
print("geometry-bad designs (val_geometry_bad):", int(M.val_geometry_bad.sum()), M[M.val_geometry_bad == 1].design_id.tolist()[:8])
bad = M[M.val_geometry_bad == 1]
if len(bad):
    print(bad[["design_id", "bond_angle_bad", "cis_nonpro", "d_amino", "breaks_A", "breaks_B", "breaks_C"]].to_string(index=False))

print("\n## 4 disulfide")
print("SG-SG median by framework:", M.groupby("framework").ss_sg_sg.median().round(2).to_dict(), "| stretched 2.4-3.0:", int(M.ss_stretched.sum()), "| not formed >3.0:", int(M.val_ss_not_formed.sum()))
print("HER2 Cys with a partner SG<2.4A:", int(M.her2_ss_formed.median()), "of", int(M.her2_cys.median()), "(median)")
print("unexpected Cys count (!=2):", M[M.val_cys_unexpected == 1][["design_id", "cdr3"]].values.tolist())

print("\n## 5 clashes and interface")
print("intra-VHH clash (<2.2A, non-adjacent) designs:", int((M.clash_intra_vhh > 0).sum()), "| VHH-HER2 <2.2A:", int((M.clash_vhh_her2_lt2p2 > 0).sum()), "| <2.5A:", int((M.clash_vhh_her2_lt2p5 > 0).sum()),
      "| val_clash_bad:", int(M.val_clash_bad.sum()), "of", n)
print(M.groupby(["framework", "caa"]).agg(n=("design_id", "size"), intra=("clash_intra_vhh", lambda x: int((x > 0).sum())), her2=("clash_vhh_her2_lt2p2", lambda x: int((x > 0).sum()))).to_string())
print("min VHH-HER2 distance: median", round(M.min_dist_vhh_her2.median(), 2), " p5", round(M.min_dist_vhh_her2.quantile(.05), 2))
print("designs with NO contact:", int(M.val_no_contact.sum()))
for c in ("vhh_res_in_contact", "contact_res_cdr1", "contact_res_cdr2", "contact_res_cdr3", "contact_res_fw", "her2_res_contact", "epi_res_contact", "epi_precision", "hot_touched", "hot_total", "cdr3_cen_to_epi_cen"):
    print(f"{c:22s} median {M[c].median():.2f}  p10 {M[c].quantile(.1):.2f}  p90 {M[c].quantile(.9):.2f}")
print("CDR3 contacts epitope residues in all designs? val_cdr3_no_epitope_contact sum:", int(M.val_cdr3_no_epitope_contact.sum()))
print("by pose: epitope residues contacted:", M.groupby("pose").epi_res_contact.median().to_dict(), "| hotspots touched/total:", (M.hot_touched.sum() / M.hot_total.sum()).round(2))

print("\n## 6 competition geometry (VHH vs trastuzumab Fv after superposing HER2 on 1N8Z)")
print(M.groupby("pose")[["fv_clash_lt2p5", "fv_atoms_lt4", "her2_ca_rmsd_to_1n8z"]].median().round(1).to_string())
print("designs with zero steric overlap with the Fv:", int((M.fv_clash_lt2p5 == 0).sum()), "| overlap >100 atoms:", int((M.fv_clash_lt2p5 > 100).sum()), "of", n)

print("\n## 9 pose diversity (framework relative to the p0 input, in the 1N8Z HER2 frame)")
print(M[M.val_motif_bad == 0].groupby("pose")[["pose_rot_deg_vs_p0", "pose_shift_vs_p0", "pose_fw_rmsd_to_p0", "pose_cen_to_trast_vh"]].agg(["median", "std", "max"]).round(2).to_string())

print("\n## 7 sequence sanity")
print("unique vhh_seq:", M.vhh_seq.nunique(), "of", n, "| unique CDR3:", M.cdr3.nunique(), "| unique (cdr1,cdr2,cdr3):", len(set(zip(M.cdr1, M.cdr2, M.cdr3))))
print("seeds unique:", M.seed.nunique(), "| duplicated vhh_seq rows:", int(M.vhh_seq.duplicated().sum()))
try:
    X = pd.read_csv(V.ROOT / "denovo" / "s1_design" / "duplicates_removed.csv")
    print("S1a duplicates_removed.csv rows:", len(X), "| unique seqs among them:", X.vhh_seq.nunique(), "| of those also kept:", int(X.vhh_seq.isin(M.vhh_seq).sum()))
except Exception as e:  # noqa: BLE001
    print("no duplicates file", e)
# near duplicates: CDR3 edit distance <= 2 among kept designs of the same framework (sample-based O(n^2) on short strings)
def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[-1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]
c3 = list(M.cdr3.unique())
near = 0
for i in range(len(c3)):
    for j in range(i + 1, len(c3)):
        if abs(len(c3[i]) - len(c3[j])) <= 2 and lev(c3[i], c3[j]) <= 2:
            near += 1
print("CDR3 pairs with edit distance <= 2:", near, "of", len(c3) * (len(c3) - 1) // 2)
# composition vs natural CDRH3
des = Counter("".join(M.cdr3.dropna().map(lambda s: s[2:-0] if False else s)))
tot = sum(des.values())
nat = None
try:
    P = pd.read_parquet(V.ROOT / "denovo" / "s2_rules" / "cache" / "abguard_cdrh3.parquet")
    col = [c for c in P.columns if "cdr" in c.lower() or "seq" in c.lower()][0]
    s = P[col].dropna().astype(str)
    s = s[s.str.len().between(8, 25)].sample(min(len(s), 200000), random_state=0)
    nat = Counter("".join(s))
    print("natural CDRH3 reference:", col, len(s))
except Exception as e:  # noqa: BLE001
    print("natural reference unavailable:", str(e)[:80])
if nat:
    nt = sum(nat.values())
    rows = []
    for a in "ACDEFGHIKLMNPQRSTVWY":
        rows.append((a, des[a] / tot, nat[a] / nt))
    R = pd.DataFrame(rows, columns=["aa", "designed", "natural"])
    R["ratio"] = R.designed / R.natural
    print(R.sort_values("ratio", ascending=False).round(3).to_string(index=False))
ends = Counter(s[-3:] for s in M.cdr3)
print("last-3 motifs (designed):", ends.most_common(8))
starts = Counter(s[2:5] for s in M.cdr3)
print("first-3 after AA:", starts.most_common(5))
print("homopolymer>=4 in a CDR:", int(M.homopoly4.sum()), "| dipeptide repeat x3:", int(M.dipep_repeat.sum()), "| low complexity flag:", int(M.val_low_complexity.sum()), "| CDR N-glyc motif:", int((M.nglyc_cdr > 0).sum()))
print("CDR lengths (designed block): cdr1", M.cdr1_len.describe()[["min", "50%", "max"]].to_dict(), "cdr2", M.cdr2_len.describe()[["min", "50%", "max"]].to_dict(), "cdr3 (block)", M.cdr3_len.describe()[["min", "50%", "max"]].to_dict())

print("\n## 8 self-consistency (NanoBodyBuilder2)")
try:
    S = pd.read_csv(OUT / "validity_selfconsist.csv")
    cols = [c for c in S.columns if c.startswith("nb2_") and ("rmsd" in c)]
    ok = S[S.nb2_fw_ca_rmsd < 3]
    print("n:", len(S), " framework-consistent:", len(ok), "(the rest are the invalid motif designs:", S[S.nb2_fw_ca_rmsd >= 3].design_id.tolist(), ")")
    print(ok[["nb2_fw_ca_rmsd", "nb2_cdr1_bb_rmsd_on_fw", "nb2_cdr2_bb_rmsd_on_fw", "nb2_cdr3_bb_rmsd_on_fw"]].describe().loc[["min", "50%", "max"]].round(2).to_string())
    oc = [c for c in S.columns if c.startswith("of3_")]
    print(ok[oc].describe().loc[["count", "min", "50%", "max"]].round(1).to_string())
except Exception as e:  # noqa: BLE001
    print("selfconsist not available", e)
