"""PyMOL session: DESIGNED structure (RFD4 output) vs RE-PREDICTED OpenDDE-abag structure for the current Pareto top 10.

Run from the repo root with the PyMOL conda env:
    micromamba run -p $PYMOL_ENV pymol -cq denovo/validation/make_pymol_dde_compare.py

Per design, two objects in ONE frame (the designed HER2, held rigid by RFD4):
  <rk>_design   designed complex: HER2 grey, VHH light blue, CDR1 blue, CDR2 green, CDR3 red
  <rk>_dde      OpenDDE-abag prediction: VHH only (violet, CDR3 magenta), superposed on the designed HER2 by HER2 residue number
                (Kabsch on common CA atoms, the same transform as denovo/validation/dde_pose_check.py); its HER2 is hidden
                (`<rk>_dde_her2` is loaded but disabled, for checking how well the two HER2 copies agree)
Reference: trastuzumab Fv (1N8Z), transparent wheat/orange; trastuzumab epitope residues as orange sticks on HER2.
Scenes: one per design (design only / prediction only / both); the RMSDs are printed and written to
denovo/validation/pymol/dde_compare.csv. Numbering: design CIF chain A = HER2 489-580 (resi 1-92), chain B = VHH, chain C = HER2 591-607 (resi 1-17);
OpenDDE CIF chain A = VHH, chain B = HER2 489-607 (resi 1-119).
"""
import glob
import json
import os

import csv

import numpy as np
from pymol import cmd

ROOT = os.getcwd()
OUT = f"{ROOT}/denovo/validation/pymol"
N = 10
rd = lambda p, key=None: list(csv.DictReader(open(f"{ROOT}/results/denovo/{p}")))  # noqa: E731
S = rd("submission_top20_pareto.csv")[:N]
D = {r["design_id"]: r for r in rd("designs_rfd4jf.csv")}
SO = {r["design_id"]: r for r in rd("score_orth.csv")}
num = lambda d, k: float(d[k]) if d.get(k) not in (None, "") else float("nan")  # noqa: E731
EPI = [int(r["auth_seq"]) for r in json.load(open(f"{ROOT}/denovo/targets/epitope.json"))["residues"]]


def read_ca(path):
    hdr, rows = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            rows.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z")}
    return [(a[ix["label_asym_id"]], int(a[ix["auth_seq_id"]]), np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])]))
            for a in rows if a[ix["label_atom_id"]] == "CA"]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc  # x -> R x + t maps P onto Q


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731


def her2(obj, n):
    return f"({obj} and chain A and resi {n - 488})" if n <= 580 else f"({obj} and chain C and resi {n - 590})"


cmd.reinitialize()
cmd.bg_color("white")
for k, v in (("ray_opaque_background", 1), ("antialias", 2), ("cartoon_side_chain_helper", 1), ("ray_shadow", 0), ("depth_cue", 0)):
    cmd.set(k, v)

cmd.load(f"{ROOT}/denovo/targets/ref_1n8z.cif", "tz_ref")
cmd.hide("everything", "tz_ref")
cmd.show("cartoon", "tz_ref and chain A+B")
cmd.color("wheat", "tz_ref and chain A")
cmd.color("orange", "tz_ref and chain B")
cmd.set("cartoon_transparency", 0.65, "tz_ref")

info, first = [], None
for i, r in enumerate(S, start=1):
    did = r["design_id"]
    rk = f"{i:02d}"
    dcif = f"{ROOT}/{D[did]['design_cif']}" if D[did]["design_cif"] else f"{ROOT}/results/denovo/designs/{did}.cif"
    pcifs = sorted(glob.glob(f"{ROOT}/results/denovo/s3_raw/dde/*/{did}__her2.cif"), key=os.path.getmtime)
    if not pcifs or not os.path.exists(dcif):
        print("MISSING", did, bool(pcifs), os.path.exists(dcif))
        continue
    pcif = pcifs[-1]
    # transform: OpenDDE HER2 -> designed HER2, by HER2 residue number
    P, Dd = read_ca(pcif), read_ca(dcif)
    ph = {488 + n: x for c, n, x in P if c == "B"}
    dh = {(488 + n if c == "A" else 590 + n): x for c, n, x in Dd if c in ("A", "C")}
    common = sorted(set(ph) & set(dh))
    Rm, t = kabsch(np.array([ph[k] for k in common]), np.array([dh[k] for k in common]))
    her_fit = rmsd(np.array([ph[k] for k in common]) @ Rm.T + t, np.array([dh[k] for k in common]))
    pv = np.array([x for c, n, x in sorted(P, key=lambda z: z[1]) if c == "A"])
    dv = np.array([x for c, n, x in sorted(Dd, key=lambda z: z[1]) if c == "B"])
    pose = rmsd(pv @ Rm.T + t, dv)
    R2, t2 = kabsch(pv, dv)
    fold = rmsd(pv @ R2.T + t2, dv)
    M = [*Rm[0], t[0], *Rm[1], t[1], *Rm[2], t[2], 0, 0, 0, 1]

    seq = D[did]["vhh_seq"]
    ob_d, ob_p, ob_h = f"{rk}_design", f"{rk}_dde", f"{rk}_dde_her2"
    cmd.load(dcif, ob_d)
    cmd.load(pcif, ob_p)
    cmd.create(ob_h, f"{ob_p} and chain B")
    cmd.remove(f"{ob_p} and chain B")
    cmd.transform_selection(ob_p, M, homogenous=0)
    cmd.transform_selection(ob_h, M, homogenous=0)
    for ob in (ob_d, ob_p, ob_h):
        cmd.hide("everything", ob)
    # designed complex
    cmd.show("cartoon", ob_d)
    cmd.color("gray70", f"{ob_d} and chain A+C")
    cmd.color("lightblue", f"{ob_d} and chain B")
    cmd.show("cartoon", ob_p)
    cmd.color("violet", ob_p)
    for nm, dcol, pcol in (("cdr1", "tv_blue", "purple"), ("cdr2", "tv_green", "purple"), ("cdr3", "tv_red", "magenta")):
        s = D[did][nm]
        k = seq.find(s)
        if k < 0:
            continue
        sel = f"resi {k + 1}-{k + len(s)}"
        cmd.color(dcol, f"{ob_d} and chain B and {sel}")
        cmd.color(pcol, f"{ob_p} and chain A and {sel}")
        if nm == "cdr3":
            cmd.show("sticks", f"{ob_d} and chain B and {sel} and not name H*")
            cmd.show("sticks", f"{ob_p} and chain A and {sel} and not name H*")
            cmd.color("magenta", f"{ob_p} and chain A and {sel} and elem C")
    cmd.set("stick_radius", 0.12, ob_d)
    cmd.set("stick_radius", 0.12, ob_p)
    cmd.select(f"{ob_d}_epi", " or ".join(her2(ob_d, n) for n in EPI))
    cmd.show("sticks", f"{ob_d}_epi and not name H*")
    cmd.color("orange", f"{ob_d}_epi and elem C")
    cmd.set("cartoon_transparency", 0.0, ob_d)
    cmd.disable(ob_h)
    iptm = num(SO.get(did, {}), "dde_iptm")
    info.append({"rank": i, "design_id": did, "framework": r["framework"], "pose_rmsd_A": round(pose, 1), "scrmsd_A": round(fold, 2),
                 "her2_fit_rmsd_A": round(her_fit, 2), "dde_iptm": round(float(iptm), 2), "dde_ipsae_chn15": round(num(SO.get(did, {}), "dde_ipsae_chn15"), 2),
                 "design_cif": os.path.relpath(dcif, ROOT), "dde_cif": os.path.relpath(pcif, ROOT)})
    print(f"{rk} {did}: pose RMSD {pose:.1f} A | scRMSD {fold:.2f} A | HER2 fit {her_fit:.2f} A | OpenDDE iPTM {iptm:.2f}")
    if first is None:
        first = rk
    cmd.align("tz_ref and chain C+D and name CA", f"{ob_d} and chain A+C and name CA") if i == 1 else None

with open(f"{OUT}/dde_compare.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(info[0]))
    w.writeheader()
    w.writerows(info)


def snap(rk, mode):
    cmd.disable("all")
    cmd.enable("tz_ref")
    if mode in ("design", "both"):
        cmd.enable(f"{rk}_design")
        cmd.enable(f"{rk}_design_epi")
    if mode in ("dde", "both"):
        cmd.enable(f"{rk}_dde")
        if mode == "dde":
            cmd.enable(f"{rk}_design")
            cmd.hide("sticks", f"{rk}_design and chain B")
            cmd.hide("cartoon", f"{rk}_design and chain B")
    cmd.orient(f"{rk}_design and chain B or {rk}_dde or {rk}_design_epi")
    cmd.zoom(f"{rk}_design and chain B or {rk}_dde or {rk}_design_epi", 6)


for row in info:
    rk = f"{row['rank']:02d}"
    snap(rk, "both")
    cmd.scene(f"{rk}_both_{row['design_id']}", "store")
    cmd.ray(1300, 950)
    cmd.png(f"{OUT}/dde_{rk}_{row['design_id']}.png", dpi=150)
    # restore the full designed VHH for the next scenes
    cmd.show("cartoon", f"{rk}_design and chain B")
    cmd.show("sticks", f"{rk}_design and chain B and resn CYS")

# first design: both objects shown, others disabled
snap(first, "both")
cmd.set("cartoon_transparency", 0.65, "tz_ref")
cmd.save(f"{OUT}/rfd4_top10_designed_vs_dde.pse")
print("saved", f"{OUT}/rfd4_top10_designed_vs_dde.pse", "| table", f"{OUT}/dde_compare.csv")
