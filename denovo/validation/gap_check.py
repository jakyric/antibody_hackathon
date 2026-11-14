"""Is the unresolved HER2 loop 581-590 (the chain break between design chains A and C) a problem for the designs?

Two questions, answered with numbers:
 1. Where is the gap relative to the epitope?  (distances from the gap termini to the 1N8Z epitope residues)
 2. Would a modelled loop clash with the designed VHHs?  We take the AlphaFold DB model of full-length HER2 (P04626), superpose its
    domain IV on the design's HER2 (chains A+C, CA atoms), keep the 581-590 loop atoms and count VHH heavy atoms close to them.
    Calibration: the same numbers for trastuzumab's Fv from 1N8Z (a real binder that sits next to the gap).
The AF loop is ONE plausible conformer of a flexible loop, so 'clash' = incompatible with that conformer, not proof of a real clash.
Run: denovo/.venv/bin/python denovo/validation/gap_check.py   -> results/denovo/gapcheck.csv
"""
import glob
import json
import os
import urllib.request

import numpy as np
import pandas as pd

ROOT = os.getcwd()
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-P04626-F1-model_v4.pdb"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def read_cif(path):
    hdr, rows = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            rows.append(l.split())
    df = pd.DataFrame(rows, columns=hdr[: len(rows[0])])
    df = df[df.type_symbol != "H"]
    return pd.DataFrame({"chain": df.label_asym_id.values, "res": df.auth_seq_id.astype(int).values, "name": df.label_atom_id.values,
                         "resn": df.label_comp_id.values, "xyz": list(df[["Cartn_x", "Cartn_y", "Cartn_z"]].astype(float).values)})


def xyz(df):
    return np.array(list(df.xyz)) if len(df) else np.zeros((0, 3))


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc  # x -> R x + t maps P onto Q


# ---- AlphaFold model, residue map by SEQUENCE (not by assumed numbering)
tgt = open(f"{ROOT}/denovo/targets/her2_domIV.a3m").read().split("\n")[1].strip()  # 119 aa domain IV (HER2 489-607)
af_path = f"{ROOT}/denovo/validation/AF-P04626-F1-model_v4.pdb"
if not os.path.exists(af_path):
    # the model version in the file name changes over time: ask the AFDB API for the current PDB URL
    meta = json.load(urllib.request.urlopen("https://alphafold.ebi.ac.uk/api/prediction/P04626", timeout=60))
    urllib.request.urlretrieve(meta[0]["pdbUrl"], af_path)
    print("AF model:", meta[0].get("pdbUrl"))
af = []
for l in open(af_path):
    if l.startswith("ATOM"):
        af.append((int(l[22:26]), l[12:16].strip(), l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66])))
af_seq = {}
for r, n, rn, x, b in af:
    af_seq[r] = THREE.get(rn, "X")
seq = "".join(af_seq[r] for r in sorted(af_seq))
s0 = sorted(af_seq)[seq.find(tgt)]
assert seq.find(tgt) >= 0, "domain IV sequence not found in the AFDB model"
off = s0 - 489  # HER2 number n (mature numbering of the crystal) -> AF residue n + off
af_by = {(r - off, n): (x, b) for r, n, rn, x, b in af}
print(f"AF model: domain IV found at AF residues {s0}-{s0 + 118} (offset {off}); mean pLDDT of loop 581-590: "
      f"{np.mean([af_by[(n, 'CA')][1] for n in range(581, 591)]):.0f}")

# ---- reference HER2 frame from one design; AF loop transformed into it
dfiles = sorted(glob.glob(f"{ROOT}/results/denovo/designs/*rfd4jf*.cif"))
ref = read_cif(dfiles[0])
her2_map = lambda df: df[(df.chain.isin(["A", "C"])) & (df.name == "CA")].assign(
    n=lambda d: np.where(d.chain == "A", d.res + 488, d.res + 590))  # chain A resi k = HER2 488+k, chain C resi k = HER2 590+k
rc = her2_map(ref)
P = np.array([af_by[(n, "CA")][0] for n in rc.n]); Q = xyz(rc)
R0, t0 = kabsch(P, Q)
print(f"AF domain IV onto design HER2: CA RMSD {np.sqrt((((P @ R0.T + t0) - Q) ** 2).sum(1).mean()):.2f} A over {len(rc)} CA")
loop = np.array([(af_by[(n, nm)][0] @ R0.T + t0) for n in range(581, 591) for (m, nm) in [(n, k[1]) for k in af_by if k[0] == n]])
t580 = af_by[(580, "CA")][0] @ R0.T + t0
t591 = af_by[(591, "CA")][0] @ R0.T + t0
print("loop atoms:", len(loop), "| termini CA-CA distance in the model: %.1f A" % np.linalg.norm(t580 - t591))

# ---- epitope vs gap (1N8Z)
ep = json.load(open(f"{ROOT}/denovo/targets/epitope.json"))["residues"]
z = read_cif(f"{ROOT}/denovo/targets/ref_1n8z.cif")
her = z[z.chain.isin(["C", "D"]) & (z.name == "CA")]
gcen = (t580 + t591) / 2
# put 1N8Z into the design frame via HER2 CA superposition
common = [(n) for n in rc.n if n in set(her.res)]
Pz = np.array([xyz(her[her.res == n])[0] for n in common]); Qz = np.array([rc[rc.n == n].iloc[0].xyz for n in common])
Rz, tz = kabsch(Pz, Qz)
dist_ep = {e["auth_seq"]: np.linalg.norm((xyz(her[her.res == e["auth_seq"]])[0] @ Rz.T + tz) - gcen) for e in ep if e["auth_seq"] in set(her.res)}
near = sorted(dist_ep.items(), key=lambda kv: kv[1])[:8]
print("epitope residues closest to the gap centre (A):", [(f"{n}", round(d, 1)) for n, d in near])
print("epitope residues within 10 A of the gap centre: %d of %d" % (sum(d < 10 for d in dist_ep.values()), len(dist_ep)))
fv = z[z.chain.isin(["A", "B"])]
fvx = xyz(fv) @ Rz.T + tz
dfv = np.linalg.norm(fvx[:, None] - loop[None], axis=2)
print("trastuzumab Fv vs AF loop: min dist %.1f A, atoms <3.0 A: %d, <4.0 A: %d  (a real binder next to the gap)" %
      (dfv.min(), (dfv < 3.0).any(1).sum(), (dfv < 4.0).any(1).sum()))

# ---- every design
rows = []
for f in dfiles:
    did = os.path.basename(f)[:-4]
    d = read_cif(f)
    h = her2_map(d)
    R, t = kabsch(xyz(h), np.array([rc[rc.n == n].iloc[0].xyz for n in h.n]))  # design HER2 -> reference frame
    v = xyz(d[d.chain == "B"]) @ R.T + t
    dd = np.linalg.norm(v[:, None] - loop[None], axis=2)
    rows.append({"design_id": did, "gap_min_dist": round(dd.min(), 2), "gap_atoms_lt3": int((dd < 3.0).any(1).sum()),
                 "gap_atoms_lt4": int((dd < 4.0).any(1).sum()), "gap_clash_lt2p5": int((dd < 2.5).sum()),
                 "gap_to_t580": round(np.linalg.norm(v - t580, axis=1).min(), 2), "gap_to_t591": round(np.linalg.norm(v - t591, axis=1).min(), 2)})
G = pd.DataFrame(rows)
G.to_csv(f"{ROOT}/results/denovo/gapcheck.csv", index=False)
print(f"\n{len(G)} designs: min distance VHH to AF loop  median {G.gap_min_dist.median():.1f} A (p10 {G.gap_min_dist.quantile(.1):.1f}, min {G.gap_min_dist.min():.1f})")
print("designs with a heavy-atom pair < 2.5 A to the AF loop:", int((G.gap_clash_lt2p5 > 0).sum()), "| with atoms < 3 A:", int((G.gap_atoms_lt3 > 0).sum()),
      "| atoms < 4 A:", int((G.gap_atoms_lt4 > 0).sum()))

# ---- stricter, loop-flexibility-aware view: atoms of the loop residues next to the anchors cannot move away
anchor = np.array([af_by[(n, nm)][0] @ R0.T + t0 for n in (581, 582, 589, 590) for nm in {k[1] for k in af_by if k[0] == n}])
anch_pts = np.vstack([anchor, t580[None], t591[None]])
dfa = np.linalg.norm(fvx[:, None] - anch_pts[None], axis=2)
print("\nANCHOR view (loop residues 581,582,589,590 + CA580/CA591):")
print("trastuzumab Fv: min distance to the anchor atoms %.1f A; atoms < 3.0 A: %d, < 3.5 A: %d" % (dfa.min(), (dfa < 3.0).any(1).sum(), (dfa < 3.5).any(1).sum()))
rows2 = []
for f in dfiles:
    d = read_cif(f); h = her2_map(d)
    R, t = kabsch(xyz(h), np.array([rc[rc.n == n].iloc[0].xyz for n in h.n]))
    v = xyz(d[d.chain == "B"]) @ R.T + t
    da = np.linalg.norm(v[:, None] - anch_pts[None], axis=2)
    rows2.append({"design_id": os.path.basename(f)[:-4], "gap_anchor_min": round(da.min(), 2), "gap_anchor_atoms_lt3": int((da < 3.0).any(1).sum())})
G2 = G.merge(pd.DataFrame(rows2), on="design_id")
G2.to_csv(f"{ROOT}/results/denovo/gapcheck.csv", index=False)
print("designs: min distance VHH to anchor atoms median %.1f A (p10 %.1f);  >=1 atom < 3.0 A: %d of %d (%.0f%%);  >=3 atoms < 3.0 A: %d" % (
    G2.gap_anchor_min.median(), G2.gap_anchor_min.quantile(.1), (G2.gap_anchor_atoms_lt3 > 0).sum(), len(G2), 100 * (G2.gap_anchor_atoms_lt3 > 0).mean(), (G2.gap_anchor_atoms_lt3 >= 3).sum()))
