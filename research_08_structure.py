#!/usr/bin/env python
"""Part 8: structural evidence. Which CDRH3 positions actually touch HER2
in 1N8Z, and does that agree with the SPR position-criticality ranking?

Parses the hackathon's own her2_fv_1n8z.cif (chains: A=VL, B=VH,
C/D=HER2 domain IV). CDRH3 = B99-108.
"""
import numpy as np
from collections import defaultdict

CIF = "/home/hackathon/rfd4/Data/12-her2/her2_fv_1n8z.cif"
WT = "WGGDGFYAMD"

atoms = []
with open(CIF) as fh:
    in_loop = False
    cols = []
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("_atom_site."):
            in_loop = True
            cols.append(line.strip().split(".")[1])
            continue
        if in_loop:
            if line.startswith("#") or line.startswith("loop_") or not line.strip():
                if atoms:
                    break
                continue
            p = line.split()
            if len(p) < len(cols):
                continue
            d = dict(zip(cols, p))
            atoms.append(d)

print(f"parsed {len(atoms)} atom rows; columns: {cols[:12]} ...")
key_chain = "label_asym_id" if "label_asym_id" in cols else "auth_asym_id"
key_seq = "label_seq_id" if "label_seq_id" in cols else "auth_seq_id"

def fnum(x):
    try: return float(x)
    except: return np.nan

recs = []
for a in atoms:
    if a.get("group_PDB") not in (None, "ATOM"):
        continue
    xyz = (fnum(a["Cartn_x"]), fnum(a["Cartn_y"]), fnum(a["Cartn_z"]))
    if any(np.isnan(v) for v in xyz):
        continue
    try:
        sid = int(a[key_seq])
    except Exception:
        continue
    recs.append((a[key_chain], sid, a["label_comp_id"], a["label_atom_id"], xyz))

ch = defaultdict(set)
for c, s, comp, at, xyz in recs:
    ch[c].add(s)
print("chains:", {c: f"{min(v)}-{max(v)} ({len(v)} res)" for c, v in sorted(ch.items())})

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

# verify the CDRH3 sequence at B99-108
seqB = {}
for c, s, comp, at, xyz in recs:
    if c == "B":
        seqB[s] = AA3.get(comp, "X")
got = "".join(seqB.get(i, "-") for i in range(99, 109))
print(f"\nchain B 99-108 = {got}   (expected {WT})  match: {got == WT}")

# antigen atoms
ag = np.array([xyz for c, s, comp, at, xyz in recs if c in ("C", "D")])
ag_id = [(c, s, comp, at) for c, s, comp, at, xyz in recs if c in ("C", "D")]
print(f"HER2 atoms: {len(ag)}")

print("\n== per-CDRH3-position contacts with HER2 (chains C+D) ==")
print(f"{'pos':>4}{'res':>5}{'n_atoms':>9}{'min_dist':>10}{'contacts<4A':>13}"
      f"{'contacts<5A':>13}  nearest HER2 residues (<5 A)")
rowsout = []
for i, p in enumerate(range(99, 109), 1):
    at_i = np.array([xyz for c, s, comp, at, xyz in recs if c == "B" and s == p])
    if len(at_i) == 0:
        print(f"{i:>4}{WT[i-1]:>5}  (no coordinates)")
        continue
    d = np.linalg.norm(at_i[:, None, :] - ag[None, :, :], axis=2)
    c4 = int((d < 4.0).sum()); c5 = int((d < 5.0).sum())
    near = sorted({f"{ag_id[j][0]}{ag_id[j][1]}{AA3.get(ag_id[j][2],'X')}"
                   for j in np.where(d.min(axis=0) < 5.0)[0]})
    rowsout.append((i, WT[i-1], len(at_i), d.min(), c4, c5))
    print(f"{i:>4}{WT[i-1]:>5}{len(at_i):>9}{d.min():>10.2f}{c4:>13}{c5:>13}  "
          f"{' '.join(near[:8])}")

# side-chain-only contacts (what a substitution actually changes)
print("\n== side-chain-only contacts (backbone N,CA,C,O excluded) ==")
bb = {"N", "CA", "C", "O"}
print(f"{'pos':>4}{'res':>5}{'sc_atoms':>10}{'sc_min_dist':>13}{'sc_<4A':>8}{'sc_<5A':>8}")
sc_stats = []
for i, p in enumerate(range(99, 109), 1):
    at_i = np.array([xyz for c, s, comp, at, xyz in recs
                     if c == "B" and s == p and at not in bb])
    if len(at_i) == 0:
        print(f"{i:>4}{WT[i-1]:>5}{0:>10}{'-':>13}{0:>8}{0:>8}  (glycine)")
        sc_stats.append((i, WT[i-1], 0, np.inf, 0, 0))
        continue
    d = np.linalg.norm(at_i[:, None, :] - ag[None, :, :], axis=2)
    sc_stats.append((i, WT[i-1], len(at_i), d.min(), int((d<4).sum()), int((d<5).sum())))
    print(f"{i:>4}{WT[i-1]:>5}{len(at_i):>10}{d.min():>13.2f}"
          f"{int((d<4).sum()):>8}{int((d<5).sum()):>8}")

# compare with the SPR criticality ranking measured in part 2
spr_delta = {1:0.91, 2:0.29, 3:0.21, 4:0.53, 5:0.49, 6:0.99, 7:1.00,
             8:0.42, 9:1.38, 10:0.37}
print("\n== structure vs measured criticality ==")
print(f"{'pos':>4}{'WT':>4}{'sc_contacts<5A':>16}{'SPR dlog10KD if mutated':>26}")
from scipy.stats import spearmanr
xs, ys = [], []
for i, aa, n, mn, c4, c5 in sc_stats:
    print(f"{i:>4}{aa:>4}{c5:>16}{spr_delta[i]:>26.2f}")
    xs.append(c5); ys.append(spr_delta[i])
print(f"\n   Spearman(side-chain contacts, SPR criticality) = "
      f"{spearmanr(xs, ys)[0]:+.3f}")
print("   (positive = residues that touch HER2 are the ones you must not change)")
print("\n   the hackathon spec's declared hotspots: D593 NZ, C573 CD1, D591 CG2, C572 CD")
