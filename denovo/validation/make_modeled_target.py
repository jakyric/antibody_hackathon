"""Modelled HER2 domain IV target WITHOUT the chain break: AlphaFold DB model of P04626 (domain IV, residues 489-607 in the crystal's
mature numbering = AF 511-629), superposed onto the 1N8Z crystal HER2 (denovo/targets/her2_domIV.cif) so the VHH framework poses
that were anchored on trastuzumab's VH stay valid. One chain A, residues 1-119 (resi k = HER2 488+k), all heavy atoms, B-factor = pLDDT.
Run: denovo/.venv/bin/python denovo/validation/make_modeled_target.py   -> denovo/targets/her2_domIV_modeled.cif
"""
import json
import os

import numpy as np

ROOT = os.getcwd()
tgt = open(f"{ROOT}/denovo/targets/her2_domIV.a3m").read().split("\n")[1].strip()
af = {}
for l in open(f"{ROOT}/denovo/validation/AF-P04626-F1-model_v4.pdb"):
    if l.startswith("ATOM"):
        af.setdefault(int(l[22:26]), []).append((l[12:16].strip(), l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), float(l[60:66])))
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L",
         "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
seq = "".join(THREE[af[r][0][1]] for r in sorted(af))
s0 = sorted(af)[seq.find(tgt)]
assert seq.find(tgt) >= 0
# crystal reference (her2_domIV.cif: chain A, label_seq contiguous over modelled residues, auth_seq = HER2 number)
cr = {}
for l in open(f"{ROOT}/denovo/targets/her2_domIV.cif"):
    t = l.split()
    if l.startswith("ATOM") and t[3] == "CA":
        cr[int(t[15])] = np.array([float(t[10]), float(t[11]), float(t[12])])  # auth_seq_id, x y z
ns = sorted(cr)
P = np.array([af[s0 + n - 489][[a[0] for a in af[s0 + n - 489]].index("CA")][2] for n in ns])
Q = np.array([cr[n] for n in ns])
pc, qc = P.mean(0), Q.mean(0)
U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T; t = qc - R @ pc
rmsd = np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean())
lines = ["data_her2_domIV_modeled", "#", "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
         "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id",
         "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy",
         "_atom_site.B_iso_or_equiv", "_atom_site.auth_seq_id", "_atom_site.auth_comp_id", "_atom_site.auth_asym_id", "_atom_site.auth_atom_id",
         "_atom_site.pdbx_PDB_model_num"]
n = 0
for k in range(1, 120):
    for name, rn, x, b in af[s0 + k - 1]:
        n += 1
        x = x @ R.T + t
        lines.append(f"ATOM {n} {name[0]} {name} . {rn} A 1 {k} ? {x[0]:.3f} {x[1]:.3f} {x[2]:.3f} 1.00 {b:.2f} {k} {rn} A {name} 1")
open(f"{ROOT}/denovo/targets/her2_domIV_modeled.cif", "w").write("\n".join(lines) + "\n")
loop_plddt = np.mean([[a[3] for a in af[s0 + n - 489] if a[0] == "CA"][0] for n in range(581, 591)])
open(f"{ROOT}/denovo/targets/her2_domIV_modeled.README.txt", "w").write(
    f"AlphaFold DB P04626 domain IV (AF residues {s0}-{s0 + 118}) superposed on the 1N8Z HER2 CA atoms (RMSD {rmsd:.2f} A over {len(ns)} CA). "
    f"One chain A, resi k = HER2 number 488+k (489-607), contiguous, loop 581-590 modelled (mean pLDDT {loop_plddt:.0f}, low = flexible).\n")
print(f"wrote her2_domIV_modeled.cif: 119 residues, {n} atoms, superposition RMSD {rmsd:.2f} A, loop pLDDT {loop_plddt:.0f}")
