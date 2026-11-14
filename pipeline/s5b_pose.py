"""Stage 5b: does the design bind HER2 where the crystal does? RMSD to 1N8Z after superposing the Fv framework.

All predictions are compared with the 1N8Z crystal on the standardised, renumbered structures written by the contacts
stage (chains H, L, G; HER2 numbered like the crystal). Each prediction is superposed on the crystal using the C-alpha
atoms of the Fv framework (VH without CDR-H3, plus VL), then:

  fv_fit_rmsd   residual of that superposition (sanity)
  h3_rmsd       C-alpha RMSD of CDR-H3 (99-108)
  her2_rmsd     C-alpha RMSD of all HER2 residues resolved in the crystal   -> where HER2 sits relative to the Fv
  epitope_rmsd  C-alpha RMSD of the crystal epitope (HER2 residues with any heavy atom within 5 A of the crystal Fv)
                                                                              -> same epitope?
HER2 residues 581-590 are unresolved in the crystal and are ignored.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import out, save_summary


def read_ca(pdb) -> tuple[dict, np.ndarray, list]:
    """CA per (chain, resnum) and the full heavy-atom array with residue keys."""
    ca, atoms, keys = {}, [], []
    for line in open(pdb):
        if line.startswith("ATOM") and line[76:78].strip() != "H":
            key = (line[21], int(line[22:26]))
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            atoms.append(xyz)
            keys.append(key)
            if line[12:16].strip() == "CA":
                ca[key] = xyz
    return ca, np.array(atoms), keys


def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t that map P onto Q (least squares)."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def rmsd(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def crystal_epitope(cry_atoms: np.ndarray, cry_keys: list, cutoff: float = 5.0) -> list:
    fv = np.array([a for a, k in zip(cry_atoms, cry_keys) if k[0] in "HL"])
    epi = []
    for key in sorted({k for k in cry_keys if k[0] == "G"}):
        g = np.array([a for a, k in zip(cry_atoms, cry_keys) if k == key])
        if np.linalg.norm(g[:, None, :] - fv[None, :, :], axis=2).min() <= cutoff:
            epi.append(key)
    return epi


def run(cfg: dict) -> None:
    pdbs = out(cfg, "contacts") / "pdb"
    cry_ca, cry_atoms, cry_keys = read_ca(pdbs / "crystal.pdb")
    epitope = [k for k in crystal_epitope(cry_atoms, cry_keys) if k in cry_ca]
    her2 = [k for k in cry_ca if k[0] == "G"]
    frame = [k for k in cry_ca if (k[0] == "L") or (k[0] == "H" and not 99 <= k[1] <= 108)]
    h3 = [k for k in cry_ca if k[0] == "H" and 99 <= k[1] <= 108]
    print(f"  crystal epitope: {len(epitope)} HER2 residues ({epitope[0][1]}-{epitope[-1][1]}); framework {len(frame)} CA")
    rows = []
    for f in sorted(pdbs.glob("*.pdb")):
        if f.stem == "crystal":
            continue
        ca, _, _ = read_ca(f)
        if not all(k in ca for k in frame + h3 + her2):
            continue
        R, t = kabsch(np.array([ca[k] for k in frame]), np.array([cry_ca[k] for k in frame]))
        mv = lambda keys: (np.array([ca[k] for k in keys]) @ R.T) + t
        ref = lambda keys: np.array([cry_ca[k] for k in keys])
        rows.append(dict(id=f.stem, fv_fit_rmsd=rmsd(mv(frame), ref(frame)), h3_rmsd=rmsd(mv(h3), ref(h3)),
                         her2_rmsd=rmsd(mv(her2), ref(her2)), epitope_rmsd=rmsd(mv(epitope), ref(epitope))))
    df = pd.DataFrame(rows)
    df.to_csv(out(cfg, "pose.csv"), index=False)
    print(df[df.id.str.startswith("ctrl_")].round(2).to_string(index=False))
    save_summary(cfg, n_epitope_residues=len(epitope))
