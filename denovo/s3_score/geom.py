"""Structure / PAE metrics for VHH (or Fv) + HER2 domain IV complexes. numpy only.

Prediction convention (everything in S3 uses it):
    VHH complex : chain A = VHH (L residues), chain B = HER2 domain IV (119 aa, label 1..119 = HER2 489..607)
    Fv complex  : chain A = VH, chain B = VL, chain C = HER2 domain IV   (the trastuzumab control)
Reference crystal: denovo/targets/ref_1n8z.cif, chains A=VL, B=VH, C,D=HER2 (label_seq_id == HER2 residue number).
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "denovo" / "targets"
HER2_OFFSET = 488            # prediction label i  <->  HER2 residue 488 + i
CONTACT_A = 4.5              # heavy-atom distance defining a residue contact
IPSAE_PAE_CUTOFF = 10.0      # A, as in Dunbrack's default
HER2_LEN = 119


# ------------------------------------------------------------------ parsing
class Atoms:
    """Heavy atoms of one structure: parallel arrays."""

    def __init__(self, chain, res, resname, atom, xyz, b):
        self.chain, self.res, self.resname, self.atom = map(np.asarray, (chain, res, resname, atom))
        self.xyz, self.b = np.asarray(xyz, float), np.asarray(b, float)

    def sel(self, mask):
        return Atoms(self.chain[mask], self.res[mask], self.resname[mask], self.atom[mask], self.xyz[mask], self.b[mask])

    def of_chain(self, c):
        return self.sel(self.chain == c)

    def chains(self):
        seen = []
        for c in self.chain:
            if c not in seen:
                seen.append(c)
        return seen

    def ca(self):
        return self.sel(self.atom == "CA")

    def seq(self, c):
        THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                 "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                 "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
        a = self.of_chain(c).ca()
        return "".join(THREE.get(r, "X") for r in a.resname)


def parse_cif(text: str) -> Atoms:
    hdr, rows = [], []
    for line in text.splitlines():
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1])
        elif line.startswith(("ATOM", "HETATM")):
            rows.append(line.split())
    if not rows:
        raise ValueError("no atoms")
    ix = {h: i for i, h in enumerate(hdr)}
    ch = ix.get("label_asym_id", ix.get("auth_asym_id"))
    rs = ix["label_seq_id"] if "label_seq_id" in ix and rows[0][ix["label_seq_id"]] != "." else ix["auth_seq_id"]
    keep = [r for r in rows if r[ix["type_symbol"]] != "H"]
    if "pdbx_PDB_model_num" in ix:
        keep = [r for r in keep if r[ix["pdbx_PDB_model_num"]] == "1"]
    return Atoms([r[ch] for r in keep], [int(r[rs]) for r in keep], [r[ix["label_comp_id"]] for r in keep],
                 [r[ix["label_atom_id"]].strip('"') for r in keep],
                 [[float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]])] for r in keep],
                 [float(r[ix["B_iso_or_equiv"]]) if "B_iso_or_equiv" in ix else 0.0 for r in keep])


def load_cif(path) -> Atoms:
    return parse_cif(Path(path).read_text())


# ------------------------------------------------------------------ geometry
def kabsch(P, Q):
    """R, t minimising |R P + t - Q| for matched rows."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def dist_matrix(A, B):
    return np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))


def her2_frame(atoms: Atoms, her_chain: str, ref_ca_xyz: dict):
    """(R, t, n) moving `atoms` into the crystal HER2 frame by CA superposition on shared residues."""
    ca = atoms.of_chain(her_chain).ca()
    ids = [r for r in ca.res if (r + HER2_OFFSET) in ref_ca_xyz]
    if len(ids) < 40:
        raise ValueError(f"only {len(ids)} HER2 residues shared with crystal")
    P = np.array([ca.xyz[ca.res == r][0] for r in ids])
    Q = np.array([ref_ca_xyz[r + HER2_OFFSET] for r in ids])
    R, t = kabsch(P, Q)
    return R, t, rmsd(P @ R.T + t, Q)


@lru_cache(maxsize=1)
def reference():
    """Crystal Fv + HER2 (heavy atoms), HER2 CA table, epitope sets in prediction numbering."""
    ref = load_cif(TARGETS / "ref_1n8z.cif")
    her = ref.sel(np.isin(ref.chain, ["C", "D"])).ca()
    her_ca = {int(r): x for r, x in zip(her.res, her.xyz)}
    fv = ref.sel(np.isin(ref.chain, ["A", "B"]))
    ep = json.loads((TARGETS / "epitope.json").read_text())
    epitope = sorted(int(r["auth_seq"]) - HER2_OFFSET for r in ep["residues"])
    hot = {k: sorted(int(h["auth_seq"]) - HER2_OFFSET for h in v) for k, v in ep["hotspot_sets"].items()}
    return {"fv": fv, "her_ca": her_ca, "epitope": epitope, "hotspots": hot}


# ------------------------------------------------------------------ sequence helpers
def find_cdrs(vhh_seq: str, cdr1="", cdr2="", cdr3=""):
    """Index ranges (0-based, end exclusive) of CDR1-3 in the VHH; use the design table when given, else heuristics."""
    out = {}
    for name, s in (("cdr1", cdr1), ("cdr2", cdr2), ("cdr3", cdr3)):
        i = vhh_seq.find(s) if s else -1
        if i >= 0:
            out[name] = (i, i + len(s))
    if "cdr3" not in out:  # between the conserved 'YxC' and 'WGxG'
        ms = list(re.finditer(r"Y[YFHW]C", vhh_seq[:115]))
        if ms:
            st = ms[-1].end()
            w = re.search(r"W[GS].G", vhh_seq[st:])
            if w:
                out["cdr3"] = (st, st + w.start())
    if "cdr1" not in out and len(vhh_seq) > 40:
        out["cdr1"] = (26, 33)
    if "cdr2" not in out and len(vhh_seq) > 70:
        out["cdr2"] = (50, 58)
    return out


# ------------------------------------------------------------------ PAE metrics
def _d0(L: float) -> float:
    L = max(L, 27.0)
    return max(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def ipsae(pae: np.ndarray, ia: slice, ib: slice, cutoff=IPSAE_PAE_CUTOFF, d0chn: bool = False) -> float:
    """ipSAE (Dunbrack 2025, d0res variant): for each residue i of chain A, average 1/(1+(PAE_ij/d0)^2) over
    residues j of chain B with PAE_ij < cutoff, d0 from the number of such residues (L>=27); take the max over i.
    Symmetrised by taking the larger of A->B and B->A."""
    d0c = _d0((ia.stop - ia.start) + (ib.stop - ib.start))   # chain-length d0 (Dunbrack's ipSAE_d0chn), ~5-6 A for a VHH + domain
    def one(M):
        best = 0.0
        for i in range(M.shape[0]):
            ok = M[i] < cutoff
            n = int(ok.sum())
            if n == 0:
                continue
            d0 = d0c if d0chn else _d0(n)
            best = max(best, float((1.0 / (1.0 + (M[i][ok] / d0) ** 2)).mean()))
        return best
    return max(one(pae[ia, ib]), one(pae[ib, ia].T))


def pae_metrics(pae: np.ndarray, n_vhh: int, contact_pairs: set, cdr_idx: list, epi_idx: list, n_her: int = HER2_LEN) -> dict:
    v, h = slice(0, n_vhh), slice(n_vhh, n_vhh + n_her)
    cross = np.concatenate([pae[v, h].ravel(), pae[h, v].ravel()])
    out = {"min_ipae": float(cross.min()), "mean_ipae": float(cross.mean()), "ipsae": ipsae(pae, v, h),
           # the 10 A cutoff gives exactly 0 when the best interface PAE is above 10 A; wider cutoffs keep a ranking signal
           "ipsae15": ipsae(pae, v, h, 15.0), "ipsae20": ipsae(pae, v, h, 20.0),
           # chain-length d0 variants: non-zero when the best interface PAE is 10-15 A (e.g. the trastuzumab control)
           "ipsae_chn10": ipsae(pae, v, h, 10.0, True), "ipsae_chn15": ipsae(pae, v, h, 15.0, True)}
    if contact_pairs:
        vals = [0.5 * (pae[i, n_vhh + j - 1] + pae[n_vhh + j - 1, i]) for i, j in contact_pairs]
        out["ipae_contacts"] = float(np.mean(vals))
        # PAE-weighted contact score: TM-style 1/(1+(PAE/8 A)^2) averaged over contacting residue pairs, 0..1, no cutoff,
        # so it keeps resolving between 10 and 25 A where ipSAE (d0 ~ 1-2 A) is flat at 0
        out["ipae_tm"] = float(np.mean([1.0 / (1.0 + (v / 8.0) ** 2) for v in vals]))
    else:
        out["ipae_contacts"] = out["ipae_tm"] = float("nan")
    if cdr_idx and epi_idx:
        out["ipae_cdr_epitope"] = float(np.mean([0.5 * (pae[i, n_vhh + j - 1] + pae[n_vhh + j - 1, i])
                                                 for i in cdr_idx for j in epi_idx]))
    else:
        out["ipae_cdr_epitope"] = float("nan")
    return out


# ------------------------------------------------------------------ complex metrics
def contact_pairs(vhh: Atoms, her: Atoms, cutoff=CONTACT_A):
    """{(vhh 0-based residue index, her2 label 1..119)} with any heavy-atom distance < cutoff, plus residue sets."""
    D = dist_matrix(vhh.xyz, her.xyz)
    ii, jj = np.where(D < cutoff)
    first = int(vhh.res.min())
    return {(int(vhh.res[i]) - first, int(her.res[j])) for i, j in zip(ii, jj)}


def complex_metrics(atoms: Atoms, cdrs: dict, vhh_chain="A", her_chain="B", extra_ab_chains=()) -> dict:
    """Interface geometry of one model: contacts, epitope overlap, trastuzumab competition proxy, pLDDT."""
    ref = reference()
    vhh, her = atoms.of_chain(vhh_chain), atoms.of_chain(her_chain)
    if extra_ab_chains:  # Fv control: pool VH + VL as the "binder"
        vhh = atoms.sel(np.isin(atoms.chain, [vhh_chain, *extra_ab_chains]))
    pairs = contact_pairs(vhh if not extra_ab_chains else atoms.of_chain(vhh_chain), her)
    if extra_ab_chains:
        for c in extra_ab_chains:
            pairs |= {(i + 10000, j) for i, j in contact_pairs(atoms.of_chain(c), her)}
    her_hit = {j for _, j in pairs}
    vhh_hit = {i for i, _ in pairs}
    epi = set(ref["epitope"])
    out = {"n_contact_pairs": len(pairs), "n_vhh_res_contact": len(vhh_hit), "n_her2_res_contact": len(her_hit),
           "epi_overlap": len(her_hit & epi) / max(len(her_hit), 1),     # of what the VHH touches, share inside trastuzumab epitope
           "epi_coverage": len(her_hit & epi) / max(len(epi), 1)}        # of the epitope, share the VHH touches
    for k, v in ref["hotspots"].items():
        out[f"cov_{k}"] = len(her_hit & set(v)) / max(len(v), 1)
    # binder pLDDT (B column, 0-100)
    ca = vhh.ca()
    bsc = 100.0 if ca.b.max() <= 1.0 else 1.0
    out["binder_plddt"] = float(ca.b.mean() * bsc)
    if "cdr3" in cdrs:
        a, b = cdrs["cdr3"]
        out["cdr3_plddt"] = float(ca.b[a:b].mean() * bsc) if b > a else float("nan")
    # competition proxy: move the binder into the crystal HER2 frame, test against trastuzumab Fv
    try:
        R, t, fit = her2_frame(atoms, her_chain, ref["her_ca"])
        X = vhh.xyz @ R.T + t
        D = dist_matrix(X, ref["fv"].xyz)
        dmin = D.min(1)
        out["her2_fit_rmsd"] = fit
        out["clash25"] = int((dmin < 2.5).sum())
        out["clash30"] = int((dmin < 3.0).sum())
        out["fv_contact40"] = int((dmin < 4.0).sum())
        cdr_mask = np.zeros(len(vhh.res), bool)
        first = int(vhh.res.min())
        for a, b in cdrs.values():
            cdr_mask |= (vhh.res - first >= a) & (vhh.res - first < b)
        out["cdr_fv_frac4"] = float((dmin[cdr_mask] < 4.0).mean()) if cdr_mask.any() else float("nan")
    except Exception as e:  # pragma: no cover
        out["comp_error"] = str(e)[:80]
    return out


def pose_vs(atoms_a: Atoms, atoms_b: Atoms, cdr3, vhh_chain="A", her_chain="B"):
    """Pose difference of two predictions of the same sequences: VHH CA RMSD and CDR3 centroid distance after
    superposing on HER2 CA (all 119 residues); also the VHH-internal RMSD (own fold agreement)."""
    ha, hb = atoms_a.of_chain(her_chain).ca(), atoms_b.of_chain(her_chain).ca()
    R, t = kabsch(ha.xyz, hb.xyz)
    va, vb = atoms_a.of_chain(vhh_chain).ca(), atoms_b.of_chain(vhh_chain).ca()
    n = min(len(va.xyz), len(vb.xyz))
    A, B = va.xyz[:n] @ R.T + t, vb.xyz[:n]
    out = {"vhh_rmsd_her2frame": rmsd(A, B)}
    R2, t2 = kabsch(va.xyz[:n], vb.xyz[:n])
    out["vhh_rmsd_own"] = rmsd(va.xyz[:n] @ R2.T + t2, vb.xyz[:n])
    if cdr3:
        a, b = cdr3
        out["cdr3_centroid_dist"] = float(np.linalg.norm(A[a:b].mean(0) - B[a:b].mean(0)))
    return out


def own_design_rmsd(design_cif, pred: Atoms, her_seq: str, vhh_chain="A", her_chain="B"):
    """Refold fidelity: VHH CA RMSD between the generator's own design and the prediction, in the HER2 frame
    (pose) and after superposing on the VHH itself (fold). HER2 / VHH chains of the design are found by sequence."""
    from difflib import SequenceMatcher
    d = load_cif(design_cif)
    seqs = {c: d.seq(c) for c in d.chains()}
    her_c = max(seqs, key=lambda c: SequenceMatcher(None, seqs[c], her_seq, autojunk=False).ratio())
    others = [c for c in seqs if c != her_c and len(seqs[c]) > 80]
    if not others:
        return {}
    vc = max(others, key=lambda c: len(seqs[c]))
    pv, ph = pred.of_chain(vhh_chain).ca(), pred.of_chain(her_chain).ca()
    dv, dh = d.of_chain(vc).ca(), d.of_chain(her_c).ca()

    def pairs(sa, sb):
        sm = SequenceMatcher(None, sa, sb, autojunk=False)
        return [(a + k, b + k) for a, b, n in sm.get_matching_blocks() for k in range(n)]
    hp, vp = pairs(dh_seq := "".join(seqs[her_c]), "".join(pred.seq(her_chain))), pairs(seqs[vc], pred.seq(vhh_chain))
    if len(hp) < 40 or len(vp) < 60:
        return {}
    R, t = kabsch(ph.xyz[[b for _, b in hp]], dh.xyz[[a for a, _ in hp]])
    P = pv.xyz[[b for _, b in vp]] @ R.T + t
    Q = dv.xyz[[a for a, _ in vp]]
    out = {"own_vhh_rmsd_her2frame": rmsd(P, Q)}
    R2, t2 = kabsch(pv.xyz[[b for _, b in vp]], Q)
    out["own_vhh_rmsd_fold"] = rmsd(pv.xyz[[b for _, b in vp]] @ R2.T + t2, Q)
    return out


def to_pdb(chains) -> str:
    """PDB text from [(chain id, Atoms)], heavy atoms only (AF2-IG wants the designed complex as coordinates)."""
    lines, n = [], 0
    for cid, a in chains:
        order = {r: k + 1 for k, r in enumerate(dict.fromkeys(a.res.tolist()))}   # renumber 1..N (design files have gaps / offsets)
        for i in range(len(a.res)):
            n += 1
            name = str(a.atom[i])
            nm = f" {name:<3}" if len(name) < 4 else name
            x, y, z = a.xyz[i]
            lines.append(f"ATOM  {n:5d} {nm} {a.resname[i]:>3} {cid}{order[int(a.res[i])]:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{min(a.b[i], 99.99):6.2f}          {name[0]:>2}")
        lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def split_complex(path, her_seq: str):
    """(HER2 Atoms pooled from every chain that is a piece of `her_seq`, VHH Atoms) of a design/prediction file.
    RFD4 designs carry HER2 as two segments (489-580, 591-607); they are pooled in sequence order."""
    from difflib import SequenceMatcher
    d = load_cif(path)
    seqs = {c: d.seq(c) for c in d.chains()}
    def frac(c):
        m = SequenceMatcher(None, seqs[c], her_seq, autojunk=False).find_longest_match(0, len(seqs[c]), 0, len(her_seq))
        return m.size / max(len(seqs[c]), 1), m.b
    her_cs = sorted([c for c in seqs if len(seqs[c]) >= 10 and frac(c)[0] >= 0.9], key=lambda c: frac(c)[1])
    others = [c for c in seqs if c not in her_cs and len(seqs[c]) > 80]
    if not her_cs or not others:
        raise ValueError("HER2 / VHH chains not found")
    her = d.sel(np.isin(d.chain, her_cs))
    order = {c: k for k, c in enumerate(her_cs)}
    idx = np.argsort([order[c] for c in her.chain], kind="stable")
    her = Atoms(her.chain[idx], her.res[idx], her.resname[idx], her.atom[idx], her.xyz[idx], her.b[idx])
    return her, d.of_chain(max(others, key=lambda c: len(seqs[c])))
