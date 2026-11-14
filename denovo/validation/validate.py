"""VAL: hard validity checks on every RFD4 (JapanFold rfd4) design. Writes results/denovo/validity.csv.

    denovo/.venv/bin/python denovo/validation/validate.py [--n N] [--workers 8]

Read-only on everything else. Checks (flag columns are val_*, 1 = problem unless named *_ok):
 1 table vs structure   2 fixed parts fixed   3 covalent geometry   4 disulfide   5 clashes + interface
 6 competition geometry vs trastuzumab Fv    7 sequence sanity     9 pose (orientation) of the VHH vs HER2
(8, self-consistency with NanoBodyBuilder2, is in selfconsist.py.)
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
S1 = ROOT / "denovo" / "s1_design"
TGT = ROOT / "denovo" / "targets"
OUT = ROOT / "results" / "denovo"

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
NHEAVY = {"G": 4, "A": 5, "S": 6, "C": 6, "V": 7, "T": 7, "P": 7, "L": 8, "I": 8, "N": 8, "D": 8, "Q": 9, "K": 9, "E": 9, "M": 8,
          "H": 10, "F": 11, "R": 11, "Y": 12, "W": 14}
STD_ATOMS = {"N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3",
             "CH2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "ND1", "ND2", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "SD",
             "SG", "OXT"}
BB = ("N", "CA", "C", "O")


# ----------------------------------------------------------------------------------------------- parsing
def read_cif(path):
    """-> dict chain -> list of residues {num, aa, atoms{name: xyz}}, plus element set and raw checks."""
    cols, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
            elif line.startswith(("ATOM", "HETATM")):
                rows.append(line.split())
    ix = {c: i for i, c in enumerate(cols)}
    seqcol = "label_seq_id" if rows and rows[0][ix["label_seq_id"]] not in (".", "?") else "auth_seq_id"
    chains, elems, bad = {}, set(), 0
    for t in rows:
        ch = t[ix["label_asym_id"]]
        num = int(t[ix["auth_seq_id"]] if seqcol == "auth_seq_id" else t[ix["label_seq_id"]])
        aa = THREE.get(t[ix["label_comp_id"]], "X")
        xyz = np.array([float(t[ix["Cartn_x"]]), float(t[ix["Cartn_y"]]), float(t[ix["Cartn_z"]])])
        elems.add(t[ix["type_symbol"]])
        if not np.isfinite(xyz).all() or not np.any(xyz):
            bad += 1
        res = chains.setdefault(ch, {})
        r = res.setdefault(num, {"num": num, "aa": aa, "atoms": {}})
        r["atoms"][t[ix["label_atom_id"]]] = xyz
    out = {ch: [res[k] for k in sorted(res)] for ch, res in chains.items()}
    return out, elems, bad


def seq_of(res):
    return "".join(r["aa"] for r in res)


def ca(res):
    return np.array([r["atoms"].get("CA", [np.nan] * 3) for r in res])


def kabsch(P, Q):
    """R, t with R @ P_i + t ~ Q_i (rows are points)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


# ----------------------------------------------------------------------------------------------- reference data
FW = {n: json.load(open(S1 / "frameworks" / f"{n}.json")) for n in ("bcii10", "capla")}
HER2_SEQ = "".join(l.strip() for l in open(TGT / "her2_domIV.fasta") if not l.startswith(">"))
HER2_A, HER2_C = HER2_SEQ[:92], HER2_SEQ[92:]
EPI = json.load(open(TGT / "epitope.json"))
EPI_RES = {r["auth_seq"] for r in EPI["residues"]}
REF, _, _ = read_cif(TGT / "ref_1n8z.cif")  # A=VL B=VH C,D=HER2
REF_HER2 = [r for r in REF["C"] if r["num"] <= 580] + [r for r in REF["D"] if r["num"] >= 591]
REF_FV = [r for ch in ("A", "B") for r in REF[ch]]
REF_HER2_CA = ca(REF_HER2)
REF_FV_XYZ = np.array([x for r in REF_FV for x in r["atoms"].values()])
REF_FV_TREE = cKDTree(REF_FV_XYZ)
INPUT_CACHE = {}


def load_input(name):
    if name not in INPUT_CACHE:
        INPUT_CACHE[name] = read_cif(S1 / "inputs" / name)[0]
    return INPUT_CACHE[name]


def input_name(fw, pose, caa):
    tag = "_caa" if (caa and fw == "capla") else ""
    return f"{fw}{tag}_on_her2" + ("" if pose == "p0" else f"_{pose}") + ".cif"


def parse_contig(contig):
    kept, bare = [], []
    for tok in contig.split(","):
        m = re.fullmatch(r"([A-Z])(\d+)-(\d+)", tok)
        if m:
            kept.append((m.group(1), int(m.group(2)), int(m.group(3))))
        elif re.fullmatch(r"\d+-\d+", tok):
            a, b = tok.split("-")
            bare.append((int(a), int(b)))
    return kept, bare


def chirality_sign(r):
    a = r["atoms"]
    if not all(k in a for k in ("N", "CA", "C", "CB")):
        return 0
    v = np.dot(np.cross(a["N"] - a["CA"], a["C"] - a["CA"]), a["CB"] - a["CA"])
    return int(np.sign(v))


# ----------------------------------------------------------------------------------------------- checks
def derive_segments(seq, fr_list, bare):
    """Locate FR1..FR4 in the designed chain from the sequence alone. Returns (ok, starts, cdrs, why)."""
    fr1, fr2, fr3, fr4 = fr_list
    if not seq.startswith(fr1):
        return False, None, None, "FR1 not at N-terminus"
    if not seq.endswith(fr4):
        return False, None, None, "FR4 not at C-terminus"
    (l1a, l1b), (l2a, l2b), (l3a, l3b) = bare
    s2 = [p for p in range(len(fr1) + l1a, len(fr1) + l1b + 1) if seq.startswith(fr2, p)]
    if len(s2) != 1:
        return False, None, None, f"FR2 found {len(s2)}x in window"
    s2 = s2[0]
    s3 = [p for p in range(s2 + len(fr2) + l2a, s2 + len(fr2) + l2b + 1) if seq.startswith(fr3, p)]
    if len(s3) != 1:
        return False, None, None, f"FR3 found {len(s3)}x in window"
    s3 = s3[0]
    s4 = len(seq) - len(fr4)
    n3 = s4 - (s3 + len(fr3))
    if not (l3a <= n3 <= l3b):
        return False, None, None, f"CDR3 block {n3} outside {l3a}-{l3b}"
    starts = [0, s2, s3, s4]
    cdr1 = seq[len(fr1):s2]
    cdr2 = seq[s2 + len(fr2):s3]
    cdr3 = seq[s3 + len(fr3):s4]
    return True, starts, (cdr1, cdr2, cdr3), ""


def analyse(args):
    row = args
    did = row["design_id"]
    out = {"design_id": did}
    P = json.loads(row["params_json"])
    fw, pose, caa = row["framework"], P["pose"], bool(P.get("caa_fixed", 0))
    try:
        ch, elems, nbad = read_cif(ROOT / row["design_cif"])
    except Exception as e:  # noqa: BLE001
        out.update(val_parse_fail=1, why=str(e)[:100])
        return out
    out["val_parse_fail"] = 0
    out["n_chains"] = len(ch)
    A, B, C = ch.get("A", []), ch.get("B", []), ch.get("C", [])
    seqB = seq_of(B)
    out["vhh_len_cif"] = len(B)
    # ---- 1 table vs structure
    out["val_chain_layout"] = int(not (sorted(ch) == ["A", "B", "C"] and len(A) == 92 and len(C) == 17))
    out["val_seq_mismatch"] = int(seqB != row["vhh_seq"])
    rec = FW[fw]
    fr_list = rec["caa"]["fr"] if caa else rec["fr"]
    kept, bare = parse_contig(P["contig"])
    kb = [k for k in kept if k[0] == "B"]
    ok, starts, cdrs, why = derive_segments(seqB, fr_list, bare)
    out["val_fr_not_intact"] = int(not ok)
    out["why_fr"] = why
    cdr_len = [np.nan] * 3
    if ok:
        cdr1, cdr2, cdr3 = cdrs
        cdr_len = [len(cdr1), len(cdr2), len(cdr3)]
        tab3 = row["cdr3"][2:] if caa else row["cdr3"]
        # the table cdr3 is the IMGT 105-117 style block: AA prefix + designed block (wave 2); wave 1 block incl. A105/A106
        out["val_cdr_table_mismatch"] = int(not (cdr1 == row["cdr1"] and cdr2 == row["cdr2"] and
                                                 ((("AA" + cdr3) == row["cdr3"]) if caa else (cdr3 == row["cdr3"]))))
        out["val_cdr_len_out_of_request"] = 0  # derive_segments enforces the windows
        out["cdr3_designed_len"] = len(cdr3)
        if not caa and fw == "bcii10":
            out["wave1_AA_changed"] = int(cdr3[:2] != "AA")
        else:
            out["wave1_AA_changed"] = 0
        # N-glyc in CDRs, homopolymers, repeats, low complexity
        out["nglyc_cdr"] = sum(len(re.findall(r"N[^P][ST]", c)) for c in (cdr1, cdr2, cdr3))
        out["homopoly4"] = int(any(re.search(r"(.)\1{3,}", c) for c in (cdr1, cdr2, cdr3)))
        out["dipep_repeat"] = int(any(re.search(r"(..)\1{2,}", c) for c in (cdr1, cdr2, cdr3)))
        out["cdr3_uniq_aa_frac"] = len(set(cdr3)) / max(1, len(cdr3))
        out["val_low_complexity"] = int(out["cdr3_uniq_aa_frac"] < 0.4 or out["homopoly4"])
        out["cdr3"] = cdr3
    else:
        out["val_cdr_table_mismatch"] = 1
    for i, n in enumerate(("cdr1_len", "cdr2_len", "cdr3_len")):
        out[n] = cdr_len[i]
    out["val_cys_unexpected"] = int(seqB.count("C") != 2)
    out["val_cys_odd"] = int(seqB.count("C") % 2 == 1)                                                # SEVERE: unpaired Cys (abguard fails it)
    # ---- 2 fixed parts fixed (HER2 sequence, HER2 + framework coordinates vs the submitted input)
    out["val_her2_seq_changed"] = int(seq_of(A) != HER2_A or seq_of(C) != HER2_C)
    inp = load_input(input_name(fw, pose, caa))
    in_her2 = inp["A"] + inp["C"]
    in_B = {r["num"]: r for r in inp["B"]}
    HA = A + C
    if len(HA) == len(in_her2) and ok:
        # map kept framework residues -> input residues via contig ranges
        pairs = []
        for (c, lo, hi), s in zip(kb, starts):
            for k in range(hi - lo + 1):
                if lo + k in in_B and s + k < len(B):
                    pairs.append((B[s + k], in_B[lo + k]))
        out["n_fw_mapped"] = len(pairs)
        out["val_fw_seq_changed"] = int(any(a["aa"] != b["aa"] for a, b in pairs))
        # kept group (HER2 + framework backbone) rigid superposition
        def bbxyz(rs):
            return np.array([[r["atoms"].get(a, [np.nan] * 3) for a in BB] for r in rs])
        o_bb = np.concatenate([bbxyz(HA).reshape(-1, 3), bbxyz([p[0] for p in pairs]).reshape(-1, 3)])
        i_bb = np.concatenate([bbxyz(in_her2).reshape(-1, 3), bbxyz([p[1] for p in pairs]).reshape(-1, 3)])
        m = np.isfinite(o_bb).all(1) & np.isfinite(i_bb).all(1)
        R, t = kabsch(o_bb[m], i_bb[m])
        out["kept_group_rmsd"] = rmsd((R @ o_bb[m].T).T + t, i_bb[m])
        # HER2 alone
        oc, ic = ca(HA), ca(in_her2)
        R2, t2 = kabsch(oc, ic)
        out["her2_ca_rmsd"] = rmsd((R2 @ oc.T).T + t2, ic)
        # framework pose relative to HER2 (superpose on HER2 only)
        ofw = bbxyz([p[0] for p in pairs]).reshape(-1, 3)
        ifw = bbxyz([p[1] for p in pairs]).reshape(-1, 3)
        mm = np.isfinite(ofw).all(1) & np.isfinite(ifw).all(1)
        out["fw_pose_rmsd_on_her2"] = rmsd(((R2 @ ofw[mm].T).T + t2), ifw[mm])
        out["val_her2_moved"] = int(out["her2_ca_rmsd"] > 0.5)
        out["val_fw_pose_moved"] = int(out["fw_pose_rmsd_on_her2"] > 1.5)
        out["val_motif_bad"] = int(out["kept_group_rmsd"] > 1.5)
    else:
        out.update(val_her2_moved=np.nan, val_fw_pose_moved=np.nan, val_motif_bad=np.nan, val_fw_seq_changed=np.nan)
    # ---- 3 covalent geometry (all chains)
    nbond_bad = nclash_cb = nomega = nchir = nmiss = nh = 0
    rama = {"fw": [0, 0], "cdr1": [0, 0], "cdr2": [0, 0], "cdr3": [0, 0]}
    for name, res in (("A", A), ("B", B), ("C", C)):
        n = len(res)
        if n < 3:
            continue
        get = lambda a: np.array([r["atoms"].get(a, [np.nan] * 3) for r in res])  # noqa: E731
        N_, CA_, C_, O_ = get("N"), get("CA"), get("C"), get("O")
        d = lambda a, b: np.linalg.norm(a - b, axis=1)  # noqa: E731
        bad = (np.abs(d(N_, CA_) - 1.46) > 0.12) | (np.abs(d(CA_, C_) - 1.52) > 0.12) | (np.abs(d(C_, O_) - 1.23) > 0.12)
        pep = d(C_[:-1], N_[1:])
        caca = d(CA_[:-1], CA_[1:])
        brk = caca > 4.3                      # a real chain break (the decoder does not enforce the C-N bond, see cn_* evidence)
        out[f"breaks_{name}"] = int(brk.sum())
        out[f"cn_sd_{name}"] = float(np.nanstd(pep))
        out[f"cn_min_{name}"] = float(np.nanmin(pep))
        out[f"cn_max_{name}"] = float(np.nanmax(pep))
        out[f"cn_outlier_frac_{name}"] = float(np.nanmean((pep < 1.25) | (pep > 1.41)))
        out[f"caca_min_{name}"] = float(np.nanmin(caca))
        out[f"caca_max_{name}"] = float(np.nanmax(caca))
        nbond_bad += int(np.nansum(bad)) + int(brk.sum())
        v1, v2 = N_ - CA_, C_ - CA_
        cosang = (v1 * v2).sum(1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
        ang = np.degrees(np.arccos(np.clip(cosang, -1, 1)))
        nbond_bad += int((np.abs(ang - 111) > 15).sum())
        omega = dihedral(CA_[:-1], C_[:-1], N_[1:], CA_[1:])
        nxt_pro = np.array([r["aa"] == "P" for r in res[1:]])
        nomega += int(((np.abs(omega) < 30) & ~nxt_pro & ~brk).sum())
        phi = np.full(n, np.nan)
        psi = np.full(n, np.nan)
        phi[1:] = dihedral(C_[:-1], N_[1:], CA_[1:], C_[1:])
        psi[:-1] = dihedral(N_[:-1], CA_[:-1], C_[:-1], N_[1:])
        aas = np.array([r["aa"] for r in res])
        region = np.array(["fw"] * n, dtype=object)
        if name == "B" and ok:
            s2, s3, s4 = starts[1], starts[2], starts[3]
            c1a, c1b = len(fr_list[0]), s2
            c2a, c2b = s2 + len(fr_list[1]), s3
            c3a, c3b = s3 + len(fr_list[2]), s4
            region[c1a:c1b], region[c2a:c2b], region[c3a:c3b] = "cdr1", "cdr2", "cdr3"
        left = (phi > 25) & (phi < 100) & (psi > -30) & (psi < 90)
        outl = ((phi > 0) & ~left & (aas != "G")) | ((aas == "P") & ((phi > -30) | (phi < -110)))
        valid = np.isfinite(phi) & np.isfinite(psi)
        for reg in rama:
            msk = (region == reg) & valid
            rama[reg][0] += int(outl[msk].sum())
            rama[reg][1] += int(msk.sum())
        for r in res:
            nh += int(any(a not in STD_ATOMS for a in r["atoms"]))
            nmiss += int(len([a for a in r["atoms"] if a != "OXT"]) < NHEAVY.get(r["aa"], 0))
            if name == "B":
                s = chirality_sign(r)
                nchir += int(s == -ref_sign_expect() and r["aa"] != "G")
    out["bond_angle_bad"] = nbond_bad
    out["cis_nonpro"] = nomega
    out["d_amino"] = nchir
    out["incomplete_res"] = nmiss
    out["odd_atoms"] = nh
    out["H_atoms"] = int("H" in elems or "D" in elems)
    out["zero_or_nan_xyz"] = nbad
    for reg, (o, t) in rama.items():
        out[f"rama_out_{reg}"] = o
        out[f"rama_n_{reg}"] = t
    out["val_geometry_bad"] = int(nomega > 0 or nchir > 0 or out["breaks_A"] + out["breaks_B"] + out["breaks_C"] > 0 or nbad > 0)   # SEVERE only; bond/angle outliers are evidence (bond_angle_bad)
    out["val_incomplete"] = int(nmiss > 0 or nh > 0 or out["H_atoms"])
    out["val_rama_cdr3_outlier_frac"] = rama["cdr3"][0] / max(1, rama["cdr3"][1])
    # ---- 4 disulfide (Cys at the end of FR1 and the Cys before CDR3)
    ss = np.nan
    if ok:
        i1 = len(fr_list[0]) - 1
        i2 = starts[2] + len(fr_list[2]) - (3 if caa else 1)
        if seqB[i1] == "C" and seqB[i2] == "C":
            ss = float(np.linalg.norm(B[i1]["atoms"]["SG"] - B[i2]["atoms"]["SG"])) if "SG" in B[i1]["atoms"] and "SG" in B[i2]["atoms"] else np.nan
    out["ss_sg_sg"] = ss
    out["ss_stretched"] = int(np.isfinite(ss) and 2.4 < ss <= 3.0)
    out["val_ss_not_formed"] = int(not (np.isfinite(ss) and ss <= 3.0))
    # HER2 disulfides formed
    cys = [r for r in HA if r["aa"] == "C" and "SG" in r["atoms"]]
    if cys:
        X = np.array([r["atoms"]["SG"] for r in cys])
        dm = np.linalg.norm(X[:, None] - X[None], axis=2)
        np.fill_diagonal(dm, 9)
        out["her2_ss_formed"] = int((dm.min(1) < 2.4).sum())
        out["her2_cys"] = len(cys)
    # ---- 5 clashes + interface
    allB = [(i, a, x) for i, r in enumerate(B) for a, x in r["atoms"].items()]
    XB = np.array([x for _, _, x in allB])
    IB = np.array([i for i, _, _ in allB])
    NB = [a for _, a, _ in allB]
    HAx = np.array([x for r in HA for x in r["atoms"].values()])
    HAi = np.array([k for k, r in enumerate(HA) for _ in r["atoms"]])
    tB = cKDTree(XB)
    pairs = tB.query_pairs(2.2, output_type="ndarray")
    if len(pairs):
        keep = np.abs(IB[pairs[:, 0]] - IB[pairs[:, 1]]) > 1
        sgsg = np.array([(NB[a] == "SG" and NB[b] == "SG") for a, b in pairs])
        out["clash_intra_vhh"] = int((keep & ~sgsg).sum())
    else:
        out["clash_intra_vhh"] = 0
    tH = cKDTree(HAx)
    dd, ii = tH.query(XB)
    out["clash_vhh_her2_lt2p5"] = int((dd < 2.5).sum())
    out["clash_vhh_her2_lt2p2"] = int((dd < 2.2).sum())
    out["min_dist_vhh_her2"] = float(dd.min())
    cont = dd < 4.5
    contacted_vhh = sorted(set(IB[cont]))
    out["vhh_res_in_contact"] = len(contacted_vhh)
    out["val_no_contact"] = int(len(contacted_vhh) == 0)
    out["val_clash_bad"] = int(out["clash_intra_vhh"] > 0 or out["clash_vhh_her2_lt2p2"] > 0)          # MILD: raw decoder output, bonds +-0.1 A
    out["val_clash_severe"] = int(out["clash_intra_vhh"] > 5 or out["min_dist_vhh_her2"] < 1.6)        # SEVERE: collapsed chain / interpenetration
    # HER2 residues contacted (auth numbering) and epitope/hotspot overlap
    pairs2 = tB.query_ball_tree(tH, 4.5)
    her2_idx = set()
    cdr_contacts = {"cdr1": 0, "cdr2": 0, "cdr3": 0, "fw": 0}
    if ok:
        bounds = {"cdr1": (len(fr_list[0]), starts[1]), "cdr2": (starts[1] + len(fr_list[1]), starts[2]),
                  "cdr3": (starts[2] + len(fr_list[2]), starts[3])}
    seen_res = {}
    for ai, hits in enumerate(pairs2):
        if hits:
            for h in hits:
                her2_idx.add(int(HAi[h]))
            seen_res[int(IB[ai])] = True
    if ok:
        for i in seen_res:
            tag = "fw"
            for k, (a, b) in bounds.items():
                if a <= i < b:
                    tag = k
            cdr_contacts[tag] += 1
        for k in ("cdr1", "cdr2", "cdr3"):
            out[f"contact_res_{k}"] = cdr_contacts[k]
        out["contact_res_fw"] = cdr_contacts["fw"]
    auth = []
    for k in her2_idx:
        n = HA[k]["num"]
        auth.append(488 + n if k < 92 else 590 + (k - 92) + 1)
    out["her2_res_contact"] = len(auth)
    out["epi_res_contact"] = len(set(auth) & EPI_RES)
    out["epi_precision"] = len(set(auth) & EPI_RES) / max(1, len(set(auth)))
    hot = {int(h[1:]) for h in P["hotspots"].split(",")} if isinstance(P.get("hotspots"), str) else set()
    out["hot_total"] = len(hot)
    out["hot_touched"] = len(hot & set(auth))
    if ok:
        c3a, c3b = bounds["cdr3"]
        cen3 = np.nanmean(ca(B[c3a:c3b]), axis=0)
        epi_ca = np.array([r["atoms"]["CA"] for k, r in enumerate(HA) if (488 + r["num"] if k < 92 else 590 + (k - 92) + 1) in EPI_RES])
        out["cdr3_cen_to_epi_cen"] = float(np.linalg.norm(cen3 - epi_ca.mean(0)))
        out["val_cdr3_no_epitope_contact"] = int(not any(a in EPI_RES for a in auth_of_cdr3(B, c3a, c3b, tH, HA)))
    # ---- 6 competition geometry: superpose design HER2 onto 1N8Z HER2, count VHH clashes with trastuzumab Fv
    R3, t3 = kabsch(ca(HA), REF_HER2_CA)
    XB_ref = (R3 @ XB.T).T + t3
    dfv, _ = REF_FV_TREE.query(XB_ref)
    out["fv_clash_lt2p5"] = int((dfv < 2.5).sum())
    out["fv_atoms_lt4"] = int((dfv < 4.0).sum())
    out["her2_ca_rmsd_to_1n8z"] = rmsd((R3 @ ca(HA).T).T + t3, REF_HER2_CA)
    # ---- 9 pose: rotation/translation of the submitted-pose framework relative to the p0 input framework, in the 1N8Z HER2 frame
    if ok and len(pairs) is not None:
        fwca = np.array([B[s + k]["atoms"]["CA"] for (c, lo, hi), s in zip(kb, starts) for k in range(hi - lo + 1) if s + k < len(B)])
        fwca_ref = (R3 @ fwca.T).T + t3
        p0 = load_input(input_name(fw, "p0", caa))
        R0, t0 = kabsch(ca(p0["A"] + p0["C"]), REF_HER2_CA)
        p0fw = np.array([r["atoms"]["CA"] for r in p0["B"] if any(lo <= r["num"] <= hi for (c, lo, hi) in kb)])
        p0fw_ref = (R0 @ p0fw.T).T + t0
        if len(p0fw_ref) == len(fwca_ref):
            R4, t4 = kabsch(p0fw_ref, fwca_ref)
            out["pose_rot_deg_vs_p0"] = float(np.degrees(np.arccos(np.clip((np.trace(R4) - 1) / 2, -1, 1))))
            out["pose_shift_vs_p0"] = float(np.linalg.norm(fwca_ref.mean(0) - p0fw_ref.mean(0)))
            out["pose_fw_rmsd_to_p0"] = rmsd(fwca_ref, p0fw_ref)
            vh = np.array([r["atoms"]["CA"] for r in REF["B"]])
            out["pose_cen_to_trast_vh"] = float(np.linalg.norm(fwca_ref.mean(0) - vh.mean(0)))
    return out


_REFSIGN = []


def ref_sign_expect():
    if not _REFSIGN:
        s = [chirality_sign(r) for r in REF_HER2 if r["aa"] != "G"]
        _REFSIGN.append(int(np.sign(np.sum(s))))
    return _REFSIGN[0]


def auth_of_cdr3(B, a, b, tH, HA):
    res = []
    for i in range(a, b):
        for x in B[i]["atoms"].values():
            for h in tH.query_ball_point(x, 4.5):
                pass
    # simple: any CDR3 atom within 4.5 A of an epitope residue atom
    xs = np.array([x for i in range(a, b) for x in B[i]["atoms"].values()])
    kk = set()
    for k, r in enumerate(HA):
        an = 488 + r["num"] if k < 92 else 590 + (k - 92) + 1
        if an in EPI_RES:
            hx = np.array(list(r["atoms"].values()))
            d = np.linalg.norm(xs[:, None] - hx[None], axis=2).min()
            if d < 4.5:
                kk.add(an)
    return kk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    D = pd.read_csv(OUT / "designs_rfd4jf.csv")
    if a.n:
        D = D.head(a.n)
    rows = D.to_dict("records")
    if a.workers > 1:
        with ProcessPoolExecutor(a.workers) as ex:
            res = list(ex.map(analyse, rows, chunksize=8))
    else:
        res = [analyse(r) for r in rows]
    V = pd.DataFrame(res)
    flags = [c for c in V.columns if c.startswith("val_") and c not in ("val_rama_cdr3_outlier_frac",)]
    hard = ["val_parse_fail", "val_chain_layout", "val_seq_mismatch", "val_fr_not_intact", "val_cdr_table_mismatch",
            "val_her2_seq_changed", "val_geometry_bad", "val_incomplete", "val_clash_severe", "val_no_contact", "val_motif_bad",
            "val_ss_not_formed", "val_cys_odd", "val_fw_seq_changed"]
    V["val_ok"] = (V[[h for h in hard if h in V]].fillna(1).sum(axis=1) == 0).astype(int)
    V.to_csv(OUT / "validity.csv" if not a.n else OUT / "validity_test.csv", index=False)
    print(len(V), "designs; val_ok =", int(V.val_ok.sum()))
    for f in flags:
        print(f"{f:32s} {int(np.nansum(V[f]))}")


if __name__ == "__main__":
    main()
