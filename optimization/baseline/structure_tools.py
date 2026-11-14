"""Compare predicted Fv-HER2 complexes with the 1N8Z crystal structure used by the Apheris HER2 example.

Predicted complex chains: A = VH (120), B = VL (107), C = HER2 domain IV (119 aa, residues 489-607).
Crystal chains (optimization/data/her2/her2_fv_1n8z.cif): B = VH, A = VL, C = HER2 489-580, D = HER2 591-607.
"""
import numpy as np
import pandas as pd

import common as C

CONTACT_CUTOFF = 4.5
CDRH3_RES = list(range(99, 109))  # VH numbering shared by crystal and predictions
PARATOPE = {"W99": 99, "G103": 103, "F104": 104, "Y105": 105}


def parse_cif(text):
    """Atom table from mmCIF text: chain, resseq, resname, atom, x, y, z, b."""
    hdr, rows = [], []
    for line in text.splitlines():
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1])
        elif line.startswith(("ATOM", "HETATM")):
            rows.append(line.split())
    df = pd.DataFrame(rows, columns=hdr[: len(rows[0])])
    seq = "label_seq_id" if "label_seq_id" in df and not (df.label_seq_id == ".").all() else "auth_seq_id"
    chain = "label_asym_id" if "label_asym_id" in df else "auth_asym_id"
    out = pd.DataFrame({
        "chain": df[chain], "res": df[seq].astype(int), "resname": df.label_comp_id,
        "atom": df.label_atom_id.str.strip('"'),
        "x": df.Cartn_x.astype(float), "y": df.Cartn_y.astype(float), "z": df.Cartn_z.astype(float),
        "b": df.B_iso_or_equiv.astype(float) if "B_iso_or_equiv" in df else np.nan,
    })
    return out[out.atom != "H"].reset_index(drop=True)


def kabsch(P, Q):
    """Rotation R and translation t minimising |R P + t - Q| (rows are points)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def _ca(df, chain, resids=None):
    d = df[(df.chain == chain) & (df.atom == "CA")]
    if resids is not None:
        d = d[d.res.isin(resids)]
    return d.set_index("res")


def load_reference():
    """Crystal reference re-labelled to the prediction convention (A=VH, B=VL, C=HER2 489-607 numbered 1-119)."""
    t = parse_cif(open(C.DATA / "her2" / "her2_fv_1n8z.cif").read())
    ref = t.copy()
    ref["chain"] = ref.chain.map({"B": "A", "A": "B", "C": "C", "D": "C"})
    ref.loc[ref.chain == "C", "res"] = ref.loc[ref.chain == "C", "res"] - 488
    return ref


def contact_pairs(df, cutoff=CONTACT_CUTOFF):
    """Set of (antibody chain, residue, HER2 residue) pairs with any heavy-atom distance below cutoff."""
    her = df[df.chain == "C"][["res", "x", "y", "z"]].to_numpy()
    pairs = set()
    for ch in ("A", "B"):
        ab = df[df.chain == ch][["res", "x", "y", "z"]].to_numpy()
        d = np.sqrt(((ab[:, None, 1:] - her[None, :, 1:]) ** 2).sum(-1))
        for i, j in zip(*np.where(d < cutoff)):
            pairs.add((ch, int(ab[i, 0]), int(her[j, 0])))
    return pairs


REF = None


def compare(pred):
    """Metrics of one predicted complex (atom table) against the crystal."""
    global REF
    if REF is None:
        REF = load_reference()
        REF.attrs["pairs"] = contact_pairs(REF)
    ref = REF
    out = {}
    # 1) align on HER2 residues present in the crystal
    rh = _ca(ref, "C")
    ph = _ca(pred, "C")
    common = sorted(set(rh.index) & set(ph.index))
    if len(common) < 50:
        return {"error": f"only {len(common)} HER2 residues matched"}
    R, t = kabsch(ph.loc[common, ["x", "y", "z"]].to_numpy(), rh.loc[common, ["x", "y", "z"]].to_numpy())
    out["her2_ca_rmsd_fit"] = rmsd(ph.loc[common, ["x", "y", "z"]].to_numpy() @ R.T + t, rh.loc[common, ["x", "y", "z"]].to_numpy())
    P = pred.copy()
    P[["x", "y", "z"]] = P[["x", "y", "z"]].to_numpy() @ R.T + t
    # 2) pose: Fv CA RMSD in the HER2 frame, and CDRH3 backbone RMSD
    fv_ref = pd.concat([_ca(ref, "A").assign(ch="A"), _ca(ref, "B").assign(ch="B")])
    fv_prd = pd.concat([_ca(P, "A").assign(ch="A"), _ca(P, "B").assign(ch="B")])
    key = lambda d: set(zip(d.ch, d.index))
    both = sorted(key(fv_ref) & key(fv_prd))
    a = np.array([fv_prd[(fv_prd.ch == c) & (fv_prd.index == r)][["x", "y", "z"]].to_numpy()[0] for c, r in both])
    b = np.array([fv_ref[(fv_ref.ch == c) & (fv_ref.index == r)][["x", "y", "z"]].to_numpy()[0] for c, r in both])
    out["fv_ca_rmsd_pose"] = rmsd(a, b)
    R2, t2 = kabsch(a, b)
    out["fv_ca_rmsd_own"] = rmsd(a @ R2.T + t2, b)

    def bb(df_, chain, resids):
        d = df_[(df_.chain == chain) & df_.res.isin(resids) & df_.atom.isin(["N", "CA", "C"])]
        return d.sort_values(["res", "atom"])[["x", "y", "z"]].to_numpy()
    out["cdrh3_bb_rmsd_pose"] = rmsd(bb(P, "A", CDRH3_RES), bb(ref, "A", CDRH3_RES))
    # 3) interface contacts reproduced (same residue pairs, 4.5 A)
    ref_pairs = ref.attrs["pairs"]
    pr_pairs = contact_pairs(P)
    hit = ref_pairs & pr_pairs
    out["contacts_crystal"] = len(ref_pairs)
    out["contacts_pred"] = len(pr_pairs)
    out["contact_recall"] = len(hit) / max(len(ref_pairs), 1)
    out["contact_precision"] = len(hit) / max(len(pr_pairs), 1)
    for name, r in PARATOPE.items():
        c_ref = {p for p in ref_pairs if p[0] == "A" and p[1] == r}
        c_prd = {p for p in pr_pairs if p[0] == "A" and p[1] == r}
        out[f"{name}_kept"] = len(c_ref & c_prd) / max(len(c_ref), 1)
    # 3b) epitope placement in the antibody frame: align the Fv, then compare the HER2 epitope residues.
    #     Robust to a mis-folded HER2 (the HER2-frame numbers above are not).
    Rf, tf = kabsch(np.array([fv_prd.loc[(fv_prd.ch == c) & (fv_prd.index == r)][["x", "y", "z"]].to_numpy()[0] for c, r in both]),
                    np.array([fv_ref.loc[(fv_ref.ch == c) & (fv_ref.index == r)][["x", "y", "z"]].to_numpy()[0] for c, r in both]))
    epi = sorted({p[2] for p in ref_pairs})
    pe = _ca(P, "C", epi)  # P: prediction already moved into the HER2 frame, the same frame fv_prd was built in
    re_ = _ca(ref, "C", epi)
    ok = sorted(set(pe.index) & set(re_.index))
    if ok:
        pe_xyz = pe.loc[ok, ["x", "y", "z"]].to_numpy() @ Rf.T + tf
        out["epitope_ca_rmsd_fvframe"] = rmsd(pe_xyz, re_.loc[ok, ["x", "y", "z"]].to_numpy())
    # 4) per-residue confidence of CDRH3 (B-factor column; pLDDT in these tools)
    h3 = pred[(pred.chain == "A") & pred.res.isin(CDRH3_RES) & (pred.atom == "CA")]
    out["cdrh3_plddt"] = float(h3.b.mean())
    return out


def interface_pae(pae, ref_pairs=None):
    """Mean PAE between crystal antibody-contact residues and HER2 contact residues (both directions).
    pae is the full matrix over chains A (120), B (107), C (119) in that order."""
    ref_pairs = ref_pairs or REF.attrs["pairs"]
    off = {"A": 0, "B": 120, "C": 227}
    vals = []
    for ch, ar, hr in ref_pairs:
        i, j = off[ch] + ar - 1, off["C"] + hr - 1
        vals.append(pae[i, j])
        vals.append(pae[j, i])
    vals = np.array(vals)
    return {"ipae_crystal_contacts": float(vals.mean()), "min_ipae": float(pae[:227, 227:].min())}


if __name__ == "__main__":
    # self-test: the crystal relabelled as a "prediction" must score perfectly
    ref = load_reference()
    r = compare(ref)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
