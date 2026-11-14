"""Geometry/sequence checks on one rfd4 design CIF (chains A,C = HER2 pieces renumbered 1.., B = VHH)."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
       "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
EPITOPE = json.load(open(HERE.parents[0] / "targets" / "epitope.json"))
EPI_RES = {r["auth_seq"] for r in EPITOPE["residues"]}


def read_cif(path):
    """-> {chain: [ {i (1-based order), aa3, atoms{name: xyz}} ]}, heavy atoms only."""
    hdr, ch = [], {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("_atom_site."):
            hdr.append(line.split(".")[1].strip())
        elif line.startswith(("ATOM", "HETATM")):
            d = dict(zip(hdr, line.split()))
            if d["type_symbol"] == "H" or d["label_comp_id"] not in AA3:
                continue
            c = d["label_asym_id"]
            r = ch.setdefault(c, {}).setdefault(int(d["label_seq_id"]), {"i": int(d["label_seq_id"]), "aa3": d["label_comp_id"], "atoms": {}})
            r["atoms"][d["label_atom_id"]] = np.array([float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])])
    return {c: [v[k] for k in sorted(v)] for c, v in ch.items()}


def seq_of(res):
    return "".join(AA3.get(r["aa3"], "X") for r in res)


def her2_auth(ch):
    """HER2 residues of an output (chain A 1..92 -> 489-580, chain C 1..17 -> 591-607) with their HER2 numbers."""
    out = []
    for r in ch.get("A", []):
        out.append({**r, "hnum": 488 + r["i"]})
    for r in ch.get("C", []):
        out.append({**r, "hnum": 590 + r["i"]})
    return out


def split_vhh(vseq, fr):
    fr1, fr2, fr3, fr4 = fr
    if not (vseq.startswith(fr1) and vseq.endswith(fr4)):
        return None
    mid = vseq[len(fr1):len(vseq) - len(fr4)]
    i = mid.find(fr2)
    if i < 0:
        return None
    c1, rest = mid[:i], mid[i + len(fr2):]
    j = rest.find(fr3)
    if j < 0:
        return None
    return {"cdr1": c1, "cdr2": rest[:j], "cdr3": rest[j + len(fr3):]}


def fr_mismatch(vseq, fr, parts):
    """number of framework positions that differ from the fixed framework (0 = held exactly)."""
    if parts is None:
        return -1
    rebuilt = fr[0] + parts["cdr1"] + fr[1] + parts["cdr2"] + fr[2] + parts["cdr3"] + fr[3]
    return 0 if rebuilt == vseq else -1


def interface(her2, vhh, vparts, hot_hnums):
    """Contacts (4.5 A heavy-atom) between HER2 and the VHH."""
    A = [(i, x) for i, r in enumerate(her2) for x in r["atoms"].values()]
    B = [(j, x) for j, r in enumerate(vhh) for x in r["atoms"].values()]
    XA = np.array([a[1] for a in A]); XB = np.array([b[1] for b in B])
    D = np.linalg.norm(XA[:, None] - XB[None], axis=2)
    close = D < 4.5
    ai = np.array([a[0] for a in A]); bj = np.array([b[0] for b in B])
    her_c = sorted({int(ai[k]) for k in np.where(close.any(1))[0]})
    vhh_c = sorted({int(bj[k]) for k in np.where(close.any(0))[0]})
    h_nums = {her2[i]["hnum"] for i in her_c}
    n1, n2, n3 = len(vparts["cdr1"]), len(vparts["cdr2"]), len(vparts["cdr3"])
    return {"contact_pairs": int(close.sum()), "clash_lt2p2": int((D < 2.2).sum()), "min_dist": round(float(D.min()), 2),
            "her2_res_contact": len(h_nums), "epi_res_contact": len(h_nums & EPI_RES),
            "epi_precision": round(len(h_nums & EPI_RES) / max(1, len(h_nums)), 3),
            "hot_touched": len(h_nums & set(hot_hnums)), "hot_total": len(hot_hnums),
            "vhh_res_contact": len(vhh_c), "cdr3_res_contact": 0, "contact_hnums": sorted(h_nums), "_vhh_c": vhh_c,
            "_cdr_lens": (n1, n2, n3)}


def geometry(cif, fr, hot_hnums, loop=False):
    """loop=True: HER2 is ONE chain A of 119 residues (resi k = HER2 488+k, modelled 581-590 loop included); else chains A (1-92) + C (1-17)."""
    ch = read_cif(cif)
    if not ({"A", "B"} if loop else {"A", "B", "C"}) <= set(ch):
        return {"ok": 0, "why": f"chains {sorted(ch)}"}
    vhh = ch["B"]
    vseq = seq_of(vhh)
    parts = split_vhh(vseq, fr)
    if parts is None:
        return {"ok": 0, "why": "framework not held", "vhh_seq": vseq}
    her2 = [{**r, "hnum": 488 + r["i"]} for r in ch["A"]] if loop else her2_auth(ch)
    g = interface(her2, vhh, parts, hot_hnums)
    n1, n2, n3 = g.pop("_cdr_lens")
    c = g.pop("_vhh_c")
    s1 = len(fr[0]); s2 = s1 + n1 + len(fr[1]); s3 = s2 + n2 + len(fr[2])
    g["cdr1_res_contact"] = sum(1 for k in c if s1 <= k < s1 + n1)
    g["cdr2_res_contact"] = sum(1 for k in c if s2 <= k < s2 + n2)
    g["cdr3_res_contact"] = sum(1 for k in c if s3 <= k < s3 + n3)
    g["contact_hnums"] = ",".join(map(str, g["contact_hnums"]))
    g.update({"ok": 1, "vhh_seq": vseq, **parts})
    return g
