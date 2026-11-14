"""Build the rfd4 inputs: HER2 domain IV (chains A and C, HER2 numbering, 1N8Z frame) + a VHH framework (chain B)
pre-placed where trastuzumab's VH sits. Writes frameworks/<fw>.json and inputs/<fw>_on_her2.cif.

Run: baseline/.venv/bin/python denovo/s1_design/prep_inputs.py [--pose-variants]

Why pre-placed: the hosted rfd4 wrapper puts HER2 and the kept framework backbone in ONE rigid coordinate group
(motif_rmsd of designs vs the submitted structure is ~0.4 A), so the framework pose relative to HER2 is whatever we
submit. We submit it where trastuzumab's VH is in 1N8Z (framework CA superposition), so CDR3 lands on the epitope.

CDR convention follows the RFD4 antibody examples on the VM (CDR = everything between the conserved Cys and the
FR2 Trp, and so on): CDR1 = IMGT 24-40, CDR2 = IMGT 56-65, CDR3 = IMGT 105-117.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "denovo"))
import immunum  # noqa: E402

FR = HERE / "frameworks"
INP = HERE / "inputs"
INP.mkdir(exist_ok=True)
T = ROOT / "denovo" / "targets"
CDR_IMGT = {"cdr1": (24, 40), "cdr2": (56, 65), "cdr3": (105, 117)}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
ANN = immunum.Annotator(chains=["H"], scheme="imgt", min_confidence=0.0)
# CDR length ranges requested from the model (min,max residues in the block between the framework pieces)
CDR_LEN = {"bcii10": {"cdr1": (13, 17), "cdr2": (8, 10)}, "capla": {"cdr1": (12, 15), "cdr2": (9, 11)}}
CDR3_RANGES = {"c3a": (12, 16), "c3b": (16, 21)}


def read_atoms(path, chain, label_seq=True):
    hdr, res = [], {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("_atom_site."):
            hdr.append(line.strip().split(".")[1])
        elif line.startswith("ATOM"):
            d = dict(zip(hdr, line.split()))
            if d["label_asym_id"] != chain or d["label_comp_id"] not in AA3:
                continue
            if d.get("label_alt_id", ".") not in (".", "A") or d.get("type_symbol") == "H":
                continue
            k = int(d["label_seq_id"]) if label_seq else int(d["auth_seq_id"])
            r = res.setdefault(k, {"num": k, "aa3": d["label_comp_id"], "atoms": {}})
            r["atoms"].setdefault(d["label_atom_id"].strip('"'), (float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])))
    return [res[k] for k in sorted(res)]


def imgt_int(lab):
    return int("".join(c for c in lab if c.isdigit()))


def imgt_map(seq):
    r = ANN.number(seq)
    assert r.error is None, r.error
    labels = list(r.numbering.keys())
    vals = "".join(r.numbering.values())
    i0 = seq.find(vals[:12])
    assert seq[i0:i0 + len(vals)] == vals, "numbering is not contiguous in the query"
    return {i0 + i: lab for i, lab in enumerate(labels)}, r


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def write_cif(path, chains, name):
    L = ["data_" + name, "#", "loop_"]
    for c in ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
              "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
              "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]:
        L.append("_atom_site." + c)
    n = 0
    for ei, (cid, residues) in enumerate(chains, 1):
        for r in residues:
            for an, (x, y, z) in r["atoms"].items():
                n += 1
                L.append(f"ATOM {n} {an[0]} {an} . {r['aa3']} {cid} {ei} {r['num']} ? {x:.3f} {y:.3f} {z:.3f} 1.00 0.00 "
                         f"{r['num']} {r['aa3']} {cid} {an} 1")
    Path(path).write_text("\n".join(L) + "\n")


def pieces(res, cuts):
    """Split residues into FR pieces at the given (start,end) label ranges of the framework. Returns list of (lo,hi)."""
    nums = [r["num"] for r in res]
    return [(lo, hi) for lo, hi in cuts if lo in nums and hi in nums]


def build_record(name, res, fr_ranges, stub_cdrs=None):
    """fr_ranges: [(lo,hi)]*4 residue-number ranges of FR1..FR4 in `res` (CDRs between them)."""
    bynum = {r["num"]: r for r in res}
    seqs = ["".join(AA3[bynum[n]["aa3"]] for n in range(lo, hi + 1)) for lo, hi in fr_ranges]
    if stub_cdrs is None:  # CDR residues exist in the structure: number the real sequence
        full = "".join(AA3[r["aa3"]] for r in res)
        m, nr = imgt_map(full)
        idx = {r["num"]: i for i, r in enumerate(res)}
        imgt_of = {n: m[idx[n]] for n in idx}
        native = {}
        for cdr, (a, b) in CDR_IMGT.items():
            ns = [n for n, lab in imgt_of.items() if a <= imgt_int(lab) <= b]
            native[cdr] = "".join(AA3[bynum[n]["aa3"]] for n in range(min(ns), max(ns) + 1))
        conf = nr.confidence
    else:  # framework-only structure: number FR1+stub+FR2+stub+FR3+stub+FR4 and read the FR labels back
        full = seqs[0] + stub_cdrs[0] + seqs[1] + stub_cdrs[1] + seqs[2] + stub_cdrs[2] + seqs[3]
        m, nr = imgt_map(full)
        imgt_of, off = {}, 0
        for k, (lo, hi) in enumerate(fr_ranges):
            for j, n in enumerate(range(lo, hi + 1)):
                imgt_of[n] = m[off + j]
            off += (hi - lo + 1) + (len(stub_cdrs[k]) if k < 3 else 0)
        native, conf = None, nr.confidence
    contig_fr = [f"B{lo}-{hi}" for lo, hi in fr_ranges]
    return {"name": name, "fr": seqs, "fr_ranges": fr_ranges, "contig_fr": contig_fr, "native_cdr": native,
            "imgt_confidence": conf, "imgt_of_resnum": {str(n): l for n, l in imgt_of.items()},
            "cdr_def": "CDR1=IMGT24-40 CDR2=IMGT56-65 CDR3=IMGT105-117 (as the RFD4 antibody examples)",
            "cdr_len": CDR_LEN[name], "cdr3_ranges": CDR3_RANGES}


def place_on_trastuzumab(res, rec, extra=None):
    ref = read_atoms(T / "ref_1n8z.cif", "B")
    rseq = "".join(AA3.get(r["aa3"], "X") for r in ref)
    rmap, _ = imgt_map(rseq)
    by_lab = {lab: ref[i] for i, lab in rmap.items()}
    P, Q = [], []
    for r in res:
        lab = rec["imgt_of_resnum"].get(str(r["num"]))
        if lab is None or lab not in by_lab:
            continue
        k = imgt_int(lab)
        if any(a <= k <= b for a, b in CDR_IMGT.values()):
            continue
        if "CA" in r["atoms"] and "CA" in by_lab[lab]["atoms"]:
            P.append(r["atoms"]["CA"]); Q.append(by_lab[lab]["atoms"]["CA"])
    P, Q = np.array(P), np.array(Q)
    R, t = kabsch(P, Q)
    rms = float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean()))
    moved = [{"num": r["num"], "aa3": r["aa3"], "atoms": {a: tuple(R @ np.array(x) + t) for a, x in r["atoms"].items()}} for r in res]
    return moved, rms, len(P)


def main():
    her2 = read_atoms(T / "her2_domIV.cif", "A", label_seq=False)
    h1 = [r for r in her2 if r["num"] <= 580]
    h2 = [r for r in her2 if r["num"] >= 591]
    out = {}
    # 1) hNbBCII10-FGLA, PDB 3EAK chain A (CDR residues present; contig skips them). FR1 1-22 FR2 39-53 FR3 62-99 FR4 118-127
    res = read_atoms(FR / "3EAK.cif", "A")
    rec = build_record("bcii10", res, [(1, 22), (39, 53), (62, 99), (118, 127)])
    rec["source"] = "RCSB 3EAK chain A 'NbBCII10 humanized (FGLA mutant)', 1.95 A"
    # 2) caplacizumab VHH framework from the RFD4 VM example (CDRs are single placeholder residues, dropped here)
    res2 = read_atoms(FR / "pmvk_caplacizumab.cif", "H")
    rec2 = build_record("capla", res2, [(1, 22), (24, 37), (39, 75), (77, 86)],
                        stub_cdrs=("AASGRTFSSYAMG", "AISWSGGST", "AAGGRYSDYGSYDY"))
    rec2["source"] = "RFD4 example Data/03-vhh/pmvk_caplacizumab.cif (framework of caplacizumab; CDRs removed)"
    for name, res_, rec_ in (("bcii10", res, rec), ("capla", res2, rec2)):
        moved, rms, n = place_on_trastuzumab(res_, rec_)
        write_cif(INP / f"{name}_on_her2.cif", [("A", h1), ("C", h2), ("B", moved)], f"{name}_on_her2")
        write_cif(FR / f"{name}_clean.cif", [("B", res_)], f"{name}_clean")
        rec_["placement_rmsd_A"] = round(rms, 2); rec_["placement_n_ca"] = n
        (FR / f"{name}.json").write_text(json.dumps(rec_, indent=1))
        out[name] = {"fr": rec_["fr"], "native_cdr": rec_["native_cdr"], "contig_fr": rec_["contig_fr"],
                     "placement_rmsd_A": rec_["placement_rmsd_A"], "n_ca": n, "imgt_conf": round(rec_["imgt_confidence"], 2)}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
