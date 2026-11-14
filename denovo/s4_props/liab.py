"""Cheap sequence-liability / physico-chemical profile of VHHs (no structure, no GPU, ~1 s for 100 sequences).

Run in an env that has anarci + biopython (micromamba env `tnp` locally, `s4` on the VM):
    python liab.py --in seqs.csv --out liab.csv
seqs.csv: design_id, vhh_seq  (cdr1,cdr2,cdr3 optional; used only if ANARCI fails to number the sequence)

CDRs are IMGT (27-38, 56-65, 105-117), the same definition TNP uses.  Motif lists are standard developability
heuristics (not a trained model).  The numeric thresholds behind the `liab_f_*` flags (README.md) are set just above
the range of the well-behaved reference nanobodies in baselines.csv, so they flag outliers only; a flag is a prompt
to look, not a rejection (caplacizumab itself carries a DG in CDR3).
"""
from __future__ import annotations

import argparse
import csv
import re

KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
STRONG_HYD = set("VILMFW")
IMGT_CDR = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}


def imgt_cdrs(seq: str):
    """(cdr1, cdr2, cdr3, canonical_cys_ok) by ANARCI IMGT numbering, or None if numbering fails."""
    from anarci import run_anarci
    res = run_anarci([("q", seq)], scheme="imgt", allowed_species=None)
    num = res[1][0]
    if not num:
        return None
    numbering = num[0][0]
    cdr = {k: "" for k in IMGT_CDR}
    pos_aa = {}
    for (pos, _ins), aa in numbering:
        if aa == "-":
            continue
        pos_aa.setdefault(pos, "")
        pos_aa[pos] += aa
        for k, (lo, hi) in IMGT_CDR.items():
            if lo <= pos <= hi:
                cdr[k] += aa
    cys_ok = pos_aa.get(23, "") == "C" and pos_aa.get(104, "") == "C"
    return cdr["cdr1"], cdr["cdr2"], cdr["cdr3"], cys_ok


def max_window(s: str, w: int = 5) -> float:
    if not s:
        return 0.0
    if len(s) <= w:
        return sum(KD[a] for a in s) / len(s)
    return max(sum(KD[a] for a in s[i:i + w]) / w for i in range(len(s) - w + 1))


def longest_run(s: str, alphabet: set) -> int:
    best = cur = 0
    for a in s:
        cur = cur + 1 if a in alphabet else 0
        best = max(best, cur)
    return best


def profile(did: str, seq: str, given=None) -> dict:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    out = {"design_id": did}
    try:
        nm = imgt_cdrs(seq)
    except Exception:  # noqa: BLE001
        nm = None
    if nm:
        c1, c2, c3, cys_ok = nm
        out["liab_numbering"] = "anarci-imgt"
    else:
        c1, c2, c3 = (given or ("", "", ""))
        cys_ok = False
        out["liab_numbering"] = "given" if any((c1, c2, c3)) else "failed"
    cdrs = [c1, c2, c3]
    allc = "".join(cdrs)
    pa = ProteinAnalysis(seq)
    out["liab_len"] = len(seq)
    out["liab_cdr1"], out["liab_cdr2"], out["liab_cdr3"] = c1, c2, c3
    out["liab_n_cys"] = seq.count("C")
    out["liab_canon_cys"] = int(cys_ok)
    out["liab_unpaired_cys"] = int(seq.count("C") % 2 == 1)  # abguard's hard rule: an odd Cys count cannot all be paired
    out["liab_no_canon_ss"] = int(not cys_ok)  # canonical IMGT 23-104 disulfide absent: stability risk, not an abguard failure
    out["liab_nglyc"] = len(re.findall(r"(?=N[^P][ST])", seq))
    out["liab_nglyc_cdr"] = sum(len(re.findall(r"(?=N[^P][ST])", c)) for c in cdrs)
    out["liab_ng"] = len(re.findall(r"(?=NG)", seq))
    out["liab_ng_cdr"] = sum(c.count("NG") for c in cdrs)
    out["liab_ns_cdr"] = sum(c.count("NS") for c in cdrs)
    out["liab_dg"] = len(re.findall(r"(?=DG)", seq))
    out["liab_dg_cdr"] = sum(c.count("DG") for c in cdrs)
    out["liab_dp_cdr"] = sum(c.count("DP") for c in cdrs)
    out["liab_met_cdr"] = allc.count("M")
    out["liab_trp_cdr"] = allc.count("W")
    out["liab_pI"] = round(pa.isoelectric_point(), 2)
    out["liab_charge74"] = round(pa.charge_at_pH(7.4), 2)
    out["liab_cdr_pos"] = sum(allc.count(a) for a in "KR")
    out["liab_cdr_neg"] = sum(allc.count(a) for a in "DE")
    out["liab_cdr_gravy"] = round(sum(KD[a] for a in allc) / len(allc), 3) if allc else ""
    out["liab_cdr3_gravy"] = round(sum(KD[a] for a in c3) / len(c3), 3) if c3 else ""
    out["liab_max_hydro_win5"] = round(max(max_window(c) for c in cdrs), 3)
    out["liab_cdr_hyd_run"] = max(longest_run(c, STRONG_HYD) for c in cdrs)
    flags = {
        "unpaired_cys": out["liab_unpaired_cys"],
        "no_canon_ss": out["liab_no_canon_ss"],
        "nglyc": int(out["liab_nglyc"] > 0),
        "ng_cdr": int(out["liab_ng_cdr"] > 0),
        "dg_cdr": int(out["liab_dg_cdr"] > 0),
        "met_trp_cdr": int(out["liab_met_cdr"] + out["liab_trp_cdr"] >= 3),
        "hydrophobic": int(out["liab_max_hydro_win5"] >= 2.5 or out["liab_cdr_hyd_run"] >= 5),
        "pI": int(out["liab_pI"] < 5.5 or out["liab_pI"] > 9.5),
        "cdr_pos": int(out["liab_cdr_pos"] >= 6),
    }
    for k, v in flags.items():
        out[f"liab_f_{k}"] = v
    out["liab_n_flags"] = sum(flags.values())
    out["liab_flags"] = ";".join(k for k, v in flags.items() if v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.inp)))
    res = [profile(r["design_id"], r["vhh_seq"], (r.get("cdr1", ""), r.get("cdr2", ""), r.get("cdr3", ""))) for r in rows]
    cols = list(res[0].keys()) if res else ["design_id"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(res)


if __name__ == "__main__":
    main()
