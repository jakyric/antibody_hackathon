"""Nine-rule (and six-rule) CDR2/CDR3 sequence score for nanobody monomeric yield.

Source: "Massively parallel quantification and prediction of monomeric protein yield in Escherichia coli" (Fig. 5h-i and
Methods, "Nanobody sequence rules"); the six-rule version is adapted there from McMahon et al. 2018.  Definitions as
implemented (checked against the Methods text):

  CDR2  the six residues at scaffold positions 55-60 = chars [1:7] of an 8-residue IMGT CDR2 (IMGT 56-65, e.g. IASMGGLT -> ASMGGL).
        Any other CDR2 length: the four CDR2 rules are NOT APPLICABLE (not counted as pass or fail).
  CDR3  the loop AFTER the framework Cys + 2 residues (C-A-A) and BEFORE FR4 W-G-x-G, i.e. IMGT 105-117 minus its first two residues.
        N1, N2 = first two residues; C1, C2, C3 = last three counted from the C-terminus (C1 = last).  Needs >= 5 residues
        to give five distinct positions, else the five CDR3 rules are not applicable.
  rule 1 CDR2 pos1 in A/D/G/N/S/T          rule 5 (broadened) C1 in Y/W     [six-rule version: C1 = Y]
  rule 2 CDR2 pos4 = G                     rule 6 C3 in F/H/L/Y
  rule 3 CDR2 pos5 in A/G/S/T              rule 7 N1 NOT in Y/W/V/F/L/I
  rule 4 CDR2 pos6 in I/N/S/T              rule 8 N2 in L/T
                                           rule 9 C2 NOT in F/L/Y
  six-rule score = rules 1-4 + (C1 = Y) + rule 6.   Cys and Met in the scored loops are recorded (`nb_has_cys_met`), not scored.

Outputs: nb9_score (passed rules, 0-9), nb9_applicable (rules that could be evaluated), nb6_score, nb6_applicable,
nb9_flags (failed rules, ';'-separated), nb_has_cys_met, nb_cdr2_core, nb_cdr3_loop.
Caveats (from the paper): the three added rules and the C1 broadening were derived on the same dataset and are exploratory
(33% of nine-rule designs > ~3 mg/L vs 29% for the six-rule version); data are from ONE fixed-framework scaffold with an
invariant CDR1, a six-residue CDR2 and 5-13 residue CDR3s, so longer CDR3s and other frameworks are extrapolation.
Use as a small yield prior, never as a gate.

    python nanobody_rules.py --in seqs.csv --out nb9.csv      (env with anarci, e.g. micromamba env `tnp`)
    python nanobody_rules.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

RULES = [  # (name, scheme membership, test on (cdr2core, cdr3loop) -> bool)
    ("r1_cdr2_p1", "c2", lambda c2, c3: c2[0] in "ADGNST"),
    ("r2_cdr2_p4", "c2", lambda c2, c3: c2[3] == "G"),
    ("r3_cdr2_p5", "c2", lambda c2, c3: c2[4] in "AGST"),
    ("r4_cdr2_p6", "c2", lambda c2, c3: c2[5] in "INST"),
    ("r5_c1", "c3", lambda c2, c3: c3[-1] in "YW"),
    ("r6_c3", "c3", lambda c2, c3: c3[-3] in "FHLY"),
    ("r7_n1", "c3", lambda c2, c3: c3[0] not in "YWVFLI"),
    ("r8_n2", "c3", lambda c2, c3: c3[1] in "LT"),
    ("r9_c2", "c3", lambda c2, c3: c3[-2] not in "FLY"),
]
SIX = ["r1_cdr2_p1", "r2_cdr2_p4", "r3_cdr2_p5", "r4_cdr2_p6", "r5_c1_strict", "r6_c3"]


def split_loops(seq: str, given_cdr2: str = "", given_cdr3: str = ""):
    """(cdr2 IMGT string, cdr3 loop after C-A-A) from ANARCI IMGT numbering, falling back to the given IMGT CDR strings."""
    cdr2 = cdr3 = ""
    try:
        from liab import imgt_cdrs
        nm = imgt_cdrs(seq)
    except Exception:  # noqa: BLE001
        nm = None
    if nm:
        _, cdr2, cdr3, _ = nm
    else:
        cdr2, cdr3 = given_cdr2 or "", given_cdr3 or ""
    return cdr2, cdr3[2:] if len(cdr3) > 2 else ""


def score(cdr2: str, cdr3loop: str) -> dict:
    core = cdr2[1:7] if len(cdr2) == 8 else ""
    loop = cdr3loop if len(cdr3loop) >= 5 else ""
    res = {}
    for name, grp, fn in RULES:
        if (grp == "c2" and not core) or (grp == "c3" and not loop):
            continue
        res[name] = bool(fn(core, loop))
    if loop:
        res["r5_c1_strict"] = loop[-1] == "Y"
    failed = [k for k, v in res.items() if not v and k != "r5_c1_strict"]
    nine = {k: v for k, v in res.items() if k != "r5_c1_strict"}
    six = {k: res[k] for k in SIX if k in res}
    return {
        "nb9_score": sum(nine.values()), "nb9_applicable": len(nine),
        "nb6_score": sum(six.values()), "nb6_applicable": len(six),
        "nb9_flags": ";".join(failed),
        "nb_has_cys_met": int(any(a in "CM" for a in core + loop)),
        "nb_cdr2_core": core, "nb_cdr3_loop": loop,
    }


def selftest():
    # parent hNbBCII10-FGLA: must give 6/9, failing N1=V, N2=R and CDR2 position 6 = L
    seq = ("QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAVRGYFMRLPSSHNFRYWGQGTLVTVSS")
    c2, c3 = split_loops(seq)
    r = score(c2, c3)
    assert (c2, c3) == ("IASMGGLT", "VRGYFMRLPSSHNFRY"), (c2, c3)
    assert r["nb9_score"] == 6 and r["nb9_applicable"] == 9, r
    assert set(r["nb9_flags"].split(";")) == {"r4_cdr2_p6", "r7_n1", "r8_n2"}, r
    assert r["nb6_score"] == 5 and r["nb6_applicable"] == 6, r  # fails only CDR2 pos 6 (C1 = Y and C3 = F pass)
    assert r["nb_cdr2_core"] == "ASMGGL" and r["nb_has_cys_met"] == 1  # M in the CDR3 loop (recorded, not scored)
    r2 = score("IASMGGL", "VRGYFMRLPSSHNFRY")  # 7-residue CDR2: CDR2 rules not applicable
    assert r2["nb9_applicable"] == 5, r2
    print("selftest ok:", r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp")
    ap.add_argument("--out")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    rows = list(csv.DictReader(open(a.inp)))
    out = []
    for r in rows:
        c2, c3 = split_loops(r["vhh_seq"], r.get("cdr2", ""), r.get("cdr3", ""))
        out.append({"design_id": r["design_id"], **score(c2, c3)})
    cols = ["design_id", "nb9_score", "nb9_applicable", "nb6_score", "nb6_applicable", "nb9_flags", "nb_has_cys_met",
            "nb_cdr2_core", "nb_cdr3_loop"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
