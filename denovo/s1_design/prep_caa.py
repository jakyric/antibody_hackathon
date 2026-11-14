"""Wave-2 inputs: keep the invariant C104-A105-A106 of CDR3 FIXED (monoyield paper: CDR3 packs on this 'CAA' and F42).

bcii10: the real A105/A106 (3EAK residues 100,101) are simply added to FR3 (FR3 = 62-101, designed block after them).
capla : the example structure has no CDR3 residues, so A105/A106 are GRAFTED: the backbone of bcii10's two residues is moved onto
        capla by superposing the 4 residues that end FR3 (the last is C104), and written as ALA 76 and 77; FR4 is renumbered +3.
Writes inputs/capla_caa_on_her2.cif and adds a "caa" record to frameworks/<fw>.json.   Run after prep_inputs.py.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import prep_inputs as PI  # noqa: E402

INP, FR = HERE / "inputs", HERE / "frameworks"


def bb(res, nums):
    return np.array([res[n]["atoms"][a] for n in nums for a in ("N", "CA", "C", "O")])


def main():
    out = {}
    # ---- bcii10
    rec = json.load(open(FR / "bcii10.json"))
    b = {r["num"]: r for r in PI.read_atoms(INP / "bcii10_on_her2.cif", "B", label_seq=False)}
    seq = lambda lo, hi: "".join(PI.AA3[b[n]["aa3"]] for n in range(lo, hi + 1))
    caa_fr = [seq(1, 22), seq(39, 53), seq(62, 101), seq(118, 127)]
    assert caa_fr[2].endswith("YYCAA"), caa_fr[2]
    rec["caa"] = {"fr": caa_fr, "contig_fr": ["B1-22", "B39-53", "B62-101", "B118-127"], "prefix": "AA", "grafted": False,
                  "input": "bcii10_on_her2.cif"}
    json.dump(rec, open(FR / "bcii10.json", "w"), indent=1)
    out["bcii10"] = caa_fr[2][-8:]
    # ---- capla (graft A105 A106 from the placed bcii10)
    rc = json.load(open(FR / "capla.json"))
    cres = PI.read_atoms(INP / "capla_on_her2.cif", "B", label_seq=False)
    c = {r["num"]: r for r in cres}
    P = bb(b, [96, 97, 98, 99])      # bcii10 ... Y Y C104 (in 1N8Z frame)
    Q = bb(c, [72, 73, 74, 75])      # capla  ... Y Y C104 (same frame)
    R, t = PI.kabsch(P, Q)
    rms = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    graft = []
    for n, newnum in ((100, 76), (101, 77)):
        atoms = {a: tuple(R @ np.array(x) + t) for a, x in b[n]["atoms"].items()}
        graft.append({"num": newnum, "aa3": b[n]["aa3"], "atoms": atoms})
    bond = float(np.linalg.norm(np.array(c[75]["atoms"]["C"]) - np.array(graft[0]["atoms"]["N"])))
    print(f"graft rmsd {rms:.2f} A over 16 backbone atoms; C104(C)-A105(N) distance {bond:.2f} A (ideal 1.33)")
    assert 1.1 < bond < 1.7, "graft does not connect to C104"
    keep = [r for r in cres if r["num"] <= 75]
    fr4 = [{"num": r["num"] + 3, "aa3": r["aa3"], "atoms": r["atoms"]} for r in cres if r["num"] >= 77]
    chain_b = keep + graft + fr4
    ha = PI.read_atoms(INP / "capla_on_her2.cif", "A", label_seq=False)
    hc = PI.read_atoms(INP / "capla_on_her2.cif", "C", label_seq=False)
    PI.write_cif(INP / "capla_caa_on_her2.cif", [("A", ha), ("C", hc), ("B", chain_b)], "capla_caa_on_her2")
    fr = rc["fr"]
    caa_fr_c = [fr[0], fr[1], fr[2] + "AA", fr[3]]
    rc["caa"] = {"fr": caa_fr_c, "contig_fr": ["B1-22", "B24-37", "B39-77", "B80-89"], "prefix": "AA", "grafted": True,
                 "graft_rmsd_A": round(rms, 2), "graft_bond_A": round(bond, 2), "input": "capla_caa_on_her2.cif"}
    json.dump(rc, open(FR / "capla.json", "w"), indent=1)
    out["capla"] = caa_fr_c[2][-8:]
    print(json.dumps(out))


if __name__ == "__main__":
    main()
