"""What drives the abguard identity-to-trastuzumab gate for VHHs? (data for gate_drivers.md)

Replicates abguard's formula (novelty.py `_identity_to_trastuzumab`):  identity = 1 - Lev(Fv_design, Fv_trastuzumab) / max(len, len)
and decomposes it with the Levenshtein alignment into matched / substituted / inserted / deleted residues, split framework vs CDR.
Tests: (1) CDR content randomised at fixed lengths, (2) CDR3 length sweep on the BCII10 framework with random CDR3 content.
Run: .venv/bin/python gate_drivers_analysis.py
"""
import csv, json, random, statistics as st, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_lib as L  # noqa
dn = L.dn
Lev = L.Levenshtein
REF = L.TRAS_FV
LREF = len(REF)
AA = "ACDEFGHIKLMNPQRSTVWY"
rng = random.Random(7)


def regions(seq):
    a = L.analyse(seq)
    reg = [("cdr" if L.imgt_region(l).startswith("cdr") else "fr") for l in a["idx_label"]]
    return a, reg


def ident(fv):
    return 1 - Lev.distance(fv, REF) / max(len(fv), LREF)


def decompose(fv, reg):
    eq_fr = eq_cdr = sub = dele = ins = 0
    for op in Lev.opcodes(fv, REF):
        n_src, n_dst = op.src_end - op.src_start, op.dest_end - op.dest_start
        if op.tag == "equal":
            for i in range(op.src_start, op.src_end):
                if reg[i] == "cdr": eq_cdr += 1
                else: eq_fr += 1
        elif op.tag == "replace":
            sub += n_src
        elif op.tag == "delete":   # residues of the design with no partner in trastuzumab (design longer here)
            dele += n_src
        elif op.tag == "insert":   # trastuzumab residues the design lacks (design shorter here)
            ins += n_dst
    return eq_fr, eq_cdr, sub, dele, ins


def swap_cdrs(a, reg, fn):
    """Fv with every CDR residue replaced by fn() (framework untouched, lengths untouched)."""
    return "".join(fn() if r == "cdr" else c for c, r in zip(a["fv"], reg))


def main():
    rules = {r["design_id"]: r for r in dn.read_csv(dn.OUT / "rules.csv")}
    rows = []
    for d in dn.load_designs(include_controls=False):
        if d.get("vhh_seq") is None or d["design_id"] not in rules:
            continue
        a, reg = regions(d["vhh_seq"])
        if not a.get("fv"):
            continue
        fv = a["fv"]
        eqf, eqc, sub, dele, ins = decompose(fv, reg)
        Lf = len(fv)
        rnd = [ident(swap_cdrs(a, reg, lambda: rng.choice(AA))) for _ in range(20)]
        floor = ident(swap_cdrs(a, reg, lambda: "X"))  # CDR residues can never match: lower bound of the CDR-content effect
        rr = rules[d["design_id"]]
        vm = rr.get("abguard_identity_trast") if rr.get("vhh_seq") == d["vhh_seq"] else ""  # only compare identical sequences
        rows.append(dict(design_id=d["design_id"], generator=d.get("generator", ""), framework=d.get("framework", ""), L=Lf, d=Lev.distance(fv, REF), identity=round(ident(fv), 4),
                         vm_identity=vm, eq_fr=eqf, eq_cdr=eqc, sub=sub, dele=dele, ins=ins, cdr1=len(a["cdr1"]), cdr2=len(a["cdr2"]),
                         cdr3=len(a["cdr3"]), n_fr=sum(r == "fr" for r in reg), identity_cdr_random=round(float(np.mean(rnd)), 4),
                         identity_cdr_random_sd=round(float(np.std(rnd)), 4), identity_cdr_nomatch=round(floor, 4),
                         seq_len=len(d["vhh_seq"])))
    dn.write_csv(HERE / "gate_drivers_data.csv", rows)
    out = {"n": len(rows), "ref_len": LREF, "ref_cdr_lengths": [len(L.analyse(L.TRAS_VH)[k]) for k in ("cdr1", "cdr2", "cdr3")]}
    chk = [abs(float(r["identity"]) - float(r["vm_identity"])) for r in rows if r["vm_identity"] not in (None, "")]
    out["max_abs_diff_vs_vm_abguard"] = max(chk) if chk else None
    groups = sorted({(r["generator"], r["framework"]) for r in rows})
    for gen, fw in groups:
        gname = f"{gen}/{fw}"
        rs = [r for r in rows if r["generator"] == gen and r["framework"] == fw]
        if len(rs) < 4:
            continue
        x = np.array([r["L"] for r in rs], float); y = np.array([r["identity"] for r in rs], float)
        m = np.array([r["L"] - r["d"] for r in rs], float)
        c3 = np.array([r["cdr3"] for r in rs], float)
        out[gname] = dict(n=len(rs), r_identity_vs_L=float(np.corrcoef(x, y)[0, 1]), r2_identity_vs_L=float(np.corrcoef(x, y)[0, 1] ** 2),
                          r_identity_vs_cdr3=float(np.corrcoef(c3, y)[0, 1]), slope_identity_per_residue=float(np.polyfit(x, y, 1)[0]),
                          m_mean=float(m.mean()), m_sd=float(m.std()), L_mean=float(x.mean()), L_sd=float(x.std()),
                          eq_fr_mean=float(np.mean([r["eq_fr"] for r in rs])), eq_cdr_mean=float(np.mean([r["eq_cdr"] for r in rs])),
                          n_fr_mean=float(np.mean([r["n_fr"] for r in rs])),
                          content_effect_real_minus_random=float(np.mean([r["identity"] - r["identity_cdr_random"] for r in rs])),
                          content_effect_real_minus_nomatch=float(np.mean([r["identity"] - r["identity_cdr_nomatch"] for r in rs])),
                          frac_denovo=float(np.mean(y < 0.70)))
        # pass rate by Fv length and by CDR3 length
        def by(key, edges):
            res = {}
            for lo, hi in edges:
                sel = [r for r in rs if lo <= r[key] <= hi]
                if sel:
                    res[f"{lo}-{hi}"] = dict(n=len(sel), de_novo=sum(r["identity"] < 0.70 for r in sel), median_identity=round(st.median(r["identity"] for r in sel), 3))
            return res
        out[gname]["by_Fv_length"] = by("L", [(0, 118), (119, 120), (121, 122), (123, 124), (125, 126), (127, 128), (129, 140)])
        out[gname]["by_cdr3_length"] = by("cdr3", [(0, 13), (14, 15), (16, 17), (18, 19), (20, 40)])
    # counterfactual: BCII10 framework + its own CDR1/CDR2, CDR3 length swept with random content
    bc = "QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAVRGYFMRLPSSHNFRYWGQGTLVTVSS"
    s = L.ANN.segment(bc)
    sweep = {}
    for l3 in range(8, 25):
        v = []
        for _ in range(30):
            cdr3 = "".join(rng.choice(AA) for _ in range(l3))
            fv = s.fr1 + s.cdr1 + s.fr2 + s.cdr2 + s.fr3 + cdr3 + s.fr4
            v.append(ident(fv))
        sweep[l3] = dict(Fv_length=len(s.fr1 + s.cdr1 + s.fr2 + s.cdr2 + s.fr3 + s.fr4) + l3, mean=round(float(np.mean(v)), 4), sd=round(float(np.std(v)), 4))
    out["bcii10_cdr3_length_sweep_random_content"] = sweep
    # same sweep but CDR3 content = trastuzumab's own CDR3 prefix/padding (best case for content): shows content can move identity by at most the matched residues
    cons = {}
    for l3 in (10, 13, 16, 18, 22):
        t = L.TRAS_CDR3
        cdr3 = (t + "A" * l3)[:l3]
        fv = s.fr1 + s.cdr1 + s.fr2 + s.cdr2 + s.fr3 + cdr3 + s.fr4
        cons[l3] = round(ident(fv), 4)
    out["bcii10_cdr3_trastuzumab_content_best_case"] = cons
    out["bcii10_wt_identity"] = round(ident(L.analyse(bc)["fv"]), 4)
    (HERE / "gate_drivers_stats.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:6000])


if __name__ == "__main__":
    main()
