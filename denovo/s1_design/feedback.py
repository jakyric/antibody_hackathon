"""Feedback from S2/S3 into the rfd4 campaign: which (framework, hotspot set, CDR3 range, pose) cells deserve jobs.

Success = the design passes the VM abguard de novo class (`denovo_vm_pass` in results/denovo/rules.csv, S2).
Among successes, quality = mean percentile (among OUR scored designs) of OF3 iptm, -min_ipae and trastuzumab-epitope overlap. S3's own
`consensus_rank_score` is NOT used: no rfd4 design clears S3's decoy-95th-percentile gates, so it is a constant (0.464) for 93% of designs.
S3's calibration gives these OF3 metrics only AUC ~0.65 on real VHH pairs, so the quality factor is clamped to 0.5-2 and shrunk hard.
  cell weight = shrunk VM pass rate  x  exp((shrunk mean OF3 quality among passes - overall mean) / (sd / 2))   [second factor only once
  >= MIN_SCORED passing designs are OF3-scored, clamped to 0.5-2], floor 0.25; 0 for a clear loser (n >= 10 scored passes and
  shrunk mean + 1.5 SE below overall mean - 0.5 sd). New cells with no data borrow the framework's pass rate among designs whose
  CDR3 block is >= 17 aa.
CLI: python feedback.py   (prints the table and writes denovo/s1_design/feedback.csv)
"""
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import dn_common as dn  # noqa: E402

MIN_SCORED = 40
CELL_KEYS = ("framework", "hotspot_set", "c3_range", "pose")


def _rows(p):
    return dn.read_csv(p) if Path(p).exists() else []


def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _truthy(x):
    return str(x).strip().lower() in ("1", "true", "1.0")


def quality(sc, ids):
    """design_id -> mean percentile over (of3_iptm high, of3_min_ipae low, of3_epi_overlap high), among the given scored designs."""
    spec = (("of3_iptm", 1), ("of3_min_ipae", -1), ("of3_epi_overlap", 1))
    pct = defaultdict(list)
    for m, sign in spec:
        xs = [(d, sign * _f(sc[d].get(m))) for d in ids if d in sc and _f(sc[d].get(m)) is not None]
        xs.sort(key=lambda t: t[1])
        n = len(xs)
        for r, (d, _) in enumerate(xs):
            pct[d].append(r / (n - 1) if n > 1 else 0.5)
    return {d: sum(v) / len(v) for d, v in pct.items() if len(v) >= 2}


def load():
    geom = {r["design_id"]: r for r in _rows(dn.OUT / "geom_rfd4jf.csv") if r.get("ok") == "1"}
    rules = {r["design_id"]: _truthy(r.get("denovo_vm_pass")) for r in _rows(dn.OUT / "rules.csv") if r.get("design_id")}
    sc = {}
    for r in _rows(dn.OUT / "score_orth.csv"):
        sc.setdefault(r["design_id"], {}).update(r)
    for r in _rows(dn.OUT / "score_of3.csv"):
        sc.setdefault(r["design_id"], {}).update({k: v for k, v in r.items() if v != ""})
    ours = [r["design_id"] for r in _rows(dn.OUT / "designs_rfd4jf.csv")]
    q = quality(sc, ours)
    for d, v in q.items():
        sc[d]["_q"] = v
    cdr3len = {r["design_id"]: len(r["cdr3"]) for r in _rows(dn.OUT / "designs_rfd4jf.csv")}
    fwof = {r["design_id"]: r["framework"] for r in _rows(dn.OUT / "designs_rfd4jf.csv")}
    return geom, rules, sc, cdr3len, fwof


def cell_table():
    geom, rules, sc, cdr3len, fwof = load()
    cells = defaultdict(lambda: {"n": 0, "rules_n": 0, "passes": 0, "scored_pass": [], "scored_all": [], "ipsae": [], "iptm": [],
                                 "epi": [], "epi_prec": [], "clash": []})
    for did, g in geom.items():
        c = cells[tuple(g[k] for k in CELL_KEYS)]
        c["n"] += 1
        c["epi_prec"].append(_f(g.get("epi_precision"))); c["clash"].append(_f(g.get("clash_lt2p2")))
        if did in rules:
            c["rules_n"] += 1
            c["passes"] += int(rules[did])
        s = sc.get(did)
        if s:
            v = s.get("_q")
            if v is not None:
                c["scored_all"].append(v)
                if rules.get(did):
                    c["scored_pass"].append(v)
            c["ipsae"].append(_f(s.get("of3_ipsae"))); c["iptm"].append(_f(s.get("of3_iptm"))); c["epi"].append(_f(s.get("of3_epi_overlap")))
    return cells


def fw_prior():
    """VM pass rate per framework among designs whose CDR3 block (incl. AA) is in the current range (capla 12-18, bcii10 15-18); else overall."""
    geom, rules, sc, cdr3len, fwof = load()
    num, den = defaultdict(int), defaultdict(int)
    for did, ok in rules.items():
        lo, hi = {"capla": (12, 18), "bcii10": (15, 18)}.get(fwof.get(did), (0, 99))
        if did in fwof and lo <= cdr3len.get(did, 0) <= hi:
            num[fwof[did]] += int(ok); den[fwof[did]] += 1
    allok = [ok for did, ok in rules.items() if did in fwof]
    overall = sum(allok) / len(allok) if allok else 0.5
    return {fw: (num[fw] / den[fw] if den[fw] >= 8 else overall) for fw in ("bcii10", "capla")}, overall


def plan(n_jobs, all_cells):
    cells = cell_table()
    prior, overall = fw_prior()
    allpass = [x for c in cells.values() for x in c["scored_pass"]]
    mu = sum(allpass) / len(allpass) if allpass else 0.0
    sd = ((sum((x - mu) ** 2 for x in allpass) / len(allpass)) ** 0.5 or 1.0) if allpass else 1.0
    use_q = len(allpass) >= MIN_SCORED
    weights, why = {}, {}
    K, K2 = 8.0, 4.0
    for cell in all_cells:
        c = cells.get(cell)
        p0 = prior[cell[0]]
        npass, nr = (c["passes"], c["rules_n"]) if c else (0, 0)
        pr = (npass + K * p0) / (nr + K)
        w, note = pr, f"VM pass {npass}/{nr} (prior {p0:.2f}) -> {pr:.2f}"
        if use_q:
            ns = len(c["scored_pass"]) if c else 0
            m = ((sum(c["scored_pass"]) if c else 0.0) + K2 * mu) / (ns + K2)
            se = sd / math.sqrt(ns + K2)
            if ns >= 10 and m + 1.5 * se < mu - 0.5 * sd:
                weights[cell], why[cell] = 0.0, note + f"; clear loser on OF3 consensus ({m:.3f} vs {mu:.3f}, n={ns})"
                continue
            mult = min(2.0, max(0.5, math.exp((m - mu) / (sd / 2))))
            w *= mult
            note += f"; OF3 quality among passes {m:.3f} vs {mu:.3f} (n={ns}) x{mult:.2f}"
        weights[cell], why[cell] = max(0.25 * overall, w) if w > 0 else 0.0, note
    tot = sum(weights.values()) or 1.0
    alloc = {cell: int(round(n_jobs * w / tot)) for cell, w in weights.items()}
    basis = (f"VM pass rate (prior by framework {', '.join(f'{k} {v:.2f}' for k, v in prior.items())}) x "
             + (f"OF3 quality among passes ({len(allpass)} scored passing designs)" if use_q
                else f"OF3 quality not yet used ({len(allpass)} scored passing designs, need {MIN_SCORED})"))
    return alloc, weights, why, basis


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def markdown():
    """per-cell markdown table for decisions.md: geometry, VM pass, OF3 quality, plus the wave-2 plan weights."""
    cells = cell_table()
    print("| framework | hotspots | CDR3 | pose | designs | VM de novo pass | OF3 iptm (mean) | OF3 epitope overlap | quality among passes (n) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for key, c in sorted(cells.items()):
        pr = f"{c['passes']}/{c['rules_n']}" if c["rules_n"] else "-"
        fm = lambda x: "-" if x is None else f"{x:.2f}"
        print(f"| {' | '.join(key)} | {c['n']} | {pr} | {fm(_mean(c['iptm']))} | {fm(_mean(c['epi']))} | {fm(_mean(c['scored_pass']))} ({len(c['scored_pass'])}) |")


if __name__ == "__main__":
    if "--md" in sys.argv:
        markdown(); sys.exit(0)
    cells = cell_table()
    out = []
    for key, c in sorted(cells.items()):
        out.append({**dict(zip(CELL_KEYS, key)), "n_designs": c["n"], "vm_pass_n": c["rules_n"], "vm_pass_rate": (c["passes"] / c["rules_n"] if c["rules_n"] else None),
                    "n_scored_pass": len(c["scored_pass"]), "mean_quality_pass": _mean(c["scored_pass"]),
                    "mean_quality_all": _mean(c["scored_all"]), "mean_of3_ipsae": _mean(c["ipsae"]), "mean_of3_iptm": _mean(c["iptm"]),
                    "mean_epi_overlap": _mean(c["epi"]), "geom_epi_precision": _mean(c["epi_prec"]), "geom_mean_clash_pairs": _mean(c["clash"])})
    if out:
        dn.write_csv(HERE / "feedback.csv", out)
    for r in out:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
