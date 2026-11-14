"""Calibrate every scorer against controls and decoys, choose per-model gates, report per generator.

  denovo/.venv/bin/python denovo/s3_score/calibrate.py            # writes calib_scores.csv, gates.json, calibration.md

Families (all scored with exactly the same pipeline as the designs):
  fv-spr   trastuzumab Fv variants vs HER2: measured SPR binders (KD <= 10 nM) vs 'no reportable KD'
  fv-facs  trastuzumab Fv variants vs HER2: FACS high vs low
  vhh-real crystal-structure VHH-antigen pairs (cognate) vs shuffled pairs and HARD shuffled pairs (AF3 ipTM >= 0.5)
  vhh-hard cognate vs the hard shuffled pairs only
  decoys   every design vs hen egg-white lysozyme, and framework-only VHHs (controls ctl-fwrand-*) vs HER2
A gate is a threshold at the 95th percentile of the decoy distribution of the model's primary metric. We keep the
metric with the best AUC on vhh-real (the family closest to the de novo task), ties broken towards ipSAE.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "denovo"))
sys.path.insert(0, str(HERE))
import dn_common as dn  # noqa: E402
import run_s3 as R  # noqa: E402

METRICS = {  # metric -> direction (high = larger is better)
    "iptm": "high", "ipsae": "high", "ipsae15": "high", "ipsae20": "high", "ipae_tm": "high", "confidence": "high", "complex_plddt": "high", "binder_plddt": "high",
    "min_ipae": "low", "mean_ipae": "low", "ipae_contacts": "low", "interface_pae": "low",
    "ipsae_chn10": "high", "ipsae_chn15": "high", "pose_rmsd_vhh": "low",
}
PRIMARY = {"of3": "ipsae_chn15", "dde": "min_ipae", "esm2": "ipsae_chn15", "af2ig": "interface_pae", "of3s1": "ipsae_chn15"}
FAMILIES = {
    "fv-spr": (["spr-bind"], ["spr-non"]),
    "fv-facs": (["facs-high"], ["facs-low"]),
    "vhh-real": (["vhhpos"], ["vhhneg", "vhhhard"]),
    "vhh-hard": (["vhhpos"], ["vhhhard"]),
}


def auc(pos, neg, direction="high"):
    pos, neg = np.array(pos, float), np.array(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) < 3 or len(neg) < 3:
        return float("nan"), float("nan"), float("nan")
    u = mannwhitneyu(pos, neg, alternative="two-sided").statistic / (len(pos) * len(neg))
    a = u if direction == "high" else 1 - u
    rng = np.random.RandomState(0)
    bs = []
    for _ in range(300):
        p, n = rng.choice(pos, len(pos)), rng.choice(neg, len(neg))
        uu = mannwhitneyu(p, n, alternative="two-sided").statistic / (len(p) * len(n))
        bs.append(uu if direction == "high" else 1 - uu)
    return float(a), float(np.percentile(bs, 5)), float(np.percentile(bs, 95))


def load_long():
    st = R.load_state()
    calib = {f"{r['calib_id']}__{r['antigen_name']}": r for r in dn.read_csv(R.CALIB)} if R.CALIB.exists() else {}
    rows = []
    for (key, mk), r in st.items():
        if r["status"] != "ok" or key not in calib:
            continue
        c = calib[key]
        rows.append({"key": key, "group": c["group"], "label": int(c["label"]), "model": mk, **{k: v for k, v in r["metrics"].items()
                                                                                         if isinstance(v, (int, float))}})
    return rows


def decoy_rows():
    """Per-model metric values for the decoy classes of the design scoring runs."""
    st = R.load_state()
    designs = {d["design_id"]: d for d in R.load_designs()}
    out = {"hewl": [], "fwrand": []}
    for (key, mk), r in st.items():
        if r["status"] != "ok" or "__" not in key:
            continue
        did, ant = key.rsplit("__", 1)
        if did not in designs:
            continue
        gen = designs[did]["generator"]
        if ant == "hewl" and gen != "control":
            out["hewl"].append({"model": mk, **r["metrics"]})
        if ant == "her2" and did.startswith("ctl-fwrand"):
            out["fwrand"].append({"model": mk, **r["metrics"]})
    return out


def main():
    st = R.load_state()
    rows = load_long()
    if not rows:
        sys.exit("no calibration results yet")
    dn.write_csv(HERE / "calib_scores.csv", rows)
    dec = decoy_rows()
    models = sorted({r["model"] for r in rows})
    md = ["# Calibration of the scorers\n",
          "AUC = P(positive scores better than negative), 90% bootstrap interval in brackets. n = positives/negatives.\n"]
    gates = {"models": {}, "note": "thr = 95th percentile of decoys (vhh negatives + lysozyme decoys + framework-only VHHs)"}
    for mk in models:
        md.append(f"\n## {mk} ({R.MODELS[mk][0]})\n")
        md.append("| family | metric | n | AUC [5%, 95%] |\n|---|---|---|---|")
        best = {}
        for fam, (pg, ng) in FAMILIES.items():
            P = [r for r in rows if r["model"] == mk and r["group"] in pg]
            N = [r for r in rows if r["model"] == mk and r["group"] in ng]
            for met, d in METRICS.items():
                a, lo, hi = auc([r.get(met, float("nan")) for r in P], [r.get(met, float("nan")) for r in N], d)
                if not math.isnan(a):
                    md.append(f"| {fam} | {met} | {len(P)}/{len(N)} | {a:.2f} [{lo:.2f}, {hi:.2f}] |")
                    if fam == "vhh-real":
                        best[met] = a
        if not best:
            continue
        # primary metric fixed a priori (13 candidate metrics on 10 positives would select noise); the best-AUC metric is only reported
        top = max(best.values())
        met = PRIMARY.get(mk) or max(best, key=best.get)
        md.append(f"\nPrimary metric for {mk}: **{met}** (chosen a priori; vhh-real AUC {best.get(met, float('nan')):.2f}; best of all metrics: "
                  f"{max(best, key=best.get)} {top:.2f}, i.e. differences are within noise).\n")
        d = METRICS[met]
        pool = [r for r in rows if r["model"] == mk and r["group"] in ("vhhneg", "vhhhard")]
        pool += [r for r in dec["hewl"] if r["model"] == mk] + [r for r in dec["fwrand"] if r["model"] == mk]
        vals = np.array([r.get(met, float("nan")) for r in pool], float)
        vals = vals[~np.isnan(vals)]
        if len(vals) < 8:
            continue
        thr = float(np.percentile(vals, 95 if d == "high" else 5))
        pos = np.array([r.get(met, float("nan")) for r in rows if r["model"] == mk and r["group"] == "vhhpos"], float)
        rec = float(np.mean(pos >= thr)) if d == "high" else float(np.mean(pos <= thr))
        gates["models"][mk] = {"metric": met, "direction": d, "thr": thr, "epi_min": 0.5, "n_decoys": int(len(vals)),
                               "recall_cognate_vhh": rec, "auc_vhh_real": best[met]}
        md.append(f"\n**Gate for {mk}**: {met} {'>=' if d == 'high' else '<='} {thr:.3f} (95th percentile of {len(vals)} decoys); "
                  f"passes {rec:.0%} of the {len(pos)} cognate crystal VHH pairs.\n")
    (HERE / "gates.json").write_text(json.dumps(gates, indent=1))
    # per-generator behaviour of designs
    calib_keys = {f"{r['calib_id']}__{r['antigen_name']}": r["group"] for r in dn.read_csv(R.CALIB)}
    designs = {d["design_id"]: d for d in R.load_designs()}
    per = {}
    for (key, mk), r in st.items():
        if r["status"] != "ok" or not key.endswith("__her2"):
            continue
        did = key[:-6]
        if did in designs and mk in gates["models"]:
            gate = gates["models"][mk]
            v = r["metrics"].get(gate["metric"], float("nan"))
            ok = (v >= gate["thr"] if gate["direction"] == "high" else v <= gate["thr"]) and r["metrics"].get("epi_overlap", 0) >= gate["epi_min"]
            per.setdefault((designs[did]["generator"], mk), []).append((v, ok, r["metrics"].get("epi_overlap", float("nan"))))
    if per:
        md.append("\n## Per generator\n\n| generator | model | n | median primary metric | share past gate (incl. epitope overlap) | median epitope overlap |\n|---|---|---|---|---|---|")
        for (gen, mk), v in sorted(per.items()):
            vs = np.array([x[0] for x in v], float)
            md.append(f"| {gen} | {mk} | {len(v)} | {np.nanmedian(vs):.3f} | {np.mean([x[1] for x in v]):.0%} | "
                      f"{np.nanmedian([x[2] for x in v]):.2f} |")
    # AF2-IG: designs (their own pose) vs controls and calibration decoys (OF3 pose)
    af = [(k, r) for (k, mk), r in st.items() if mk == "af2ig" and r["status"] == "ok"]
    if af:
        def grp(k):
            did = k.rsplit("__", 1)[0]
            if k in calib_keys:
                return "calib:" + calib_keys[k]
            if did in designs:
                return "control(OF3 pose)" if designs[did]["generator"] == "control" else f"{designs[did]['generator']}(design pose)"
            return None
        groups = {}
        for k, r in af:
            gname = grp(k)
            if gname:
                groups.setdefault(gname, []).append(r["metrics"])
        md.append("\n## AF2-IG (re-scores a given pose; deterministic, no MSA)\n\n| group | n | iptm median [max] | interface PAE median [min] | min_ipae median | ipsae_chn15 median | pose RMSD vs input, VHH median [A] |\n|---|---|---|---|---|---|---|")
        for gname, ms in sorted(groups.items()):
            f = lambda k: [m[k] for m in ms if isinstance(m.get(k), (int, float)) and not math.isnan(m[k])]
            md.append(f"| {gname} | {len(ms)} | {np.median(f('iptm')):.2f} [{max(f('iptm')):.2f}] | "
                      f"{(np.median(f('interface_pae')) if f('interface_pae') else float('nan')):.1f} [{(min(f('interface_pae')) if f('interface_pae') else float('nan')):.1f}] | "
                      f"{(np.median(f('min_ipae')) if f('min_ipae') else float('nan')):.1f} | {(np.median(f('ipsae_chn15')) if f('ipsae_chn15') else float('nan')):.2f} | "
                      f"{(np.median(f('pose_rmsd_vhh')) if f('pose_rmsd_vhh') else float('nan')):.1f} |")
        md.append("\nCaveat: designs are scored on the generator's own pose, controls and calibration decoys on an OpenFold3 pose, so a gap between the "
                  "groups mixes pose quality with scorer behaviour; the calibration rows (cognate vs shuffled, both on OF3 poses) isolate the scorer.\n")
    (HERE / "calibration.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
