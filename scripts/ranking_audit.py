"""Audit the ranking signal: what are we actually ordering our designs by?

THREE QUESTIONS

1. The report quotes both +0.16 and -0.04 for the FACS net against measured
   affinity in the legal zone. Only one run can be right. Recompute both the
   net and the baselines on the SPR table, with intervals, and state the sign
   convention explicitly so the number cannot be misread again.

2. What is the SIGN of edit distance against affinity? `baseline_edit_distance`
   returns -hamming, so a HIGHER score means CLOSER to trastuzumab. The target
   is -log10(KD), where HIGHER means tighter binding. A positive Spearman
   therefore reads: closer to trastuzumab binds tighter. Worth stating,
   because the double negation is easy to get backwards.

3. THE ONE THAT MATTERS. All 20 of our designs sit at exactly 4 CDRH3
   mutations from trastuzumab, and exactly 3 from the nearest provided
   sequence. In the ranking formula

       rank = -2.0*n_mut + 0.8*cons + 0.4*(pct/100) + 0.1*min(d_prov, 6)

   both distance terms are therefore CONSTANT across all twenty. The only
   signal that transferred to measured affinity contributes nothing to the
   order of our own list. The ordering comes entirely from
   `conservative_fraction` and the net's percentile -- and the net is the
   thing we measured does not transfer.

   So: restrict the SPR data to variants at exactly 4 mutations, the stratum
   our designs live in, and test every scorer there. If nothing works inside
   that stratum, our ranking is noise dressed as a number, and we have to say
   so out loud rather than ship it quietly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abx.benchmark import (  # noqa: E402
    MIN_LEGAL_MUTATIONS, baseline_edit_distance, baseline_random,
    evaluate_on_spr,
)
from abx.data import TZ_CDRH3_10, load_spr  # noqa: E402
from abx.gate import hamming  # noqa: E402
from abx.metrics import bootstrap_ci, spearman  # noqa: E402
from abx import model as model_mod  # noqa: E402
from abx.model import MLP, onehot  # noqa: E402
from abx.rank import conservative_fraction  # noqa: E402


def load_net() -> MLP | None:
    for p in (model_mod.MODEL_PATH, Path("results/binder_model.npz")):
        if Path(p).exists():
            return model_mod.load(Path(p))
    return None


def net_scorer(net: MLP):
    def score(seqs):
        return net.logits(onehot(list(seqs))).ravel()
    return score


def report(name: str, s: np.ndarray, target: np.ndarray) -> None:
    if len(s) < 3 or np.allclose(s, s[0]):
        print(f"  {name:<28} CONSTANT across this stratum - no ordering at all")
        return
    rho = spearman(s, target)
    lo, hi = bootstrap_ci(spearman, s, target)
    verdict = "clear of zero" if lo > 0 or hi < 0 else "includes zero"
    print(f"  {name:<28} rho {rho:+.3f}  CI [{lo:+.3f}, {hi:+.3f}]  {verdict}")


def main() -> int:
    spr = load_spr()
    keep = spr.has_kd
    seqs = [s for s, k in zip(spr.cdrh3, keep) if k]
    kd = np.asarray(spr.kd_nM)[keep]
    target = -np.log10(kd)          # higher = tighter binding
    muts = np.array([hamming(s, TZ_CDRH3_10) if len(s) == len(TZ_CDRH3_10)
                     else 99 for s in seqs])

    net = load_net()
    print(f"SPR variants with measured KD: {len(seqs)}")
    print("target = -log10(KD); HIGHER = TIGHTER BINDING")
    print("edit-distance scorer = -hamming; HIGHER = CLOSER TO TRASTUZUMAB")
    print("so a POSITIVE rho means: closer to trastuzumab binds tighter\n")

    print("distribution of CDRH3 mutations among measured variants:")
    for d in sorted(set(muts.tolist())):
        n = int((muts == d).sum())
        tag = "  <- our designs live here" if d == 4 else ""
        print(f"  {d:>2} mutation(s): {n:>4}{tag}")
    print()

    scorers = {
        "edit distance": baseline_edit_distance,
        "conservative fraction": lambda ss: np.array(
            [conservative_fraction(s) for s in ss]),
        "random": baseline_random(0),
    }
    if net is not None:
        scorers["FACS net"] = net_scorer(net)
    else:
        print("!! no trained net found; skipping the net rows\n")

    # ---- question 1 and 2: the headline numbers ---------------------------
    for legal in (False, True):
        tag = (f"LEGAL ZONE (d >= {MIN_LEGAL_MUTATIONS})" if legal
               else "ALL MEASURED VARIANTS")
        print(f"=== {tag} ===")
        for nm, sc in scorers.items():
            r = evaluate_on_spr(sc, nm, spr, legal_only=legal)
            lo, hi = r.ci
            flag = "" if (lo > 0 or hi < 0) else "   (includes zero)"
            print(f"  {nm:<28} rho {r.value:+.3f}  "
                  f"CI [{lo:+.3f}, {hi:+.3f}]  n={r.n}{flag}")
        print()

    # ---- question 3: inside the stratum our designs occupy ----------------
    print("=== THE STRATUM OUR 20 DESIGNS ACTUALLY OCCUPY ===")
    for d in (4, 5):
        sel = muts == d
        n = int(sel.sum())
        print(f"\nexactly {d} CDRH3 mutation(s), n = {n}")
        if n < 5:
            print("  too few measured variants to say anything")
            continue
        sub = [s for s, k in zip(seqs, sel) if k]
        t = target[sel]
        for nm, sc in scorers.items():
            report(nm, np.asarray(sc(sub), dtype=float), t)

    # combined 4-5 for a little more power
    sel = (muts >= 4) & (muts <= 5)
    if int(sel.sum()) >= 8:
        print(f"\n4-5 mutations pooled, n = {int(sel.sum())}")
        sub = [s for s, k in zip(seqs, sel) if k]
        t = target[sel]
        for nm, sc in scorers.items():
            report(nm, np.asarray(sc(sub), dtype=float), t)

    print("\n=== WHAT ORDERS OUR OWN 20 ===")
    dpath = Path("results/designs.csv")
    if dpath.exists():
        import csv
        with dpath.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        nmut = {int(r["n_mutations_from_trastuzumab"]) for r in rows}
        dprov = {int(r["min_mutations_from_provided"]) for r in rows}
        cons = np.array([float(r["conservative_substitution_fraction"])
                         for r in rows])
        pct = np.array([float(r["model_percentile"]) for r in rows])
        print(f"  n_mutations_from_trastuzumab: {sorted(nmut)}  "
              f"-> term -2.0*n_mut contributes {'NOTHING' if len(nmut)==1 else 'something'}")
        print(f"  min_mutations_from_provided:  {sorted(dprov)}  "
              f"-> term 0.1*d_prov contributes {'NOTHING' if len(dprov)==1 else 'something'}")
        print(f"  conservative fraction spread: {cons.min():.2f}-{cons.max():.2f}"
              f"  -> term 0.8*cons spans {0.8*(cons.max()-cons.min()):.3f}")
        print(f"  model percentile spread:      {pct.min():.0f}-{pct.max():.0f}"
              f"   -> term 0.4*pct spans {0.4*(pct.max()-pct.min())/100:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
