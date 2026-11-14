"""The whole pipeline: proposals → gate → filtering → ranking → designs.csv

Run:  make designs
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .data import TZ_CDRH3_10
from .gate import ProvidedIndex
from .generate import Generator
from .model import load as load_model
from .rank import diversify, stage1_filter, stage2_rank
from .sequences import TRASTUZUMAB_VL

RESULTS = Path(__file__).resolve().parents[2] / "results"
OUT = RESULTS / "designs.csv"


def run(n_candidates: int = 20_000, n_submit: int = 20, seed: int = 0,
        keep_percentile: float = 70.0, verbose: bool = True,
        out: Path | None = None) -> list:
    out = Path(out) if out else OUT
    def say(*a):
        if verbose:
            print(*a)

    say("=" * 70)
    say("DESIGN RUN — trastuzumab, optimization track")
    say("=" * 70)

    say(f"\n1. Proposing {n_candidates:,} variants with 4–5 mutations …")
    gen = Generator(seed=seed)
    cands = gen.sample(n_candidates)
    say(f"   {len(cands):,} unique")

    say("\n2. The gate: novelty against all provided sequences …")
    provided = ProvidedIndex()
    say(f"   index over {provided.n:,} provided CDR-H3s")
    legal = [c for c in cands if provided.is_novel(c)]
    say(f"   {len(legal):,} passed ({100*len(legal)/len(cands):.1f} %), "
        f"{len(cands)-len(legal):,} were too close to something provided")

    if not legal:
        raise SystemExit("no legal candidates — adjust the generator")

    say(f"\n3. Filtering out those the network thinks are unlikely to bind "
        f"(keeping the top {100-keep_percentile:.0f} %) …")
    net = load_model()
    kept, z_all, pct_all = stage1_filter(legal, net, keep_percentile)
    idx = {c: i for i, c in enumerate(legal)}
    z = np.array([z_all[idx[c]] for c in kept])
    pct = np.array([pct_all[idx[c]] for c in kept])
    say(f"   {len(kept):,} remaining")

    say("\n4. Ranking by what actually carried over to SPR …")
    ranked = stage2_rank(kept, z, pct, provided)
    ok = [d for d in ranked if d.gate_ok]
    say(f"   {len(ok):,} passed the whole gate")

    say(f"\n5. Choosing {n_submit} that are not variants of each other …")
    chosen = diversify(ok, n_submit)
    for i, d in enumerate(chosen, 1):
        d.design_id = f"tz-opt-{i:02d}"
    say(f"   {len(chosen)} chosen")

    RESULTS.mkdir(parents=True, exist_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "design_id", "cdrh3", "n_mutations_from_trastuzumab",
                    "min_mutations_from_provided", "model_percentile",
                    "conservative_substitution_fraction", "format",
                    "vh_sequence", "vl_sequence", "scfv_sequence",
                    "total_length", "rationale"])
        for i, d in enumerate(chosen, 1):
            w.writerow([i, d.design_id, d.cdrh3, d.n_mut, d.min_dist_provided,
                        f"{d.model_percentile:.1f}", f"{d.conservative:.2f}",
                        "scFv", d.vh, TRASTUZUMAB_VL, d.scfv,
                        len(d.scfv), d.rationale()])

    say(f"\nWritten: {out}")
    say("\nNEXT STEPS, and they are not optional:")
    say("  · run the official checker. It is the ground truth, not this gate.")
    say("  · co-fold against HER2 and look at the interface, not just the sequence.")
    say("  · thermostability lives in the framework, which the novelty filter does not")
    say("    touch. It is the cheapest of the three axes and has not been done here.")
    return chosen


def main() -> None:
    p = argparse.ArgumentParser(description="Build a ranked design set")
    p.add_argument("--candidates", type=int, default=20_000)
    p.add_argument("--submit", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--keep-percentile", type=float, default=70.0)
    p.add_argument("--out", type=Path, default=None,
                   help="where the designs are written (default: results/designs.csv)")
    a = p.parse_args()
    run(a.candidates, a.submit, a.seed, a.keep_percentile, out=a.out)


if __name__ == "__main__":
    main()
