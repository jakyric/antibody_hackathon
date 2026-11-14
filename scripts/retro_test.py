#!/usr/bin/env python3
"""Hide the measured winners, re-derive the scorer, and see whether it finds them.

Every other test we have says how well a scorer ranks variants we have already
measured. This task is different: picking a few candidates out of a space we
have not measured. This test imitates the task.

    the pool     the 150 measured 10-mers in the legal region (>= 3 mutations)
    the winners  the 10 of them that are stronger than 1 nM
    the scorer   is trained on FACS without the winners and without anything
                 within 1 mutation of them
    the question where do the 10 end up when all 150 are ranked?

Against an expectation under random chance, not a correlation with no null
hypothesis. The idea is borrowed from the `jakob` branch of the shared repo.

    python scripts/retro_test.py --out results/retro.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10, load_facs, load_spr  # noqa: E402
from abx.model import onehot  # noqa: E402
from abx.precedent import blosum45_similarity  # noqa: E402
from abx.retro import (  # noqa: E402
    capture,
    edit_distances,
    leakage_mask,
    legal_pool,
    winner_mask,
    winner_ranks,
)

TOPK = (10, 20, 40)


def train_without(hidden: list[str], radius: int, seed: int) -> object:
    """Train the binder model on FACS without the winners and without their neighbours."""
    from abx.model import _fit  # noqa: PLC0415
    from abx.splits import distance_split, report_split  # noqa: PLC0415

    facs = load_facs()
    keep = facs.is_high | facs.is_low
    seqs = list(facs.cdrh3[keep])
    y = facs.is_high[keep].astype(np.float32)

    drop = leakage_mask(seqs, hidden, radius=radius)
    print(f"  removed {int(drop.sum()):,} FACS sequences within {radius} "
          f"mutation of a winner (of {len(seqs):,})")
    seqs = [s for s, d in zip(seqs, drop, strict=True) if not d]
    y = y[~drop]

    tr, te = distance_split(seqs, n_test=6000, min_dist=3, seed=seed)
    rep = report_split(seqs, tr, te)
    assert rep["min_dist"] >= 3, "the split leaks"
    X = onehot(seqs)
    net, log = _fit(X[tr], y[tr], X[te], y[te], seed=seed, verbose=False)
    print(f"  test AUC on FACS after the removal: {log['test_auroc']:.3f}")
    return net


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-edit", type=int, default=3)
    ap.add_argument("--threshold", type=float, default=1.0, help="nM, cutoff for winner")
    ap.add_argument("--radius", type=int, default=1, help="how much of the neighbourhood is removed")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "results/retro.json")
    args = ap.parse_args()

    spr = load_spr()
    pool, kd = legal_pool(spr.cdrh3, spr.kd_nM, spr.binding_result, args.min_edit)
    is_w = winner_mask(kd, args.threshold)
    winners = [s for s, w in zip(pool, is_w, strict=True) if w]

    print(f"pool: {len(pool)} measured variants with >= {args.min_edit} mutations")
    print(f"winners (< {args.threshold} nM): {len(winners)}")
    for s, k in sorted(zip(winners, kd[is_w], strict=True), key=lambda t: t[1]):
        print(f"    {s}  {k:.4f} nM  ({edit_distances([s])[0]} mutations)")

    print("\ntraining the binder model without them:")
    net = train_without(winners, args.radius, args.seed)

    scorers = {
        "edit_distance": -edit_distances(pool).astype(float),
        "blosum45": np.array([blosum45_similarity(s) for s in pool], dtype=float),
        "binder_net_retrained": np.asarray(net.logits(onehot(np.array(pool)))).ravel(),
    }

    out: dict[str, object] = {
        "pool_size": len(pool), "n_winners": int(is_w.sum()),
        "min_edit": args.min_edit, "threshold_nM": args.threshold,
        "leak_radius": args.radius, "winners": winners, "scorers": {},
    }
    width = max(len(k) for k in scorers)
    head = "  ".join(f"top{k:<3}" for k in TOPK)
    print(f"\n{'scorer':<{width}}  median rank  {head}")
    for name, sc in scorers.items():
        ranks = winner_ranks(sc, is_w)
        caps = {k: capture(sc, is_w, k) for k in TOPK}
        cells = "  ".join(
            f"{caps[k][0]}/{int(is_w.sum())} ({caps[k][0] / caps[k][1]:.1f}x)"
            for k in TOPK
        )
        print(f"{name:<{width}}  {float(np.median(ranks)):>10.1f}  {cells}")
        out["scorers"][name] = {  # type: ignore[index]
            "winner_ranks": [int(r) for r in ranks],
            "median_rank": float(np.median(ranks)),
            "capture": {str(k): {"hits": caps[k][0], "expected": round(caps[k][1], 2),
                                 "enrichment": round(caps[k][0] / caps[k][1], 2)}
                        for k in TOPK},
        }
    print(f"\n(random ranking gives a median rank of ~{len(pool) / 2:.0f} "
          f"and an enrichment of 1.0x)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
