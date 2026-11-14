"""The benchmark: does our ranking work on answers we already know?

WHY THIS EXISTS
  The three properties the task is judged on — affinity, epitope
  conservation and thermostability — are measured in the lab, weeks after
  we have delivered. We cannot measure them today. The only honest thing
  we can do is take data where the answer has already been measured, keep
  it away from the model, and see whether our ranking finds the answer
  again.

  Without this, "the model gave a high score" is just a claim.

THE CRUCIAL DETAIL
  The novelty filter requires at least 4 mutations in CDR-H3. The share of
  good binders falls from 97% at 2 mutations to 75% at 4 and 34% at 7.
  Ranking on the whole dataset is therefore far too easy a task: edit
  distance alone solves it. What we need to know is whether we can rank
  WITHIN the legal range. That is why everything here reports two
  numbers: all, and d >= 4.

BEATING THE BASELINE
  Two baselines are included. If your model cannot beat
  `edit_distance`, then use `edit_distance` and say so plainly. It is a
  perfectly good answer, and far better than dressing up a model that
  does not work.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from .data import TZ_CDRH3_10, Facs, Spr, load_facs, load_spr
from .gate import hamming, min_cdrh3_mutations_required
from .metrics import auroc, bootstrap_ci, spearman, top_k_precision

# A scorer takes CDR-H3 sequences and gives one number per sequence.
# Higher should mean a better binder.
Scorer = Callable[[Sequence[str]], np.ndarray]

MIN_LEGAL_MUTATIONS = min_cdrh3_mutations_required(10)   # = 4


# ── baselines ───────────────────────────────────────────────────────────────
def baseline_random(seed: int = 0) -> Scorer:
    def score(seqs: Sequence[str]) -> np.ndarray:
        rng = np.random.default_rng(seed)
        return rng.random(len(seqs))
    return score


def baseline_edit_distance(seqs: Sequence[str]) -> np.ndarray:
    """Closer to trastuzumab = better. Surprisingly strong, and the one you must beat."""
    return np.array([
        -(hamming(s, TZ_CDRH3_10) if len(s) == len(TZ_CDRH3_10) else 99)
        for s in seqs
    ], dtype=float)


# ── group split ─────────────────────────────────────────────────────────────
def grouped_split(seqs: Sequence[str], frac_test: float = 0.3,
                  seed: int = 0, prefix: int = 5) -> np.ndarray:
    """Split into train and test so that similar sequences end up in the same place.

    We group on the first `prefix` residues. Sequences that share a
    beginning are related, and end up in the same split. Without this the
    dataset leaks: nearly identical variants would sit on both sides,
    and any model would look brilliant.

    THIS IS AN APPROXIMATION. The Chinery/Mehta data has a real
    clonotype split, made precisely for this, in the bucket under
        zenodo-10831512-trastuzumab-ngs/.../clonotype_split/0.7_0.15_0.15/
    Use it once it has been downloaded. This one is for getting started.

    Returns a boolean mask: True = test split.
    """
    out = np.zeros(len(seqs), dtype=bool)
    for i, s in enumerate(seqs):
        key = f"{seed}:{s[:prefix]}".encode()
        h = int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big")
        out[i] = (h % 10_000) < frac_test * 10_000
    return out


# ── reports ─────────────────────────────────────────────────────────────────
@dataclass
class Result:
    name: str
    n: int
    value: float
    ci: tuple[float, float]
    extra: dict[str, float]


def _fmt(r: Result) -> str:
    lo, hi = r.ci
    ex = "  ".join(f"{k} {v:.3f}" for k, v in r.extra.items())
    return (f"  {r.name:<22} n={r.n:>6}  {r.value:+.3f}  "
            f"[{lo:+.3f}, {hi:+.3f}]  {ex}")


def evaluate_on_spr(scorer: Scorer, name: str, spr: Spr | None = None,
                    legal_only: bool = False) -> Result:
    """Does the scorer rank the 182 SPR variants correctly by measured KD?

    The target is -log10(KD): higher is tighter binding. Spearman, because
    we care about the order, not the absolute values — and because the guide
    is clear that order is more reliable than absolute numbers.
    """
    spr = spr if spr is not None else load_spr()
    keep = spr.has_kd & (spr.cdrh3 != TZ_CDRH3_10)
    if legal_only:
        muts = np.array([
            hamming(s, TZ_CDRH3_10) if len(s) == len(TZ_CDRH3_10) else 99
            for s in spr.cdrh3
        ])
        keep &= muts >= MIN_LEGAL_MUTATIONS

    seqs = spr.cdrh3[keep]
    target = -np.log10(spr.kd_nM[keep])
    s = np.asarray(scorer(list(seqs)), dtype=float)

    rho = spearman(s, target)
    ci = bootstrap_ci(spearman, s, target)
    # share of true top binders among the 20 we would deliver
    top20 = top_k_precision(s, target >= np.quantile(target, 0.80), 20)
    return Result(name, int(keep.sum()), rho, ci, {"top20_hits": top20})


def evaluate_on_facs(scorer: Scorer, name: str, facs: Facs | None = None,
                     legal_only: bool = False, test_only: bool = True,
                     subsample: int | None = 40_000, seed: int = 0) -> Result:
    """Does the scorer separate "high" from "low" on data it has not seen?

    'medium' is left out on purpose: the class is the most ambiguous, and
    including it makes the number harder to interpret without making it
    more true.
    """
    facs = facs if facs is not None else load_facs()
    keep = facs.is_high | facs.is_low
    if legal_only:
        keep &= facs.edit_distance >= MIN_LEGAL_MUTATIONS
    if test_only:
        keep &= grouped_split(facs.cdrh3, seed=seed)

    idx = np.flatnonzero(keep)
    if subsample is not None and len(idx) > subsample:
        idx = np.random.default_rng(seed).choice(idx, subsample, replace=False)

    seqs = facs.cdrh3[idx]
    pos = facs.is_high[idx]
    s = np.asarray(scorer(list(seqs)), dtype=float)

    a = auroc(s, pos)
    ci = bootstrap_ci(auroc, s, pos)
    return Result(name, len(idx), a, ci,
                  {"top20_hits": top_k_precision(s, pos, 20),
                   "share_high": float(pos.mean())})


def run_all(scorers: dict[str, Scorer] | None = None, seed: int = 0) -> None:
    """Run the whole benchmark and print the table."""
    scorers = scorers or {
        "random": baseline_random(seed),
        "edit_distance": baseline_edit_distance,
    }

    facs, spr = load_facs(), load_spr()
    print("=" * 78)
    print("BENCHMARK — do we rank correctly on answers we already know?")
    print("=" * 78)
    print(f"FACS {len(facs)} variants · SPR {int(spr.has_kd.sum())} with measured KD")
    print(f"The novelty filter requires at least {MIN_LEGAL_MUTATIONS} mutations in CDR-H3.")
    print("Numbers in square brackets are 95% intervals from 1000 resamples.\n")

    for legal in (False, True):
        tag = f"LEGAL ONLY (d >= {MIN_LEGAL_MUTATIONS})" if legal else "ALL"
        print(f"── {tag} " + "─" * (72 - len(tag)))

        print(" SPR · Spearman against -log10(KD), higher is better")
        for nm, sc in scorers.items():
            print(_fmt(evaluate_on_spr(sc, nm, spr, legal_only=legal)))

        print(" FACS · AUROC high against low, held out by group")
        for nm, sc in scorers.items():
            print(_fmt(evaluate_on_facs(sc, nm, facs, legal_only=legal, seed=seed)))
        print()

    print("Read it like this: a number near 0 (Spearman) or 0.5 (AUROC) means")
    print("that the scorer knows nothing. If you do not beat edit_distance in the")
    print("legal zone, then submit edit_distance and say why.")


if __name__ == "__main__":
    run_all()
