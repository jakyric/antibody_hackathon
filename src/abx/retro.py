"""The winner-hidden test: would our recipe have found the known winners?

The idea is borrowed from the `jakob` branch in the shared repo. Every other
number we report says how well a scorer ranks *measured* variants. That is not
the task. The task is to pick a few candidates out of a space we have not
measured, and hope that the best ones are among them.

The test that imitates this: take the measured variants in the legal region,
hide the best ones, derive the scorer from the residue, and see where the
hidden ones end up when the whole region is ranked again. That gives a number
with an expectation under random chance to measure against, instead of a
correlation with no null hypothesis.

The easy thing to cheat with here is leakage: the FACS table contains
sequences that lie one mutation from the winners, and a model that has seen
them has in practice seen the winner. `leakage_mask` removes them.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from abx.data import TZ_CDRH3_10


def edit_distances(seqs: Sequence[str], reference: str = TZ_CDRH3_10) -> np.ndarray:
    """Hamming distance to the reference, for sequences of the same length as it."""
    return np.array(
        [sum(a != b for a, b in zip(reference, s, strict=True)) for s in seqs],
        dtype=int,
    )


def legal_pool(
    cdrh3: Sequence[str],
    kd_nM: Sequence[float],
    binding_result: Sequence[str],
    min_edit: int = 3,
) -> tuple[list[str], np.ndarray]:
    """The measured variants that lie where our designs lie.

    Only quantified KDs, only 10-mers, only distance >= `min_edit` from
    trastuzumab, and each sequence once (best measurement if measured several times).
    """
    best: dict[str, float] = {}
    for seq, kd, res in zip(cdrh3, kd_nM, binding_result, strict=True):
        if len(seq) != len(TZ_CDRH3_10) or res != "quantified_binding":
            continue
        if not np.isfinite(kd):
            continue
        if sum(a != b for a, b in zip(TZ_CDRH3_10, seq, strict=True)) < min_edit:
            continue
        if seq not in best or kd < best[seq]:
            best[seq] = float(kd)
    seqs = sorted(best)
    return seqs, np.array([best[s] for s in seqs], dtype=float)


def winner_mask(kd_nM: np.ndarray, threshold_nM: float = 1.0) -> np.ndarray:
    """Which of the variants in the pool we call winners. Stronger than the threshold."""
    return np.asarray(kd_nM, dtype=float) < threshold_nM


def leakage_mask(
    train_seqs: Sequence[str], hidden: Sequence[str], radius: int = 1
) -> np.ndarray:
    """True for training sequences that must be removed because they reveal a hidden winner.

    A sequence that lies `radius` mutations or closer to a winner carries
    almost the whole answer. Hiding the winner but keeping its neighbours is
    hiding nothing.
    """
    hidden = [h for h in hidden]
    if not hidden:
        return np.zeros(len(train_seqs), dtype=bool)
    n = len(TZ_CDRH3_10)
    hid = np.frombuffer("".join(hidden).encode(), dtype="S1").reshape(len(hidden), n)
    out = np.zeros(len(train_seqs), dtype=bool)
    for i, s in enumerate(train_seqs):
        if len(s) != n:
            continue
        q = np.frombuffer(s.encode(), dtype="S1")
        if int((hid != q).sum(axis=1).min()) <= radius:
            out[i] = True
    return out


def winner_ranks(scores: np.ndarray, is_winner: np.ndarray) -> np.ndarray:
    """The position of each winner when the whole pool is sorted. 1 = best.

    Equal scores get the same, pessimistic position: the worst they could have.
    Otherwise a scorer that gives everything the same number could look like
    it found the winners.
    """
    scores = np.asarray(scores, dtype=float)
    is_winner = np.asarray(is_winner, dtype=bool)
    # Number that are strictly better, plus number that are equal (including itself).
    ranks = np.array(
        [int((scores > s).sum()) + int((scores == s).sum()) for s in scores]
    )
    return np.sort(ranks[is_winner])


def capture(scores: np.ndarray, is_winner: np.ndarray, k: int) -> tuple[int, float]:
    """(winners among the top k, how many a random selection would give)."""
    is_winner = np.asarray(is_winner, dtype=bool)
    n, w = len(is_winner), int(is_winner.sum())
    ranks = winner_ranks(scores, is_winner)
    return int((ranks <= k).sum()), w * k / n


def enrichment(scores: np.ndarray, is_winner: np.ndarray, k: int) -> float:
    """How many times better than random. 1.0 = no information."""
    hit, exp = capture(scores, is_winner, k)
    return float("inf") if exp == 0 else hit / exp
