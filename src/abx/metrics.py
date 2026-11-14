"""Small, exact metrics. No scipy, so there is nothing extra to install.

All of them are tested against known answers in tests/test_metrics.py.
"""

from __future__ import annotations

import numpy as np


def rankdata(x: np.ndarray) -> np.ndarray:
    """Average ranking, like scipy.stats.rankdata(method='average')."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    # distribute the average within each group of equal values
    sx = x[order]
    i = 0
    while i < len(sx):
        j = i
        while j + 1 < len(sx) and sx[j + 1] == sx[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return ranks


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman's rank correlation. 0 is random, 1 is a perfect ordering."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra, rb = rankdata(a[ok]), rankdata(b[ok])
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def auroc(score: np.ndarray, positive: np.ndarray) -> float:
    """Area under the ROC curve. The Mann–Whitney form, which handles ties correctly.

    0.5 is random. Below 0.5 means the score points the wrong way.
    """
    score = np.asarray(score, float)
    positive = np.asarray(positive, bool)
    ok = np.isfinite(score)
    score, positive = score[ok], positive[ok]
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = rankdata(score)
    return float((r[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def bootstrap_ci(fn, *arrays, n: int = 1000, alpha: float = 0.05,
                 seed: int = 0) -> tuple[float, float]:
    """95% interval by resampling. Without this, a single number is just a number."""
    rng = np.random.default_rng(seed)
    arrays = [np.asarray(a) for a in arrays]
    m = len(arrays[0])
    vals = []
    for _ in range(n):
        idx = rng.integers(0, m, m)
        v = fn(*[a[idx] for a in arrays])
        if np.isfinite(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    lo, hi = np.quantile(vals, [alpha / 2, 1 - alpha / 2])
    return (float(lo), float(hi))


def top_k_precision(score: np.ndarray, positive: np.ndarray, k: int) -> float:
    """Fraction of hits among the k highest-ranked.

    This is the measure that actually counts for us: we deliver 20 designs, not
    a ranking of half a million. AUROC can look good while the top is junk.
    """
    score = np.asarray(score, float)
    positive = np.asarray(positive, bool)
    k = min(k, len(score))
    top = np.argsort(-score, kind="mergesort")[:k]
    return float(positive[top].mean())
