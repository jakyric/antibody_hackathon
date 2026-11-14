"""Honest splits of the FACS library.

WHY THIS FILE EXISTS
  The first attempt split on the first five residues. It looked right and
  was not. The library is a near-exhaustive mutagenesis: all ten
  positions carry all twenty amino acids. Two sequences that differ at
  position 1 end up in separate groups, yet lie one mutation apart
  from each other.

  Measured: 194 of 300 test sequences lay 1 mutation from a training
  sequence, and another 95 lay 2 away. Median distance 1. A network got
  AUROC 0.995 on such a split and Spearman -0.04 against real SPR
  measurements. The number was memorization of neighborhood.

  `distance_split` instead removes everything in the training part that
  lies closer than `min_dist` to a test sequence. The test then measures
  what it should: whether the model can say something about sequences it
  has not seen the like of.

  The organizers' own clonotype split does the same thing, and does it
  better; see `load_official_split`.
"""

from __future__ import annotations

import numpy as np

from .model import AA, AA_INDEX, L


def _onehot_u8(seqs) -> np.ndarray:
    X = np.zeros((len(seqs), L * len(AA)), dtype=np.float32)
    for i, s in enumerate(seqs):
        for p, ch in enumerate(s[:L]):
            j = AA_INDEX.get(ch)
            if j is not None:
                X[i, p * len(AA) + j] = 1.0
    return X


def min_hamming_to_set(queries, reference, chunk: int = 400) -> np.ndarray:
    """Smallest Hamming distance from each of `queries` to anything in `reference`.

    Computed as 10 minus the number of hits, where hits come from an
    inner product between one-hot encodings. BLAS does the work, so it
    takes seconds instead of hours.
    """
    Q = _onehot_u8(list(queries))
    R = _onehot_u8(list(reference))
    out = np.empty(len(Q), dtype=np.int16)
    for s in range(0, len(Q), chunk):
        blokk = Q[s:s + chunk] @ R.T          # (chunk, N) number of identical positions
        out[s:s + chunk] = (L - blokk.max(axis=1)).astype(np.int16)
    return out


def distance_split(
    seqs,
    n_test: int = 6000,
    min_dist: int = 3,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Train/test where no test sequence has a close relative in training.

    Returns (train mask, test mask). The training part shrinks: that is
    the price of a number you can trust.
    """
    seqs = list(seqs)
    n = len(seqs)
    rng = np.random.default_rng(seed)
    test_idx = rng.choice(n, size=min(n_test, n), replace=False)
    test_mask = np.zeros(n, dtype=bool)
    test_mask[test_idx] = True

    rest = np.flatnonzero(~test_mask)
    d = min_hamming_to_set([seqs[i] for i in rest],
                           [seqs[i] for i in test_idx])
    train_mask = np.zeros(n, dtype=bool)
    train_mask[rest[d >= min_dist]] = True
    return train_mask, test_mask


def report_split(seqs, train_mask, test_mask, sample: int = 300,
                 seed: int = 0) -> dict:
    """Measure the actual distance between the splits. A number, not an assurance."""
    seqs = list(seqs)
    rng = np.random.default_rng(seed)
    ti = np.flatnonzero(test_mask)
    sub = rng.choice(ti, size=min(sample, len(ti)), replace=False)
    d = min_hamming_to_set([seqs[i] for i in sub],
                           [seqs[i] for i in np.flatnonzero(train_mask)])
    return {
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "min_dist": int(d.min()),
        "median_dist": float(np.median(d)),
        "frac_within_2": float((d <= 2).mean()),
    }

