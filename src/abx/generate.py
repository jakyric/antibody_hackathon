"""Proposals for legal CDR-H3s.

THE STRATEGY, AND WHY
  The novelty filter requires at least 4 mutations from trastuzumab. The share
  of good binders in the FACS library falls steeply with distance:

      2 mutations  97 %      4 mutations  75 %      7 mutations  34 %

  We are not allowed to stay close, and it does not pay to go far. We
  therefore settle on 4 to 5 mutations: as close to the original as the rules
  allow.

  Which positions tolerate change is not guessed. It is calculated
  from 524,346 measured variants, as a position weight matrix: for each
  position and each amino acid, the share of "high" among those that carry it.
"""

from __future__ import annotations

import numpy as np

from .data import TZ_CDRH3_10, load_facs
from .model import AA, AA_INDEX, L


def position_weight_matrix(min_count: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """(10, 20) share of "high", and (10, 20) number of observations.

    This is the interpretable version of what the network has learned, and the one we can
    show in a pitch: each cell is a measured share, not a weight.
    """
    f = load_facs()
    keep = f.is_high | f.is_low
    seqs, hi = f.cdrh3[keep], f.is_high[keep]
    arr = np.frombuffer("".join(seqs).encode(), np.uint8).reshape(len(seqs), L)
    codes = np.array([ord(a) for a in AA], dtype=np.uint8)

    frac = np.full((L, len(AA)), np.nan)
    cnt = np.zeros((L, len(AA)), dtype=int)
    for p in range(L):
        col = arr[:, p]
        for j, c in enumerate(codes):
            m = col == c
            n = int(m.sum())
            cnt[p, j] = n
            if n >= min_count:
                frac[p, j] = float(hi[m].mean())
    return frac, cnt


def tolerant_positions(frac: np.ndarray) -> np.ndarray:
    """How much a position tolerates being swapped, measured as spread in outcomes.

    A position where almost every swap gives low binding is intolerant, and
    should be left alone. One where many swaps go well is a place to look.
    """
    return np.nanmean(frac, axis=1)


class Generator:
    """Draws legal variants, weighted by what the data actually shows."""

    def __init__(self, frac: np.ndarray | None = None, seed: int = 0):
        if frac is None:
            frac, _ = position_weight_matrix()
        self.frac = frac
        self.rng = np.random.default_rng(seed)
        # per position: probability of each amino acid, based on the measured share of "high"
        self.aa_p = np.zeros_like(frac)
        for p in range(L):
            row = np.nan_to_num(frac[p], nan=0.0)
            # trastuzumab's own residue is not a mutation
            row = row.copy()
            row[AA_INDEX[TZ_CDRH3_10[p]]] = 0.0
            row = np.clip(row, 0, None) ** 3        # sharpen towards the good ones
            s = row.sum()
            self.aa_p[p] = row / s if s > 0 else np.full(len(AA), 1 / len(AA))
        self.pos_p = tolerant_positions(frac)
        self.pos_p = np.nan_to_num(self.pos_p, nan=0.0)
        self.pos_p = self.pos_p / self.pos_p.sum()

    def draw(self, n_mut: int) -> str:
        pos = self.rng.choice(L, size=n_mut, replace=False, p=self.pos_p)
        out = list(TZ_CDRH3_10)
        for p in pos:
            out[p] = AA[self.rng.choice(len(AA), p=self.aa_p[p])]
        return "".join(out)

    def sample(self, n: int, n_mut_choices=(4, 5), max_tries: int = 200_000) -> list[str]:
        """n unique variants, each with 4 or 5 mutations from trastuzumab."""
        seen: set[str] = set()
        tries = 0
        while len(seen) < n and tries < max_tries:
            tries += 1
            k = int(self.rng.choice(n_mut_choices))
            s = self.draw(k)
            # the draw may by chance give back an original residue
            if sum(a != b for a, b in zip(s, TZ_CDRH3_10)) < min(n_mut_choices):
                continue
            seen.add(s)
        return sorted(seen)

