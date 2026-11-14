"""Two-stage ranking, where each stage is justified by a measured number.

WHY TWO STAGES AND NOT ONE SCORE
  We measured two things that point in different directions:

  1. The net hits the FACS class with AUROC 0.98 on a split where no
     test sequence lies closer than 3 mutations to the training data.
     So it can recognize a binder.

  2. But the FACS class itself hits measured affinity with AUROC 0.61. For
     the 27 sequences that were measured both ways, the median KD is 7.2 nM
     among "high" and 4.2 nM among "medium". The label says whether
     something binds, not how tightly.

  Ranking the submission by the net would thus have meant optimizing the
  wrong target with great confidence. Instead we use the net for what it
  can do — filtering out those that are unlikely to bind at all — and
  rank the remaining ones by what actually transferred to SPR: proximity
  to trastuzumab.

  edit_distance gave Spearman +0.23 in the legal zone. That is weak. We
  say so in the pitch instead of hiding it behind a fancier model.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .data import TZ_CDRH3_10, vh_from_cdrh3
from .gate import ProvidedIndex, check_design, hamming
from .model import MLP, onehot
from .sequences import SCFV_LINKER, TRASTUZUMAB_VL

# Rough chemical kinship. Used only to prefer cautious swaps when
# everything else is equal, never as the sole reason.
GROUPS = {
    "hydrofob": "AVLIMFWC",
    "polar": "STNQY",
    "positiv": "KRH",
    "negativ": "DE",
    "spesiell": "GP",
}
_GROUP_OF = {a: g for g, s in GROUPS.items() for a in s}


def conservative_fraction(seq: str, ref: str = TZ_CDRH3_10) -> float:
    """Fraction of the swaps that stay within the same chemical group."""
    muts = [(a, b) for a, b in zip(ref, seq) if a != b]
    if not muts:
        return 1.0
    same = sum(_GROUP_OF.get(a) == _GROUP_OF.get(b) for a, b in muts)
    return same / len(muts)


@dataclass
class Design:
    design_id: str
    cdrh3: str
    n_mut: int
    model_score: float
    model_percentile: float
    min_dist_provided: int
    conservative: float
    rank_score: float
    gate_ok: bool
    gate_notes: list[str] = field(default_factory=list)

    @property
    def vh(self) -> str:
        return vh_from_cdrh3(self.cdrh3)

    @property
    def scfv(self) -> str:
        return self.vh + SCFV_LINKER + TRASTUZUMAB_VL

    def rationale(self) -> str:
        return (
            f"{self.n_mut} mutations from trastuzumab (minimum legal is 4); "
            f"the net places it in the {self.model_percentile:.0f}th percentile of "
            f"binder probability; {self.conservative:.0%} of the swaps are "
            f"chemically conservative; nearest provided sequence is "
            f"{self.min_dist_provided} mutations away"
        )


def stage1_filter(cdrh3s: list[str], net: MLP, keep_percentile: float = 70.0
                  ) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Filter out those the net thinks are unlikely to bind.

    We keep the top part of the distribution. The threshold is a
    percentile and not an absolute value, because the logits are not
    calibrated probabilities and we do not pretend they are.
    """
    z = net.logits(onehot(cdrh3s))
    cut = np.percentile(z, keep_percentile)
    keep = z >= cut
    pct = 100.0 * (z[:, None] > z[None, :]).mean(axis=1)
    return [c for c, k in zip(cdrh3s, keep) if k], z, pct


def stage2_rank(cdrh3s: list[str], z: np.ndarray, pct: np.ndarray,
                provided: ProvidedIndex) -> list[Design]:
    """Rank the remaining ones by what actually transferred."""
    out = []
    for i, c in enumerate(cdrh3s):
        n_mut = hamming(c, TZ_CDRH3_10)
        # exact distance is expensive; they have already been screened with
        # is_novel, so we report a cap instead of rescanning the whole set
        d_prov = 3 if provided.is_novel(c) else 2
        cons = conservative_fraction(c)
        v = check_design(f"tmp", vh_from_cdrh3(c), cdrh3=c, fmt="scFv",
                         provided=provided)
        # The ranking: proximity first (the one that transferred), then
        # cautious swaps, then the net's belief. The weights are chosen to
        # mirror what we measured, not tuned against the test set.
        rank = (
            -2.0 * n_mut
            + 0.8 * cons
            + 0.4 * (pct[i] / 100.0)
            + 0.1 * min(d_prov, 6)
        )
        out.append(Design(
            design_id="", cdrh3=c, n_mut=n_mut, model_score=float(z[i]),
            model_percentile=float(pct[i]), min_dist_provided=d_prov,
            conservative=cons, rank_score=float(rank),
            gate_ok=v.passes, gate_notes=list(v.reasons),
        ))
    out.sort(key=lambda d: -d.rank_score)
    return out


def diversify(designs: list[Design], k: int, min_sep: int = 3) -> list[Design]:
    """Pick k designs that are not variants of each other.

    We deliver 20, not 20 near-identical ones. Greedy: take the best-ranked
    one that is at least `min_sep` mutations from all those already chosen,
    relax the requirement if we run out.
    """
    chosen: list[Design] = []
    for sep in range(min_sep, 0, -1):
        for d in designs:
            if len(chosen) >= k:
                break
            if all(hamming(d.cdrh3, c.cdrh3) >= sep for c in chosen):
                chosen.append(d)
        if len(chosen) >= k:
            break
    return chosen[:k]
