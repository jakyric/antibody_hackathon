"""The gate: the hard filters, mirrored locally.

The official checker is the source of truth. This is a local copy that lets us
discard impossible designs before we spend GPU time on them, and that makes the
rules testable. Where we cannot reproduce the official checker, we say so instead of guessing.

The rules, from the task description:

  common       format VHH, scFv or Fab
               no unpaired cysteines
               total length at most 250 amino acids
               numbering with at least 0.1 confidence in `immunum`

  optimisation below 70 % CDR-H3 identity to trastuzumab
               AND strictly more than 2 point mutations from every single
               data point that has been provided

  de novo      below 70 % identity to any known antibody
               AND at least 2 CDR-H3 mutations from trastuzumab

WHAT WE CANNOT CHECK HERE
  - The confidence of the `immunum` numbering. Requires their library.
  - The de novo track's "known antibody" database. We have no such database
    locally, so de novo novelty can only be decided by the official checker.
  Both are reported as UNKNOWN, never as passed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .data import TZ_CDRH3_10, provided_cdrh3s

MAX_TOTAL_LENGTH = 250
CDRH3_IDENTITY_CEILING = 0.70      # must be STRICTLY below
MIN_MUTATIONS_FROM_PROVIDED = 3    # "strictly more than 2"

AA = set("ACDEFGHIKLMNPQRSTVWY")


def levenshtein(a: str, b: str) -> int:
    """Exact edit distance. Used when the lengths differ."""
    if a == b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("hamming requires equal length")
    return sum(x != y for x, y in zip(a, b))


def identity(a: str, b: str) -> float:
    """Fraction of identical positions. Unequal lengths are handled with edit distance."""
    if not a or not b:
        return 0.0
    if len(a) == len(b):
        return 1.0 - hamming(a, b) / len(a)
    return 1.0 - levenshtein(a, b) / max(len(a), len(b))


def count_unpaired_cysteines(seq: str) -> int:
    """An unpaired cysteine gives a free thiol: scrambling and aggregation.

    Without a structure we cannot see which ones actually pair, so we
    use the rule everyone uses: an odd number of cysteines means at least one is free.
    """
    n = seq.count("C")
    return n % 2


class ProvidedIndex:
    """All provided CDR-H3s, with fast minimum distance.

    525,063 sequences of length 9–10.

    Two ways in:

    `min_distance`  exact, but scans the whole set. Fine for a handful
                    of queries, hopeless for twenty thousand.

    `is_novel`      answers only what we need — "is this more than
                    2 mutations from everything provided?" — and uses
                    the pigeonhole principle: if a 10-mer is split into three blocks,
                    two mutations cannot possibly touch all three. Two sequences
                    within distance 2 must therefore share at least one block
                    exactly. We look up the three blocks and compare
                    only against those that share one. Exact, and a thousand times
                    faster.

    That the two give the same answer is tested in test_gate.py.
    """

    BLOCKS = ((0, 4), (4, 7), (7, 10))     # three blocks of a 10-mer

    def __init__(self, sequences: set[str] | None = None):
        seqs = sequences if sequences is not None else provided_cdrh3s()
        self.n = len(seqs)
        self._by_len: dict[int, np.ndarray] = {}
        self._raw: dict[int, list[str]] = {}
        by_len: dict[int, list[str]] = {}
        for s in seqs:
            by_len.setdefault(len(s), []).append(s)
        for L, group in by_len.items():
            arr = np.frombuffer("".join(group).encode(), dtype=np.uint8)
            self._by_len[L] = arr.reshape(len(group), L)
            self._raw[L] = group

        # block index, only for 10-mers (the rest are few enough to scan)
        self._blocks: list[dict[str, list[str]]] = [{} for _ in self.BLOCKS]
        for s in self._raw.get(10, []):
            for bi, (a, b) in enumerate(self.BLOCKS):
                self._blocks[bi].setdefault(s[a:b], []).append(s)
        self._short = [s for L, g in self._raw.items() if L != 10 for s in g]

    def is_novel(self, query: str, min_mutations: int = MIN_MUTATIONS_FROM_PROVIDED) -> bool:
        """Is `query` at least `min_mutations` from everything provided?

        Fast path for 10-mers and the requirement "more than 2". Falls back to
        the exact scan for everything else.
        """
        if len(query) != 10 or min_mutations != 3:
            return self.min_distance(query) >= min_mutations

        seen: set[str] = set()
        for bi, (a, b) in enumerate(self.BLOCKS):
            for cand in self._blocks[bi].get(query[a:b], ()):
                if cand in seen:
                    continue
                seen.add(cand)
                if sum(x != y for x, y in zip(query, cand)) < 3:
                    return False
        # shorter provided sequences: a 9-mer is at least 1 away from a 10-mer,
        # so only the very close ones can be a problem
        for cand in self._short:
            if levenshtein(query, cand) < 3:
                return False
        return True

    def min_distance(self, query: str) -> int:
        """Smallest number of point mutations from `query` to anything provided.

        Equal length gives Hamming. Unequal length gives edit distance, but
        only against lengths that can lie close enough to matter.
        """
        best = len(query)
        q = np.frombuffer(query.encode(), dtype=np.uint8)
        for L, mat in self._by_len.items():
            if L == len(query):
                d = int((mat != q).sum(axis=1).min())
                best = min(best, d)
            elif abs(L - len(query)) < best:
                for cand in self._raw[L]:
                    d = levenshtein(query, cand)
                    if d < best:
                        best = d
                        if best == 0:
                            return 0
        return best


@dataclass
class Verdict:
    design_id: str
    passes: bool
    track: str                       # 'optimisation' | 'de novo' | 'rejected'
    checks: dict[str, object] = field(default_factory=dict)
    reasons: list[str] = field(default_factory=list)
    unknown: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        mark = "OK " if self.passes else "REJECTED"
        out = f"{mark} {self.design_id} [{self.track}]"
        for r in self.reasons:
            out += f"\n    rejected: {r}"
        for u in self.unknown:
            out += f"\n    unknown: {u}"
        return out


def check_design(
    design_id: str,
    sequence: str,
    cdrh3: str | None = None,
    fmt: str = "scFv",
    provided: ProvidedIndex | None = None,
) -> Verdict:
    """Run all the filters we can run locally.

    `sequence` is the whole construct (for Fab: VH+VL merged).
    `cdrh3` in the dataset's 10-mer convention; see data.py.
    """
    v = Verdict(design_id=design_id, passes=True, track="optimisation")
    seq = sequence.upper().strip()

    # ── common ──────────────────────────────────────────────────────────
    v.checks["length"] = len(seq)
    if len(seq) > MAX_TOTAL_LENGTH:
        v.passes = False
        v.reasons.append(f"length {len(seq)} > {MAX_TOTAL_LENGTH}")

    bad = sorted(set(seq) - AA)
    if bad:
        v.passes = False
        v.reasons.append(f"unknown characters: {''.join(bad)}")

    unpaired = count_unpaired_cysteines(seq)
    v.checks["cysteines"] = seq.count("C")
    v.checks["unpaired_cysteine"] = bool(unpaired)
    if unpaired:
        v.passes = False
        v.reasons.append(f"odd number of cysteines ({seq.count('C')}): free thiol")

    if fmt not in ("VHH", "scFv", "Fab"):
        v.passes = False
        v.reasons.append(f"format {fmt!r} is not VHH, scFv or Fab")
    v.checks["format"] = fmt

    v.unknown.append("confidence of the immunum numbering — only the official checker can decide")

    # ── novelty, optimisation track ─────────────────────────────────────
    if cdrh3:
        ident = identity(cdrh3, TZ_CDRH3_10)
        v.checks["cdrh3"] = cdrh3
        v.checks["cdrh3_identity_to_trastuzumab"] = round(ident, 4)
        v.checks["cdrh3_mutations_from_trastuzumab"] = (
            hamming(cdrh3, TZ_CDRH3_10) if len(cdrh3) == len(TZ_CDRH3_10)
            else levenshtein(cdrh3, TZ_CDRH3_10)
        )
        if ident >= CDRH3_IDENTITY_CEILING:
            v.passes = False
            v.reasons.append(
                f"CDR-H3 identity {ident:.0%} to trastuzumab is not below 70 %"
            )

        if provided is not None:
            # the pigeonhole path: answers what the requirement asks without scanning
            # the whole set. Exact distance is only fetched when something fails.
            if provided.is_novel(cdrh3):
                v.checks["min_mutations_from_provided"] = ">=3"
            else:
                d = provided.min_distance(cdrh3)
                v.checks["min_mutations_from_provided"] = d
                v.passes = False
                v.reasons.append(
                    f"only {d} point mutation(s) from a provided data point, "
                    f"the requirement is more than 2"
                )
        else:
            v.unknown.append("distance to provided data — no index given")
    else:
        v.unknown.append("no CDR-H3 given, novelty not assessed")

    if not v.passes:
        v.track = "rejected"
    return v


def min_cdrh3_mutations_required(cdrh3_length: int = 10) -> int:
    """How many mutations are needed to get below 70 % identity.

    For a 10-mer: 7 hits gives 70 %, which is not below. 6 hits gives 60 %.
    So at least 4 mutations. This number decides the whole design strategy,
    because the share of good binders falls steeply with distance.
    """
    for k in range(cdrh3_length + 1):
        if (cdrh3_length - k) / cdrh3_length < CDRH3_IDENTITY_CEILING:
            return k
    return cdrh3_length
