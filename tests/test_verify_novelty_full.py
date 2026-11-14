"""Tests for the exhaustive novelty verification.

The script's whole value is that it does not take a shortcut, so the tests
check the two things a shortcut would get wrong: a nearer neighbour of a
different length, and the loop being read in the tables' 10-residue
convention rather than Kabat's 13.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from verify_novelty_full import cdrh3_from_heavy, min_distance_exhaustive  # noqa: E402


def pack(seqs: list[str]) -> np.ndarray:
    if not seqs:
        return np.empty((0, 0), dtype=np.uint8)
    L = len(seqs[0])
    return np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)


def test_reads_the_ten_mer_not_the_kabat_thirteen_mer():
    heavy = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    assert cdrh3_from_heavy(heavy) == "WGGDGFYAMD"
    assert len(cdrh3_from_heavy(heavy)) == 10


def test_identical_sequence_is_distance_zero():
    d, nearest = min_distance_exhaustive("WGGDGFYAMD", pack(["WGGDGFYAMD"]), [])
    assert d == 0
    assert nearest == "WGGDGFYAMD"


def test_hamming_over_equal_lengths_finds_the_minimum():
    provided = ["WGGDGFYAMD", "AAADGFYAMD", "WGGDGFYAMA"]
    d, nearest = min_distance_exhaustive("WGGDGFYAMG", pack(provided), [])
    assert d == 1
    assert nearest in {"WGGDGFYAMD", "WGGDGFYAMA"}


def test_a_shorter_provided_sequence_can_be_the_nearest_neighbour():
    """The case Hamming cannot see.

    A 9-mer that is the query minus one residue sits at edit distance 1, while
    a Hamming-only comparison would score it 99 and declare the design novel.
    """
    query = "WGGDGFYAMD"
    far_same_length = pack(["AAAAAAAAAA"])
    d, nearest = min_distance_exhaustive(query, far_same_length, ["WGGDGFYAM"])
    assert d == 1
    assert nearest == "WGGDGFYAM"


def test_length_gap_shortcut_never_hides_a_closer_neighbour():
    """Skipping on |len difference| >= best must not change the answer."""
    query = "WGGDGFYAMD"
    others = ["WGGDGFYAM", "WGG", "WGGDGFYAMDXY", "WGGDGFYAMDX"]
    d, _ = min_distance_exhaustive(query, pack([]), others)
    brute = min(_levenshtein(query, o) for o in others)
    assert d == brute == 1


def test_equal_length_neighbour_reached_by_a_shift_is_an_edit_distance():
    """The case positional comparison cannot see, and the authority can.

    `YPGRYFYTYD` and `YPAGRFYTYD` are both 10-mers and differ at three
    positions, so Hamming calls them 3 apart -- exactly the floor, which reads
    as a pass. Aligned, the second is the first
    with an A inserted and a Y deleted: two edits. The official checker aligns,
    so a design we score at 3 can be rejected at 2, and a rejected design is
    never tested.
    """
    query = "YPGRYFYTYD"
    neighbour = "YPAGRFYTYD"
    assert sum(a != b for a, b in zip(query, neighbour)) == 3
    assert _levenshtein(query, neighbour) == 2

    d, nearest = min_distance_exhaustive(query, pack([neighbour]), [])
    assert d == 2, "equal lengths must still be compared by edit distance"
    assert nearest == neighbour


def test_equal_length_edit_distance_matches_the_slow_reference():
    rng = random.Random(0)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    query = "".join(rng.choice(alphabet) for _ in range(10))
    provided = ["".join(rng.choice(alphabet) for _ in range(10)) for _ in range(200)]
    provided += [query[:3] + query[4:] + rng.choice(alphabet)]  # a shift neighbour
    d, _ = min_distance_exhaustive(query, pack(provided), [])
    assert d == min(_levenshtein(query, p) for p in provided)


def test_control_that_is_in_the_data_is_not_mistaken_for_novel():
    provided = ["WGGDGFYAMD", "MGPEPFYAMD"]
    d, _ = min_distance_exhaustive("WGGDGFYAMD", pack(provided), [])
    assert d == 0, "trastuzumab must collide with the provided data"


def _levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]
