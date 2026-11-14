"""Tests for the exhaustive novelty verification.

The script's whole value is that it does not take a shortcut, so the tests
check the two things a shortcut would get wrong: a nearer neighbour of a
different length, and the loop being read in the tables' 10-residue
convention rather than Kabat's 13.
"""

from __future__ import annotations

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
