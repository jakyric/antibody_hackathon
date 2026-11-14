"""Tests the gate against answers we know in advance.

The point is not that the code runs. The point is that the rules are
correctly understood. Each test below is a claim about the task text that
may be wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from abx.data import (TZ_CDRH3_10, cdrh3_from_vh, load_facs, load_spr,
                      provided_cdrh3s, vh_from_cdrh3)
from abx.gate import (ProvidedIndex, check_design, hamming, identity,
                      levenshtein, min_cdrh3_mutations_required)
from abx.metrics import auroc, spearman, top_k_precision
from abx.sequences import TRASTUZUMAB_CDRS, TRASTUZUMAB_VH, TRASTUZUMAB_VL


# ── the sequences match the structure they came from ────────────────────────
def test_cdrs_exist_in_the_sequences():
    for navn, cdr in TRASTUZUMAB_CDRS.items():
        kjede = TRASTUZUMAB_VH if navn.startswith("H") else TRASTUZUMAB_VL
        assert cdr in kjede, f"{navn} ({cdr}) is not in the chain"


def test_dataset_cdrh3_is_one_residue_shorter():
    """The dataset drops the final Y. If the conventions are mixed up,
    every distance is off by one."""
    assert TRASTUZUMAB_CDRS["H3"] == "WGGDGFYAMDY"      # Kabat, 11
    assert TZ_CDRH3_10 == "WGGDGFYAMD"                   # dataset, 10
    assert TRASTUZUMAB_CDRS["H3"] == TZ_CDRH3_10 + "Y"


def test_cdrh3_can_be_extracted_and_reinserted():
    assert cdrh3_from_vh(TRASTUZUMAB_VH) == TZ_CDRH3_10
    for loop in ["ABCDEFGHIK", "WGHCCLYQLG"]:
        assert cdrh3_from_vh(vh_from_cdrh3(loop)) == loop


# ── the threshold that drives the whole design strategy ─────────────────────
def test_four_mutations_is_the_minimum():
    """70 % of 10 is 7. Seven hits is not BELOW 70 %, six is. So four."""
    assert min_cdrh3_mutations_required(10) == 4
    assert identity("WGGDGFYAMD", "WGGDGFYAMD") == 1.0
    # TZ = W G G D G F Y A M D
    assert identity("WGGDGFYKRE", TZ_CDRH3_10) == pytest.approx(0.7)   # 3 mut, TOO SIMILAR
    assert identity("WGGDGFKRES", TZ_CDRH3_10) == pytest.approx(0.6)   # 4 mut, legal


def test_three_mutations_rejected_four_pass():
    tre = "WGGDGFYKRE"
    fire = "WGGDGFKRES"
    assert hamming(tre, TZ_CDRH3_10) == 3
    assert hamming(fire, TZ_CDRH3_10) == 4

    v = check_design("tre", vh_from_cdrh3(tre), cdrh3=tre, fmt="scFv")
    assert not v.passes
    assert any("70" in r for r in v.reasons)

    v = check_design("fire", vh_from_cdrh3(fire), cdrh3=fire, fmt="scFv")
    assert "cdrh3_identity_to_trastuzumab" in v.checks
    assert v.checks["cdrh3_identity_to_trastuzumab"] < 0.70


# ── the hard filters ────────────────────────────────────────────────────────
def test_lone_cysteine_rejected():
    med = "WGGDGFKREC"
    v = check_design("odd", vh_from_cdrh3(med), cdrh3=med, fmt="scFv")
    seq = vh_from_cdrh3(med)
    if seq.count("C") % 2 == 1:
        assert not v.passes
        assert any("thiol" in r for r in v.reasons)


def test_too_long_rejected():
    v = check_design("lang", "A" * 251, cdrh3=None, fmt="scFv")
    assert not v.passes
    assert any("length" in r for r in v.reasons)


def test_unknown_format_rejected():
    v = check_design("x", "ACDEFGHIKLM", cdrh3=None, fmt="IgG")
    assert not v.passes


def test_immunum_reported_as_unknown_not_passed():
    """We cannot reproduce the official checker's numbering. That should be stated, not hidden."""
    v = check_design("x", vh_from_cdrh3("WGGDGFKRES"), cdrh3="WGGDGFKRES")
    assert any("immunum" in u for u in v.unknown)


# ── novelty against the provided data ───────────────────────────────────────
@pytest.fixture(scope="module")
def provided():
    return ProvidedIndex()


def test_number_of_provided_matches_readme(provided):
    """The data README says 525,063 unique CDR-H3s."""
    assert provided.n == 525_063


def test_trastuzumab_itself_is_not_novel(provided):
    assert provided.min_distance(TZ_CDRH3_10) == 0


def test_everything_provided_has_distance_zero(provided):
    facs = load_facs()
    rng = np.random.default_rng(0)
    for i in rng.choice(len(facs), 200, replace=False):
        assert provided.min_distance(str(facs.cdrh3[i])) == 0


def test_a_provided_sequence_is_rejected_by_the_gate(provided):
    facs = load_facs()
    # find one that is far enough from trastuzumab to pass the identity requirement,
    # but that is provided, and therefore should fail the distance requirement
    i = int(np.flatnonzero(facs.edit_distance >= 5)[0])
    loop = str(facs.cdrh3[i])
    v = check_design("provided", vh_from_cdrh3(loop), cdrh3=loop,
                     fmt="scFv", provided=provided)
    assert not v.passes
    assert any("provided" in r for r in v.reasons)


def test_distance_calculation_matches_raw_pass(provided):
    """The fast numpy path must give the same answer as the slow, explicit one."""
    facs = load_facs()
    alle = facs.cdrh3
    rng = np.random.default_rng(1)
    for _ in range(5):
        q = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 10))
        rask = provided.min_distance(q)
        qb = np.frombuffer(q.encode(), np.uint8)
        mat = np.frombuffer("".join(alle).encode(), np.uint8).reshape(len(alle), 10)
        treg = int((mat != qb).sum(axis=1).min())
        assert rask == min(treg, rask), (q, rask, treg)
        assert rask <= treg


# ── the metrics ─────────────────────────────────────────────────────────────
def test_auroc_known_values():
    assert auroc(np.array([1, 2, 3, 4]), np.array([0, 0, 1, 1])) == 1.0
    assert auroc(np.array([4, 3, 2, 1]), np.array([0, 0, 1, 1])) == 0.0
    assert auroc(np.array([1, 2, 3, 4]), np.array([0, 1, 0, 1])) == 0.75
    # equal values should give 0.5, not an arbitrary ordering
    assert auroc(np.array([1, 1, 1, 1]), np.array([0, 0, 1, 1])) == 0.5


def test_spearman_known_values():
    assert spearman(np.array([1, 2, 3]), np.array([1, 2, 3])) == pytest.approx(1.0)
    assert spearman(np.array([1, 2, 3]), np.array([3, 2, 1])) == pytest.approx(-1.0)
    # monotonic but not linear: Spearman should still be 1
    assert spearman(np.array([1, 2, 3, 4]),
                    np.array([1, 10, 100, 1000])) == pytest.approx(1.0)


def test_top_k():
    s = np.array([5.0, 4, 3, 2, 1])
    assert top_k_precision(s, np.array([1, 1, 0, 0, 0], bool), 2) == 1.0
    assert top_k_precision(s, np.array([0, 0, 0, 1, 1], bool), 2) == 0.0


# ── the data is as the README promises ──────────────────────────────────────
def test_facs_table_matches_readme():
    f = load_facs()
    assert len(f) == 524_346
    assert {str(x) for x in np.unique(f.binding_class)} == {"high", "medium", "low"}
    assert all(len(s) == 10 for s in f.cdrh3[:1000])
    assert f.edit_distance.min() >= 1 and f.edit_distance.max() <= 10


def test_facs_edit_distance_matches_our_own_calculation():
    f = load_facs()
    rng = np.random.default_rng(0)
    for i in rng.choice(len(f), 500, replace=False):
        assert hamming(str(f.cdrh3[i]), TZ_CDRH3_10) == int(f.edit_distance[i])


def test_facs_vh_contains_its_own_cdrh3():
    f = load_facs()
    rng = np.random.default_rng(0)
    for i in rng.choice(len(f), 300, replace=False):
        assert cdrh3_from_vh(str(f.vh[i])) == str(f.cdrh3[i])


def test_spr_table_matches_readme():
    s = load_spr(drop_other_antibody=False)
    assert len(s) == 768
    assert int(s.has_kd.sum()) == 190
    kd = s.kd_nM[s.has_kd]
    assert kd.min() == pytest.approx(0.023, abs=1e-3)
    assert kd.max() == pytest.approx(13863, rel=1e-3)


def test_trastuzumab_is_the_reference_in_spr():
    s = load_spr()
    tz = s.kd_nM[(s.cdrh3 == TZ_CDRH3_10) & s.has_kd]
    assert len(tz) == 8
    assert 0.02 < tz.min() < 0.19 and tz.max() < 0.19


def test_levenshtein():
    assert levenshtein("abc", "abc") == 0
    assert levenshtein("abc", "abd") == 1
    assert levenshtein("abc", "ab") == 1
    assert levenshtein("", "abc") == 3


# ── the fast novelty path must give the same answer as the slow one ─────────
def test_is_novel_matches_min_distance(provided):
    """The pigeonhole index is an optimisation. It must not change a single answer."""
    rng = np.random.default_rng(7)
    facs = load_facs()

    proever = []
    # random sequences (almost always novel)
    for _ in range(25):
        proever.append("".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 10)))
    # real provided ones (always distance 0)
    for i in rng.choice(len(facs), 15, replace=False):
        proever.append(str(facs.cdrh3[i]))
    # neighbours of provided ones, 1 and 2 mutations away
    for i in rng.choice(len(facs), 15, replace=False):
        s = list(str(facs.cdrh3[i]))
        for p in rng.choice(10, 2, replace=False):
            s[p] = str(rng.choice(list("ACDEFGHIKLMNPQRSTVWY")))
        proever.append("".join(s))

    for q in proever:
        rask = provided.is_novel(q)
        treg = provided.min_distance(q) >= 3
        assert rask == treg, (q, rask, treg, provided.min_distance(q))
