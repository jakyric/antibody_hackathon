"""The winner-hidden test: the pool, the leakage and the rankings."""

from __future__ import annotations

import numpy as np
import pytest

from abx.data import TZ_CDRH3_10
from abx.retro import (
    capture,
    edit_distances,
    enrichment,
    leakage_mask,
    legal_pool,
    winner_mask,
    winner_ranks,
)


def test_edit_distances_counts_positions() -> None:
    assert edit_distances([TZ_CDRH3_10])[0] == 0
    assert edit_distances(["AGGDGFYAMD"])[0] == 1
    assert edit_distances(["AAGDGFYAMD"])[0] == 2


# ── the pool ────────────────────────────────────────────────────────────────
def test_legal_pool_keeps_only_quantified_far_enough_variants() -> None:
    seqs = [TZ_CDRH3_10, "AAADGFYAMD", "AAGDGFYAMD", "AAAAGFYAMD"]
    kd = [0.1, 2.0, 3.0, 4.0]
    res = ["quantified_binding"] * 4
    pool, vals = legal_pool(seqs, kd, res)
    # trastuzumab (0) and the 2-mutant drop out; the 3- and 4-mutants are kept.
    assert pool == ["AAAAGFYAMD", "AAADGFYAMD"]
    assert list(vals) == [4.0, 2.0]


def test_legal_pool_drops_rows_without_a_quantified_kd() -> None:
    pool, _ = legal_pool(["AAADGFYAMD"], [2.0], ["no_reportable_kd"])
    assert pool == []


def test_legal_pool_keeps_the_best_measurement_of_a_repeat() -> None:
    pool, vals = legal_pool(
        ["AAADGFYAMD", "AAADGFYAMD"], [5.0, 2.0], ["quantified_binding"] * 2
    )
    assert pool == ["AAADGFYAMD"]
    assert vals[0] == pytest.approx(2.0)


def test_legal_pool_ignores_sequences_of_another_length() -> None:
    pool, _ = legal_pool(["ARDRIVDDYYYYGMDV"], [1.0], ["quantified_binding"])
    assert pool == []


def test_winner_mask_is_strictly_below_the_threshold() -> None:
    assert list(winner_mask(np.array([0.5, 1.0, 2.0]))) == [True, False, False]


# ── leakage ─────────────────────────────────────────────────────────────────
def test_leakage_mask_removes_the_winner_and_its_neighbours() -> None:
    winner = "AAADGFYAMD"
    train = [winner, "AAADGFYAMA", "AAADGFYAAA", TZ_CDRH3_10]
    #         distance 0, 1,          2,          3
    mask = leakage_mask(train, [winner], radius=1)
    assert list(mask) == [True, True, False, False]


def test_leakage_mask_radius_zero_removes_only_the_winner() -> None:
    winner = "AAADGFYAMD"
    mask = leakage_mask([winner, "AAADGFYAMA"], [winner], radius=0)
    assert list(mask) == [True, False]


def test_leakage_mask_without_hidden_sequences_removes_nothing() -> None:
    assert not leakage_mask([TZ_CDRH3_10], []).any()


# ── rankings ────────────────────────────────────────────────────────────────
def test_winner_ranks_puts_the_best_score_first() -> None:
    scores = np.array([10.0, 5.0, 1.0])
    assert list(winner_ranks(scores, np.array([True, False, False]))) == [1]
    assert list(winner_ranks(scores, np.array([False, False, True]))) == [3]


def test_winner_ranks_is_pessimistic_about_ties() -> None:
    """A scorer that gives everything the same number has found nothing at all."""
    scores = np.zeros(4)
    assert list(winner_ranks(scores, np.array([True, False, False, False]))) == [4]


def test_capture_compares_with_chance() -> None:
    # 10 in the pool, 2 winners, both at the top.
    scores = np.arange(10, dtype=float)[::-1]
    is_w = np.zeros(10, dtype=bool)
    is_w[:2] = True
    hit, exp = capture(scores, is_w, k=2)
    assert hit == 2
    assert exp == pytest.approx(0.4)  # 2 * 2 / 10
    assert enrichment(scores, is_w, k=2) == pytest.approx(5.0)


def test_enrichment_of_one_means_no_information() -> None:
    scores = np.arange(10, dtype=float)[::-1]
    is_w = np.zeros(10, dtype=bool)
    is_w[:1] = True
    # one winner, top 10 of 10: captured 1, expected 1.
    assert enrichment(scores, is_w, k=10) == pytest.approx(1.0)
