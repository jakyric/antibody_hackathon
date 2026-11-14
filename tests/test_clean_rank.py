"""Tests for the rank aggregation used on the liability-clean picks."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from clean_rank import average_ranks  # noqa: E402


def test_higher_is_better_puts_the_largest_first():
    assert average_ranks([0.9, 1.0, 0.5], higher_is_better=True) == [2.0, 1.0, 3.0]


def test_ties_share_the_average_of_the_places_they_fill():
    """Three equal values occupy places 1, 2 and 3, so each ranks 2."""
    assert average_ranks([1.0, 1.0, 1.0, 0.5], higher_is_better=True) == [
        2.0, 2.0, 2.0, 4.0
    ]


def test_a_tie_does_not_depend_on_the_order_the_rows_arrive_in():
    a = average_ranks([1.0, 0.95, 1.0], higher_is_better=True)
    b = average_ranks([1.0, 1.0, 0.95], higher_is_better=True)
    assert a == [1.5, 3.0, 1.5]
    assert b == [1.5, 1.5, 3.0]


def test_ranks_always_sum_to_the_same_total_whatever_the_ties():
    """Averaging must not create or destroy rank mass, or the aggregation
    would quietly favour the axis with more ties."""
    for values in ([3, 2, 1], [1, 1, 1], [2, 2, 1]):
        assert sum(average_ranks(values, higher_is_better=True)) == 6
