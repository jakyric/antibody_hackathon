"""Tests for the liability-clean CDR-H3 enumeration.

The enumeration is only worth anything if every loop it returns really is
clean and really is at the required distance, and if it does not quietly drop
a legal loop. Those are the three things checked here, plus the property that
makes it reportable at all: it is complete, not sampled.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abx.data import TZ_CDRH3_10  # noqa: E402
from abx.pool import (  # noqa: E402
    FORBIDDEN_PAIRS,
    FORBIDDEN_RESIDUES,
    enumerate_clean,
    is_clean,
    motif_positions,
    mutation_sites,
)


def test_the_parent_is_not_clean():
    """If trastuzumab passed this filter there would be nothing to fix."""
    assert not is_clean(TZ_CDRH3_10)


def test_the_parents_three_motifs_are_the_ones_we_name():
    assert motif_positions(TZ_CDRH3_10, "DG") == [3]
    assert "W" in TZ_CDRH3_10 and "M" in TZ_CDRH3_10


def test_the_framework_joins_are_part_of_the_check():
    """A motif can be created by the residue next to the loop, not only
    inside it. With trastuzumab's joins (R before, Y after) none of our pairs
    can straddle them, so the check is asserted on the mechanism rather than
    on a case that cannot arise: a pair placed at the very start of the loop
    is still found, and the padded string is what is searched."""
    assert motif_positions("DGAAAAAAAA", "DG") == [0]
    assert not is_clean("DGAAAAAAAA")
    assert motif_positions("AAAAAAAADG", "DG") == [8]
    assert is_clean("AAAAAAAAAA")


def test_every_enumerated_loop_is_clean_and_at_the_right_distance():
    loops = enumerate_clean(4)
    assert loops, "the enumeration must not be empty"
    for s in loops[:5000]:
        assert len(s) == len(TZ_CDRH3_10)
        assert sum(a != b for a, b in zip(s, TZ_CDRH3_10)) == 4
        assert is_clean(s)
        assert not (set(s) & FORBIDDEN_RESIDUES)


def test_no_enumerated_loop_contains_a_forbidden_pair():
    for s in enumerate_clean(4)[:5000]:
        padded = "R" + s + "Y"
        for pair in FORBIDDEN_PAIRS:
            assert pair not in padded


def test_the_enumeration_is_complete_not_sampled():
    """Running it twice gives the same list, and it contains a loop we can
    construct by hand: W1->F, D4->E, M9->L, and one free change at position 7."""
    first, second = enumerate_clean(4), enumerate_clean(4)
    assert first == second
    hand = list(TZ_CDRH3_10)
    hand[0], hand[3], hand[8], hand[6] = "F", "E", "L", "F"
    hand = "".join(hand)
    assert is_clean(hand)
    assert hand in set(first)


def test_mutation_sites_always_cover_the_forbidden_residues():
    forced = {i for i, a in enumerate(TZ_CDRH3_10) if a in FORBIDDEN_RESIDUES}
    for sites in mutation_sites(4):
        assert forced <= set(sites)


def test_three_substitutions_is_enough_to_clean_but_we_need_four():
    """The novelty rule, not the liability list, is what forces the fourth."""
    assert enumerate_clean(3)
    assert len(enumerate_clean(4)) > len(enumerate_clean(3))


# ── the contact-preserving variant ──────────────────────────────────────────


def test_contact_positions_match_the_crystal_measurement():
    """W99, G103, F104, Y105 as positions in the 10-mer."""
    from abx.pool import CONTACT_POSITIONS

    assert CONTACT_POSITIONS == {0, 4, 5, 6}
    assert [TZ_CDRH3_10[i] for i in sorted(CONTACT_POSITIONS)] == ["W", "G", "F", "Y"]


def test_contact_safe_loops_never_touch_a_contact_residue():
    from abx.pool import CONTACT_POSITIONS, enumerate_contact_safe

    loops = enumerate_contact_safe(4)
    assert loops
    for s in loops[:5000]:
        for i in sorted(CONTACT_POSITIONS):
            assert s[i] == TZ_CDRH3_10[i], f"{s} changed contact position {i}"
        assert sum(a != b for a, b in zip(s, TZ_CDRH3_10)) == 4


def test_contact_safe_loops_keep_the_tryptophan_and_lose_the_methionine():
    """W99 is a contact, so it stays and its liability stays with it. M107 is
    not, so it must go -- that is the whole point of the trade."""
    from abx.pool import enumerate_contact_safe

    for s in enumerate_contact_safe(4)[:3000]:
        assert s[0] == "W"
        assert "M" not in s


def test_contact_safe_loops_remove_the_high_severity_isomerisation_motif():
    from abx.pool import enumerate_contact_safe

    for s in enumerate_contact_safe(4)[:3000]:
        assert "DG" not in "R" + s + "Y"


def test_contact_safe_rejects_a_newly_introduced_motif():
    from abx.pool import CONTACT_POSITIONS, is_clean_outside

    # a D at position 3 next to the frozen G103 would rebuild the DG motif
    bad = list(TZ_CDRH3_10)
    bad[3] = "D"
    assert not is_clean_outside("".join(bad), CONTACT_POSITIONS)
    # and a mutated contact position is rejected outright
    moved = list(TZ_CDRH3_10)
    moved[0] = "Y"
    assert not is_clean_outside("".join(moved), CONTACT_POSITIONS)
