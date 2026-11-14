"""Tests for the one-command requirement check.

The point of the check is that it is trustworthy, so the tests are mostly
about it REJECTING things. A check that only ever passes is indistinguishable
from no check at all.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from requirement_check import FAIL, PASS, UNKNOWN, check_row  # noqa: E402

TZ_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
         "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
TZ_VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
         "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
PROVIDED = ["WGGDGFYAMD", "MGPEPFYAMD", "AAADGFYAMD"]


def row(sid="d", vh=TZ_VH, vl=TZ_VL, kind="paired"):
    return {"sequence_id": sid, "kind": kind,
            "heavy_sequence": vh, "light_sequence": vl}


def swap_loop(vh: str, loop: str) -> str:
    start = vh.index("YYCSR") + 5
    end = vh.index("YWGQGT", start)
    return vh[:start] + loop + vh[end:]


def test_trastuzumab_itself_must_fail_the_novelty_rule():
    """The control. A check that lets trastuzumab through is broken."""
    res = check_row(row("trastuzumab"), PROVIDED)
    assert res["checks"]["cdrh3_identity_under_70pct"]["verdict"] == FAIL
    assert res["checks"]["mutations_from_provided_data"]["verdict"] == FAIL


def test_a_legal_design_passes_every_checkable_requirement():
    res = check_row(row("ok", vh=swap_loop(TZ_VH, "FGGWPWYAMD")), PROVIDED)
    c = res["checks"]
    assert c["format"]["verdict"] == PASS
    assert c["length_max_250"]["verdict"] == PASS
    assert c["no_unpaired_cysteines"]["verdict"] == PASS
    assert c["cdrh3_identity_under_70pct"]["verdict"] == PASS
    assert c["mutations_from_provided_data"]["verdict"] == PASS


def test_a_design_one_mutation_from_provided_data_fails():
    # one away from MGPEPFYAMD, which is in PROVIDED
    res = check_row(row("near", vh=swap_loop(TZ_VH, "MGPEPFYAMA")), PROVIDED)
    c = res["checks"]["mutations_from_provided_data"]
    assert c["verdict"] == FAIL
    assert "MGPEPFYAMD" in c["detail"]


def test_margin_is_reported_as_a_number_not_a_tick():
    res = check_row(row("ok", vh=swap_loop(TZ_VH, "FGGWPWYAMD")), PROVIDED)
    assert "margin" in res["checks"]["mutations_from_provided_data"]["margin"]


def test_odd_cysteines_in_one_chain_are_caught_even_when_the_total_is_even():
    """The per-domain rule.

    One extra cysteine in each chain makes the molecule-wide count even. A
    check that counts across the whole construct would call this clean.
    """
    vh = swap_loop(TZ_VH, "FGGWPWYAMD") + "C"
    vl = TZ_VL + "C"
    res = check_row(row("twoC", vh=vh, vl=vl), PROVIDED)
    assert res["checks"]["no_unpaired_cysteines"]["verdict"] == FAIL


def test_a_bad_format_fails():
    res = check_row(row("x", vh=swap_loop(TZ_VH, "FGGWPWYAMD"), kind="nanobody-ish"),
                    PROVIDED)
    assert res["checks"]["format"]["verdict"] == FAIL


def test_over_length_fails():
    long_vh = swap_loop(TZ_VH, "FGGWPWYAMD") + "G" * 200
    res = check_row(row("long", vh=long_vh), PROVIDED)
    assert res["checks"]["length_max_250"]["verdict"] == FAIL


def test_the_two_organiser_only_checks_are_unknown_never_pass():
    res = check_row(row("ok", vh=swap_loop(TZ_VH, "FGGWPWYAMD")), PROVIDED)
    assert res["checks"]["immunum_numbering_confidence"]["verdict"] == UNKNOWN
    assert res["checks"]["identity_to_known_antibodies"]["verdict"] == UNKNOWN
