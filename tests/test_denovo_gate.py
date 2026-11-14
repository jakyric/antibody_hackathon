"""Tests for the de novo assessment.

The whole point of this script is to separate two readings of the same
requirement: identity over the whole sequence versus identity over the CDRs
only. If they are mixed up, we answer the wrong question about which track
the designs belong in.
"""

import importlib.util
import sys

import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dg = load("denovo_gate", "scripts/denovo_gate.py")

from abx.sequences import TRASTUZUMAB_VH, TRASTUZUMAB_VL  # noqa: E402


def swap(seq: str, spans, filler: str = "A") -> str:
    """Replace each CDR with the same number of alanines. The framework stays."""
    out = list(seq)
    for a, b in spans:
        for i in range(a - 1, b):
            out[i] = filler
    return "".join(out)


def test_chain_string_is_split_correctly():
    chains = dg.parse_sequences("A:DIQ;B:EVQ;C:CHQ;D:IWK")
    assert chains == {"A": "DIQ", "B": "EVQ", "C": "CHQ", "D": "IWK"}


def test_empty_string_gives_no_chains():
    assert dg.parse_sequences("") == {}


def test_trastuzumab_against_itself_is_one_hundred_percent():
    r = dg.evaluate(TRASTUZUMAB_VL, TRASTUZUMAB_VH, None)
    assert r["identity_whole_pct"] == 100.0
    assert r["identity_cdr_pct"] == 100.0
    assert r["cdrh3"] == "WGGDGFYAMD"
    assert r["unpaired_cys"] == 0


def test_six_new_cdrs_give_low_cdr_identity_but_high_overall():
    """This is the question itself, and the numbers should show why it is open."""
    vl = swap(TRASTUZUMAB_VL, dg.CDR_LIGHT)
    vh = swap(TRASTUZUMAB_VH, dg.CDR_HEAVY)
    r = dg.evaluate(vl, vh, None)
    assert r["identity_cdr_pct"] < 30, "The CDRs have been replaced, so the identity should drop"
    assert r["identity_whole_pct"] > 70, (
        "the framework remains, so the whole sequence is still similar to trastuzumab "
        "— that is exactly what makes the reading of the requirement decisive"
    )


def test_odd_number_of_cysteines_is_caught():
    vh = TRASTUZUMAB_VH[:40] + "C" + TRASTUZUMAB_VH[41:]
    r = dg.evaluate(TRASTUZUMAB_VL, vh, None)
    assert r["unpaired_cys"] == 1


@pytest.mark.skipif(not (ROOT / "configs" / "her2_all6cdr.yaml").exists(),
                    reason="the design specification does not ship with the code")
def test_cdr_boundaries_are_the_same_as_in_the_specification():
    """The boundaries live in two places. If they drift apart, we design something
    other than what we measure."""
    spec = (ROOT / "configs" / "her2_all6cdr.yaml").read_text()
    # the specification lists the framework, i.e. what is NOT CDR
    assert "[1-23, 35-49, 57-88, 98-107]" in spec
    assert "[1-30, 36-49, 67-98, 109-120]" in spec
    assert dg.CDR_LIGHT == [(24, 34), (50, 56), (89, 97)]
    assert dg.CDR_HEAVY == [(31, 35), (50, 66), (99, 108)]
