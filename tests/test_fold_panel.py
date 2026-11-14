"""Tests for the panel we fold and for the queries.json we build for the folding models.

What can go wrong without anyone noticing:
  - the panel is not deterministic, and two models get different sequences
  - trastuzumab is missing, and then none of the numbers mean anything
  - chain B is replaced, but chain C/D or residue_ids are lost
  - a VH of the wrong length slips through and does not fit the template
"""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fold_panel = load("fold_panel", "scripts/fold_panel.py")
fold_queries = load("fold_queries", "scripts/fold_queries.py")

# The panel is read from the SPR table. It does not ship with the code, so the
# tests that need it are skipped where it is missing.
needs_data = pytest.mark.skipif(
    not (ROOT / "data" / "trastuzumab_her2_spr_all_768.csv").exists(),
    reason="The SPR table is not here",
)

VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
    "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)


def reference() -> dict:
    return {
        "queries": [
            {
                "query_id": "example_0",
                "chains": [
                    {"chain_id": "A", "sequence": VL, "msa_path": None},
                    {"chain_id": "B", "sequence": VH, "msa_path": None},
                    {"chain_id": "C", "sequence": "CHQLC", "residue_ids": ["489", "490", "491", "492", "493"], "msa_path": None},
                    {"chain_id": "D", "sequence": "IWK", "residue_ids": ["591", "592", "593"], "msa_path": None},
                ],
            }
        ]
    }


def row(ident: str, cdrh3: str) -> dict:
    vh = VH.replace("WGGDGFYAMD", cdrh3)
    return {"id": ident, "group": "design", "cdrh3": cdrh3, "vh_sequence": vh, "vl_sequence": VL}


# ── panel ───────────────────────────────────────────────────────────────────
@needs_data
def test_panel_has_trastuzumab_first_and_is_deterministic():
    a = fold_panel.calibration_panel(12)
    b = fold_panel.calibration_panel(12)
    assert a == b, "two calls give different panels; then the models cannot be compared"
    assert a[0]["id"] == "trastuzumab" and a[0]["group"] == "baseline"
    assert a[0]["cdrh3"] == "WGGDGFYAMD"


@needs_data
def test_panel_spans_the_whole_affinity_range():
    panel = fold_panel.calibration_panel(12)
    kds = sorted(p["kd_nM"] for p in panel if p["group"] == "calibration")
    assert kds[0] < 1.0, "no strong binder in the calibration"
    assert kds[-1] > 100.0, "no weak binder in the calibration"
    assert len({p["cdrh3"] for p in panel}) == len(panel), "same CDRH3 twice"


@needs_data
def test_only_ten_residue_long_loops_are_included():
    # The template omits exactly ten residues. A shorter loop gives a VH that does not
    # fit, and would have mixed two geometries in the same comparison.
    for p in fold_panel.calibration_panel(400):
        assert len(p["cdrh3"]) == 10, p


@needs_data
def test_each_entry_has_its_own_vh_for_its_cdrh3():
    for p in fold_panel.calibration_panel(8):
        assert p["cdrh3"] in p["vh_sequence"]
        assert len(p["vh_sequence"]) == len(VH)


# ── queries.json ────────────────────────────────────────────────────────────
def test_only_chain_b_is_replaced():
    out = fold_queries.build(reference(), [row("d1", "APGDGFQIMD")])
    chains = {c["chain_id"]: c for c in out["queries"][0]["chains"]}
    assert set(chains) == {"A", "B", "C", "D"}
    assert chains["B"]["sequence"].count("APGDGFQIMD") == 1
    assert chains["C"]["residue_ids"][0] == "489"
    assert chains["D"]["sequence"] == "IWK"
    assert out["queries"][0]["query_id"] == "d1"


def test_reference_is_not_changed_between_entries():
    ref = reference()
    out = fold_queries.build(ref, [row("d1", "APGDGFQIMD"), row("d2", "MGPEPFYAMD")])
    assert ref["queries"][0]["chains"][1]["sequence"] == VH, "the reference was overwritten"
    seqs = [c["sequence"] for q in out["queries"] for c in q["chains"] if c["chain_id"] == "B"]
    assert seqs[0] != seqs[1]


def test_vh_with_wrong_length_is_stopped():
    bad = row("d1", "APGDGFQIMD")
    bad["vh_sequence"] = bad["vh_sequence"][:-3]
    # The message text gets translated when the code is shared, so we match on the numbers.
    with pytest.raises(ValueError, match=r"117.*120"):
        fold_queries.build(reference(), [bad])


def test_interface_iptm_takes_only_antibody_against_antigen():
    fold_analyse = load("fold_analyse", "scripts/fold_analyse.py")
    # A=0 B=1 (antibody), C=2 D=3 (HER2). Only the four cross pairs should count.
    m = [
        [1.0, 0.9, 0.5, 0.5],
        [0.9, 1.0, 0.5, 0.5],
        [0.5, 0.5, 1.0, 0.8],
        [0.5, 0.5, 0.8, 1.0],
    ]
    out = fold_analyse.interface_iptm({"chain_chain_interface_predicted_tm_score": m})
    assert out["iptm_ab_ag"] == pytest.approx(0.5), "VH-VL or HER2-internal pair included"
    assert out["iptm_vh_ag"] == pytest.approx(0.5)


def test_interface_iptm_is_asymmetric_mean():
    fold_analyse = load("fold_analyse", "scripts/fold_analyse.py")
    m = [
        [1.0, 1.0, 0.2, 0.2],
        [1.0, 1.0, 0.6, 0.4],
        [0.4, 0.8, 1.0, 1.0],
        [0.2, 0.6, 1.0, 1.0],
    ]
    out = fold_analyse.interface_iptm({"chain_chain_interface_predicted_tm_score": m})
    # VH(1) against C(2): (0.6+0.8)/2 = 0.7 ; VH against D(3): (0.4+0.6)/2 = 0.5
    assert out["iptm_vh_ag"] == pytest.approx(0.6)


def test_missing_chain_matrix_gives_no_interface_numbers():
    fold_analyse = load("fold_analyse", "scripts/fold_analyse.py")
    assert fold_analyse.interface_iptm({"predicted_tm_score": 0.9}) == {}


def test_reference_without_all_chains_is_stopped():
    ref = reference()
    ref["queries"][0]["chains"] = ref["queries"][0]["chains"][:2]
    with pytest.raises(ValueError, match=r"\['C', 'D'\]"):
        fold_queries.build(ref, [row("d1", "APGDGFQIMD")])
