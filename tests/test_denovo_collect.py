"""Tests for the collection and filtering of the de novo designs.

Filtering decides what gets to use GPU time further on, and "known
binder" is the mistake that costs the most: a design that is already in
the FACS screen is never made in the lab. That is why each reason for
discarding a design is tested separately.
"""

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


dc = load("denovo_collect", "scripts/denovo_collect.py")


def base_row(**over):
    row = {
        "unpaired_cys": 0,
        "scfv_len": 242,
        "identity_cdr_pct": 55.0,
        "cdrh3": "ALKDGWYPQE",
        "cdrh3_min_dist_provided": 5,
    }
    row.update(over)
    return row


def test_a_clean_design_passes():
    ok, why = dc.verdict(base_row())
    assert ok and why == ""


def test_unpaired_cysteine_is_discarded():
    ok, why = dc.verdict(base_row(unpaired_cys=1))
    assert not ok and "cysteine" in why


def test_too_long_scfv_is_discarded():
    ok, why = dc.verdict(base_row(scfv_len=251))
    assert not ok and "250" in why


def test_the_limit_at_250_is_included():
    assert dc.verdict(base_row(scfv_len=250))[0]


def test_seventy_percent_cdr_identity_is_not_below_seventy():
    assert not dc.verdict(base_row(identity_cdr_pct=70.0))[0]
    assert dc.verdict(base_row(identity_cdr_pct=69.9))[0]


def test_known_binder_is_discarded_at_distance_two():
    """The de novo requirement is 2; we require 3, for a margin against the lab's base."""
    ok, why = dc.verdict(base_row(cdrh3_min_dist_provided=2))
    assert not ok and "known binder" in why
    assert dc.verdict(base_row(cdrh3_min_dist_provided=3))[0]


def test_missing_cdrh3_is_discarded():
    assert not dc.verdict(base_row(cdrh3=""))[0]


def test_empty_fields_are_discarded_and_do_not_crash():
    for key in ("scfv_len", "identity_cdr_pct", "cdrh3_min_dist_provided"):
        assert not dc.verdict(base_row(**{key: ""}))[0]


def test_weight_and_seed_are_read_from_run_json(tmp_path):
    run = tmp_path / "20261114T113615Z-denovo-antibody-2026-10-30-s1"
    run.mkdir()
    (run / "run.json").write_text(json.dumps(
        {"name": "denovo-antibody-2026-10-30-s1", "seed": 1,
         "git_sha": "abcdef1234567890"}))
    meta = dc.run_metadata(run)
    assert meta["weight"] == "antibody-2026-10-30"
    assert meta["seed"] == "1"
    assert meta["run_id"] == run.name
    assert meta["commit"] == "abcdef123456"


def test_the_chain_string_is_split_into_light_and_heavy_chain():
    assert dc.parse_chains("A:DIQ;B:EVQ;C:X;D:Y")["B"] == "EVQ"


def test_sequence_is_read_from_a_generated_structure(tmp_path):
    """The fallback when a run is interrupted: the structure file is the ground truth."""
    cif = tmp_path / "d_generated.cif"
    cif.write_text(
        "loop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"
        "ATOM 1 N ASP 1 A\n"
        "ATOM 2 CA ASP 1 A\n"
        "ATOM 3 N ILE 2 A\n"
        "ATOM 4 N GLU 1 B\n"
        "#\n"
    )
    assert dc.sequences_from_cif(cif) == "A:DI;B:E"


# ── VHH track ───────────────────────────────────────────────────────────────

CAPLA_VHH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGRSTFGAYGMGWFRQAPGKGRELVAAISGSGLNTSYPDSVEGRF"
    "TISRDNAKRMVYLQMNSLRAEDTAVYYCARYGGSGPGIFESSYDAWGQGTQVTVS"
)


def test_vhh_cdr3_is_extracted_between_the_anchors():
    r = dc.evaluate_vhh(CAPLA_VHH)
    assert r["cdrh3"] == "ARYGGSGPGIFESSYDA"
    assert r["vl_len"] == 0 and r["scfv_len"] == len(CAPLA_VHH)


def test_vhh_is_far_below_seventy_percent_against_trastuzumab():
    """This is the whole reason the VHH track exists: the official checker
    measures the identity against trastuzumab over the whole sequence."""
    r = dc.evaluate_vhh(CAPLA_VHH)
    assert r["identity_whole_pct"] < 70
    assert dc.verdict_vhh(r)[0]


def test_vhh_on_trastuzumabs_own_vh_is_discarded():
    from abx.sequences import TRASTUZUMAB_VH
    r = dc.evaluate_vhh(TRASTUZUMAB_VH)
    ok, why = dc.verdict_vhh(r)
    assert not ok and "70" in why


def test_vhh_with_unpaired_cysteine_is_discarded():
    r = dc.evaluate_vhh(CAPLA_VHH + "C")
    assert not dc.verdict_vhh(r)[0]
