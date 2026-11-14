"""Tests for the collection and filtering of the de novo designs.

The filtering decides what gets GPU time going forward, and "known
binder" is the mistake that costs the most: a design that is already in
the FACS screen never gets made in the lab. So each reason for failing
is tested on its own.
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


def test_lone_cysteine_fails():
    ok, why = dc.verdict(base_row(unpaired_cys=1))
    assert not ok and "cysteine" in why


def test_too_long_scfv_fails():
    ok, why = dc.verdict(base_row(scfv_len=251))
    assert not ok and "250" in why


def test_limit_of_250_is_included():
    assert dc.verdict(base_row(scfv_len=250))[0]


def test_seventy_percent_cdr_identity_is_not_below_seventy():
    assert not dc.verdict(base_row(identity_cdr_pct=70.0))[0]
    assert dc.verdict(base_row(identity_cdr_pct=69.9))[0]


def test_known_binder_fails_at_distance_two():
    """The de novo requirement is 2; we require 3, for a margin against the lab's base."""
    ok, why = dc.verdict(base_row(cdrh3_min_dist_provided=2))
    assert not ok and "known binder" in why
    assert dc.verdict(base_row(cdrh3_min_dist_provided=3))[0]


def test_missing_cdrh3_fails():
    assert not dc.verdict(base_row(cdrh3=""))[0]


def test_empty_fields_fail_and_do_not_crash():
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


def test_chain_string_is_split_into_light_and_heavy_chain():
    assert dc.parse_chains("A:DIQ;B:EVQ;C:X;D:Y")["B"] == "EVQ"


def test_sequence_is_read_from_a_generated_structure(tmp_path):
    """The fallback route when a run is interrupted: the structure file is the ground truth."""
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


def test_vhh_cdr3_is_taken_between_the_anchors():
    r = dc.evaluate_vhh(CAPLA_VHH)
    assert r["cdrh3"] == "ARYGGSGPGIFESSYDA"
    assert r["vl_len"] == 0 and r["scfv_len"] == len(CAPLA_VHH)


def test_vhh_is_far_below_seventy_percent_against_trastuzumab():
    """This is the whole reason the VHH track exists: the official checker
    measures identity against trastuzumab over the whole sequence."""
    r = dc.evaluate_vhh(CAPLA_VHH)
    assert r["identity_whole_pct"] < 70
    assert dc.verdict_vhh(r)[0]


def test_vhh_on_trastuzumabs_own_vh_fails():
    from abx.sequences import TRASTUZUMAB_VH
    r = dc.evaluate_vhh(TRASTUZUMAB_VH)
    ok, why = dc.verdict_vhh(r)
    assert not ok and "70" in why


def test_vhh_with_lone_cysteine_fails():
    r = dc.evaluate_vhh(CAPLA_VHH + "C")
    assert not dc.verdict_vhh(r)[0]


def test_subset_filter_matches_exactly(tmp_path):
    """"her2_vhh_denovo_3" must not let in "..._39"."""
    import importlib.util as iu
    spec = iu.spec_from_file_location("subset_queries", ROOT / "scripts/subset_queries.py")
    sq = iu.module_from_spec(spec)
    sys.modules["subset_queries"] = sq
    spec.loader.exec_module(sq)

    q = tmp_path / "q.json"
    q.write_text(json.dumps({"queries": [{"query_id": "d_3"}, {"query_id": "d_39"}]}))
    keep = tmp_path / "keep.csv"
    keep.write_text("design_id\n20261114T115955-d_3\n")
    out = tmp_path / "out.json"
    sq.main(["--queries", str(q), "--keep", str(keep), "--out", str(out)])
    got = json.loads(out.read_text())["queries"]
    assert [x["query_id"] for x in got] == ["d_3"]
