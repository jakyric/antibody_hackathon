"""Tests for the ranking of the VHH designs.

What is easy to get backwards: min_ipae is BETTER when it is low, while the
percentile should be high for what is good. And with several seeds per design,
the row that is kept must be the best one, not the last one read.
"""

import importlib.util
import pathlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


vr = load("denovo_vhh_rank", "scripts/denovo_vhh_rank.py")


def test_low_min_ipae_gives_high_percentile():
    pct = vr.percentiles({"good": 5.0, "middle": 8.0, "poor": 12.0}, higher_is_better=False)
    assert pct["good"] == 100.0
    assert pct["poor"] == 0.0
    assert 0 < pct["middle"] < 100


def test_empty_values_are_not_counted():
    pct = vr.percentiles({"a": 5.0, "b": float("nan")}, higher_is_better=False)
    assert list(pct) == ["a"]


def test_single_design_alone_gets_a_value_and_does_not_crash():
    assert vr.percentiles({"a": 1.0}, higher_is_better=False) == {"a": 100.0}


def test_best_seed_per_design_is_kept(tmp_path):
    csv_path = tmp_path / "refold_metrics.csv"
    csv_path.write_text(
        "query_id,seed,min_ipae,iptm\n"
        "d1,42,11.0,0.5\n"
        "d1,43,7.0,0.8\n"
        "d2,42,9.0,0.6\n"
    )
    rows = vr.read_metrics(csv_path)
    assert rows["d1"]["min_ipae"] == "7.0"
    assert rows["d1"]["_n"] == 2
    assert rows["d2"]["_n"] == 1


def write_metrics(path, qid, ipae):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"query_id,seed,min_ipae,iptm\n{qid},42,{ipae},0.5\n")


def test_designs_from_different_runs_do_not_merge(tmp_path):
    """Each run numbers from zero, so her2_vhh_denovo_0 exists in all of them.
    Without the run in the key, two different designs become one."""
    a = tmp_path / "runs/20261114T115955Z-denovo-vhh-s11/refold/scored-boltz2/refold_metrics.csv"
    b = tmp_path / "runs/20261114T120004Z-denovo-vhh-s12/refold/scored-boltz2/refold_metrics.csv"
    write_metrics(a, "her2_vhh_denovo_0", 3.0)
    write_metrics(b, "her2_vhh_denovo_0", 9.0)
    merged = vr.read_all([a, b])
    assert len(merged) == 2, "two different designs should be counted as two"
    assert set(merged) == {"20261114T115955-her2_vhh_denovo_0",
                           "20261114T120004-her2_vhh_denovo_0"}


def test_run_id_is_extracted_from_the_path(tmp_path):
    p = tmp_path / "runs/20261114T115955Z-denovo-vhh-s11/refold/scored-boltz2/m.csv"
    assert vr.run_prefix(p) == "20261114T115955"
    assert vr.run_prefix(pathlib.Path("/tmp/without/run/m.csv")) == ""
    # the epitope files have the run ID inside the file name, not as a separate folder
    assert vr.run_prefix(pathlib.Path(
        "results/denovo/vhh/epitope_20261114T115955Z-denovo-vhh-s11_boltz2.csv"
    )) == "20261114T115955"


def test_invalid_numbers_become_nan_and_do_not_stop_the_reading():
    assert vr.num("") != vr.num("")      # nan is not equal to itself
    assert vr.num("3.5") == 3.5
