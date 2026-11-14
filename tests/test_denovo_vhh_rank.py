"""Tests for the ranking of the VHH designs.

The easy thing to get backwards: min_ipae is BETTER when it is low, while
the percentile should be high for what is good. And with several seeds per design,
the row that is kept must be the best one, not the last one that was read.
"""

import importlib.util
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
    pct = vr.percentiles({"god": 5.0, "midt": 8.0, "daarlig": 12.0}, higher_is_better=False)
    assert pct["god"] == 100.0
    assert pct["daarlig"] == 0.0
    assert 0 < pct["midt"] < 100


def test_empty_values_are_not_counted():
    pct = vr.percentiles({"a": 5.0, "b": float("nan")}, higher_is_better=False)
    assert list(pct) == ["a"]


def test_a_single_design_alone_gets_a_value_and_does_not_crash():
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


def test_invalid_numbers_become_nan_and_do_not_stop_the_reading():
    assert vr.num("") != vr.num("")      # nan is not equal to itself
    assert vr.num("3.5") == 3.5
