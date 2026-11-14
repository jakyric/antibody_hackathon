"""Tests for the selection.

The rule decides what is actually delivered to the lab. It should be readable
as code and testable as code: ranking, epitope floor and spread each on their own.
"""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ds = load("denovo_select", "scripts/denovo_select.py")


def row(ident, pct, cdrh3, epi=1.0, group="denovo"):
    return {"design_id": ident, "group": group, "cdrh3": cdrh3,
            "binding_percentile_mean": str(pct), "epitope_jaccard": str(epi)}


def test_order_follows_the_binding_percentile():
    rows = [row("a", 10, "AAAAAAAAAA"), row("b", 90, "CCCCCCCCCC"),
            row("c", 50, "DDDDDDDDDD")]
    picked = ds.select(rows, n=3, epitope_floor=0.0)
    assert [r["design_id"] for r in picked] == ["b", "c", "a"]


def test_near_cdrh3_counts_as_the_same_idea():
    """Two loops one mutation apart are the same idea, not two proposals."""
    rows = [row("best", 90, "WGGDGFYAMD"), row("near", 80, "WGGDGFYAME"),
            row("other", 70, "KKPLTSRIVN")]
    picked = ds.select(rows, n=2, epitope_floor=0.0)
    assert [r["design_id"] for r in picked] == ["best", "other"]
    assert "cluster" in [r for r in rows if r["design_id"] == "near"][0]["decision"]


def test_epitope_floor_stops_designs_that_bind_elsewhere():
    rows = [row("moved", 99, "AAAAAAAAAA", epi=0.2),
            row("holds", 40, "CCCCCCCCCC", epi=0.9)]
    picked = ds.select(rows, n=2, epitope_floor=0.8)
    assert [r["design_id"] for r in picked] == ["holds"]
    assert "epitope" in rows[0]["decision"]


def test_design_without_epitope_measurement_passes_the_floor():
    """If the measurement is missing, it must not be silently counted as passed
    OR discarded — it is ranked, and the fact that the number is missing is
    shown in the table."""
    rows = [row("unmeasured", 50, "AAAAAAAAAA", epi="")]
    assert [r["design_id"] for r in ds.select(rows, n=1, epitope_floor=0.8)] == ["unmeasured"]


def test_anchors_and_baseline_are_not_candidates():
    rows = [row("tz", 99, "WGGDGFYAMD", group="baseline"),
            row("design", 10, "CCCCCCCCCC")]
    assert [r["design_id"] for r in ds.select(rows, n=4, epitope_floor=0.0)] == ["design"]


def test_count_is_limited():
    # six completely different loops, so the clustering is not what limits the count
    loops = ["AAAAAAAAAA", "CCCCCCCCCC", "DDDDDDDDDD", "EEEEEEEEEE",
             "FFFFFFFFFF", "GGGGGGGGGG"]
    rows = [row(f"d{i}", i, loops[i]) for i in range(6)]
    assert len(ds.select(rows, n=2, epitope_floor=0.0)) == 2


def test_clustering_is_greedy_and_takes_the_highest_ranked_first():
    clusters = ds.cluster([("a", "WGGDGFYAMD"), ("b", "WGGDGFYAME"), ("c", "KKPLTSRIVN")])
    assert clusters["a"] == clusters["b"] != clusters["c"]
