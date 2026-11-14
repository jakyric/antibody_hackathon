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
    rows = [row("best", 90, "WGGDGFYAMD"), row("naer", 80, "WGGDGFYAME"),
            row("annen", 70, "KKPLTSRIVN")]
    picked = ds.select(rows, n=2, epitope_floor=0.0)
    assert [r["design_id"] for r in picked] == ["best", "annen"]
    assert "cluster" in [r for r in rows if r["design_id"] == "naer"][0]["decision"]


def test_epitope_floor_stops_a_design_that_binds_somewhere_else():
    rows = [row("flyttet", 99, "AAAAAAAAAA", epi=0.2),
            row("holder", 40, "CCCCCCCCCC", epi=0.9)]
    picked = ds.select(rows, n=2, epitope_floor=0.8)
    assert [r["design_id"] for r in picked] == ["holder"]
    assert "epitope" in rows[0]["decision"]


def test_design_without_epitope_measurement_passes_the_floor():
    """If the measurement is missing, it must not silently count as passed
    OR discarded — it is ranked, and the fact that the number is missing is
    shown in the table."""
    rows = [row("umaalt", 50, "AAAAAAAAAA", epi="")]
    assert [r["design_id"] for r in ds.select(rows, n=1, epitope_floor=0.8)] == ["umaalt"]


def test_anchors_and_baseline_are_not_candidates():
    rows = [row("tz", 99, "WGGDGFYAMD", group="baseline"),
            row("design", 10, "CCCCCCCCCC")]
    assert [r["design_id"] for r in ds.select(rows, n=4, epitope_floor=0.0)] == ["design"]


def test_the_count_is_limited():
    # six completely different loops, so the clustering is not what limits the count
    loops = ["AAAAAAAAAA", "CCCCCCCCCC", "DDDDDDDDDD", "EEEEEEEEEE",
             "FFFFFFFFFF", "GGGGGGGGGG"]
    rows = [row(f"d{i}", i, loops[i]) for i in range(6)]
    assert len(ds.select(rows, n=2, epitope_floor=0.0)) == 2


def test_the_clustering_is_greedy_and_takes_the_highest_ranked_first():
    clusters = ds.cluster([("a", "WGGDGFYAMD"), ("b", "WGGDGFYAME"), ("c", "KKPLTSRIVN")])
    assert clusters["a"] == clusters["b"] != clusters["c"]


def test_a_design_that_refolds_far_from_the_designed_pose_is_dropped():
    """min_ipae can be good somewhere else entirely on HER2. Then it is not
    the binding we designed."""
    rows = [dict(row("flyttet", 99, "AAAAAAAAAA"), boltz2_rmsd_to_designed_pose="41.0"),
            dict(row("blir", 50, "CCCCCCCCCC"), boltz2_rmsd_to_designed_pose="2.7")]
    picked = ds.select(rows, n=2, epitope_floor=0.0, max_rmsd=10.0)
    assert [r["design_id"] for r in picked] == ["blir"]
    assert "designed pose" in rows[0]["decision"]


def test_without_an_rmsd_limit_all_are_included():
    rows = [dict(row("a", 10, "AAAAAAAAAA"), boltz2_rmsd_to_designed_pose="41.0")]
    assert len(ds.select(rows, n=1, epitope_floor=0.0)) == 1


def test_a_missing_rmsd_does_not_stop_the_design():
    rows = [dict(row("a", 10, "AAAAAAAAAA"), boltz2_rmsd_to_designed_pose="")]
    assert len(ds.select(rows, n=1, epitope_floor=0.0, max_rmsd=10.0)) == 1


def test_a_vhh_does_not_get_an_scfv_linker(tmp_path):
    """An empty light chain means a VHH. If the linker tags along, the sequence is wrong."""
    import importlib.util as iu
    spec = iu.spec_from_file_location("denovo_table", ROOT / "scripts/denovo_table.py")
    dt = iu.module_from_spec(spec)
    sys.modules["denovo_table"] = dt
    spec.loader.exec_module(dt)

    cands = tmp_path / "c.csv"
    cands.write_text("design_id,vh_sequence,vl_sequence,cdrh3,weight,seed,run_id\n"
                     "d1,EVQLVESGGG,,AAA,antibody-2026-10-30,11,r1\n")
    out = tmp_path / "t.csv"
    dt.main(["--designs", str(cands), "--out", str(out)])
    import csv as _csv
    row = next(iter(_csv.DictReader(out.open())))
    assert row["kind"] == "vhh"
    assert row["heavy_sequence"] == "EVQLVESGGG"
    assert "GGGGS" not in row["heavy_sequence"]
