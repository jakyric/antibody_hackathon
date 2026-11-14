"""Data-free sanity tests: python -m pytest tests  (or: python tests/test_basic.py)."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.common import VERNIER, WT_H3, WT_VH, WT_VL, build_fv, mutations, onehot  # noqa: E402
from pipeline.s1_enumerate import liabilities, substitution_table  # noqa: E402
from pipeline.s2_panel import vernier_options  # noqa: E402
from pipeline.s6_score import d0, ipsae, ilis  # noqa: E402
from pipeline.s5b_pose import kabsch, rmsd  # noqa: E402
from pipeline.s7_rank import apply_salt_bridge_rule  # noqa: E402
import pandas as pd  # noqa: E402


def test_mutations_and_fv():
    assert mutations("WGSQGFYAYD", WT_H3) == ["G3S", "D4Q", "M9Y"]
    vh, vl = build_fv("WGSQGFYAYD", {"L46": "I", "H94": "K"})
    assert vh[97] == "K" and vl[45] == "I" and "WGSQGFYAYD" in vh
    assert len(vh) == len(WT_VH) and len(vl) == len(WT_VL)
    for name, (chain, idx, wt, _) in VERNIER.items():          # Vernier numbering matches the sequences
        assert (WT_VH if chain == "H" else WT_VL)[idx - 1] == wt, name


def test_substitution_table_is_a_distribution_minus_stay():
    t = substitution_table("WGSQGFYAYD", [2, 4, 8, 10], 6, 0.25)
    assert [b for b, _ in t[2]][:3] == ["A", "D", "S"]           # G -> A, D, S (paper recipe)
    for p, subs in t.items():
        assert len(subs) == 6 and all(0 < v < 1 for _, v in subs) and sum(v for _, v in subs) < 1


def test_vernier_options_blosum45():
    o = vernier_options(1, ["H94", "L36"])
    assert set(o["H94"]) == {"K", "Q"} and set(o["L36"]) == {"F", "H", "W"}


def test_liabilities_wt():
    assert liabilities(WT_H3) == 3


def test_ipsae_ilis_extremes():
    n = 30
    good = np.full((n, n), 0.1)
    bad = np.full((n, n), 30.0)
    X, Y = np.arange(0, 15), np.arange(15, 30)
    xyz = np.random.default_rng(0).normal(size=(n, 3)) * 2
    assert ipsae(good, X, Y) > 0.9 and ipsae(bad, X, Y) == 0.0
    assert ilis(good, X, Y, xyz)[2] > 0.8 and ilis(bad, X, Y, xyz)[2] == 0.0
    assert math.isclose(d0(10), d0(27))                            # d0 is floored at 27 residues


def test_onehot_shape():
    assert onehot(["WGGDGFYAMD"]).sum() == 10


def test_salt_bridge_rule():
    bridge = pd.Series([True, True, False, False, False, False])
    p10 = pd.Series([False, True, True, False, True, False])       # position 10 changed
    h94 = pd.Series([False, False, True, False, False, True])      # H94 swapped
    r = apply_salt_bridge_rule(bridge, p10, h94).tolist()
    assert r == ["salt bridge present", "salt bridge present", "both partners mutated", "excluded", "excluded", "excluded"]


def test_kabsch_recovers_a_rigid_motion():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(40, 3))
    th = 0.7
    R0 = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Q = P @ R0.T + np.array([3.0, -2.0, 1.0])
    R, t = kabsch(P, Q)
    assert rmsd(P @ R.T + t, Q) < 1e-9


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v()
            print("ok", k)
