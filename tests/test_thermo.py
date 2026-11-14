"""Tests for the disulfide scan.

The bug that was actually here: the scan required a Cβ atom in the structure, and
glycine has none. Every single glycine dropped out — including the pair
the literature uses most (Kabat VH44, which is Gly in trastuzumab). The tests
below pin down both the geometry and that pair.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from abx import thermo  # noqa: E402

PDB = ROOT / "data" / "1n8z.pdb"
needs_structure = pytest.mark.skipif(not PDB.exists(), reason="1N8Z is not here")


# ── the computed Cβ ─────────────────────────────────────────────────────────
@needs_structure
def test_computed_cb_matches_the_real_one():
    """Compute Cβ for residues that HAVE one, and compare. Otherwise the formula is faithful."""
    deviations = []
    for r in thermo.load_residues(PDB):
        if r.cb is None or r.cb_is_virtual or r.resn == "GLY":
            continue
        atoms = _backbone(r)
        if atoms is None:
            continue
        deviations.append(float(np.linalg.norm(thermo.virtual_cb(*atoms) - r.cb)))
    assert len(deviations) > 200, "found too few residues to check against"
    assert max(deviations) < 0.5, f"largest deviation {max(deviations):.2f} Å"
    assert float(np.median(deviations)) < 0.2


@needs_structure
def test_glycine_gets_a_cb_and_it_is_marked():
    gly = [r for r in thermo.load_residues(PDB) if r.resn == "GLY"]
    assert gly, "no glycine in the structure?"
    assert all(r.cb is not None and r.cb_is_virtual for r in gly)


@needs_structure
def test_residues_with_real_cb_are_not_marked_as_computed():
    res = [r for r in thermo.load_residues(PDB) if r.resn not in ("GLY",)]
    assert any(r.cb is not None and not r.cb_is_virtual for r in res)
    assert not any(r.cb_is_virtual for r in res if r.resn == "PHE")


# ── the pair the literature uses ────────────────────────────────────────────
@needs_structure
def test_the_known_pair_vh44_vl100_is_found():
    """Kabat VH44 + VL100 is the most widely used engineered bridge between VH and VL.

    In trastuzumab both are sequential positions: heavy 44 (Gly), light 100
    (Gln). If the scan does not find that pair, it is looking in the wrong place.
    """
    residues = [r for r in thermo.load_residues(PDB) if thermo.in_variable_domain(r)]
    hits = thermo.disulfide_candidates(residues)
    pairs = {(h["a"].resi, h["b"].resi) for h in hits}  # (light, heavy)
    assert (100, 44) in pairs, f"only found {sorted(pairs)}"


@needs_structure
def test_known_bridges_are_still_found():
    """The control that shows the method works: the real cysteine bridges in 1N8Z."""
    residues = thermo.load_residues(PDB)
    pairs = {(a.chain + str(a.resi), b.chain + str(b.resi))
             for a, b, _ in thermo.existing_disulfides(residues)}
    assert ("A23", "A88") in pairs
    assert ("B22", "B96") in pairs


def _backbone(r):
    """N, CA and C for one residue, re-read from the file."""
    want = {"N": None, "CA": None, "C": None}
    for line in PDB.read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != r.chain:
            continue
        if int(line[22:26]) != r.resi or line[26].strip() != r.icode:
            continue
        name = line[12:16].strip()
        if name in want and want[name] is None:
            want[name] = np.array([float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])])
    if any(v is None for v in want.values()):
        return None
    return want["N"], want["CA"], want["C"]

