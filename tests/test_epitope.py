"""Tests for the epitope measurement.

What can go wrong without anyone noticing: the numbering. HER2 is numbered 1.. in the
predicted files and 489.. in the crystal. If they are mixed up, the overlap is zero without
anything actually being wrong, and we would think that all the designs lost the epitope.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ep = load("epitope", "scripts/epitope.py")


def atoms(**chains):
    return {c: [(i + 1, np.array(xyz, dtype=float)) for i, xyz in enumerate(v)]
            for c, v in chains.items()}


def test_only_residues_within_the_cutoff_are_counted():
    a = atoms(B=[(0, 0, 0)], C=[(0, 0, 3.0), (0, 0, 10.0)])
    assert ep.contact_residues(a, ("B",), ("C",), cutoff=4.5) == {1}


def test_numbering_is_shifted_to_the_crystal_numbering():
    a = atoms(B=[(0, 0, 0)], C=[(0, 0, 3.0)], D=[(0, 0, 3.0)])
    hits = ep.contact_residues(a, ("B",), ("C", "D"), ep.AG_OFFSET, cutoff=4.5)
    # C residue 1 is HER2 489, D residue 1 is HER2 591
    assert hits == {489, 591}


def test_without_the_offset_the_answer_is_different():
    """The control that shows the offset actually does something."""
    a = atoms(B=[(0, 0, 0)], C=[(0, 0, 3.0)])
    assert ep.contact_residues(a, ("B",), ("C",), cutoff=4.5) == {1}


def test_empty_antibody_gives_no_epitope():
    assert ep.contact_residues(atoms(C=[(0, 0, 0)]), ("B",), ("C",)) == set()


def test_block_splitting_does_not_change_the_answer():
    """Distances are computed 512 atoms at a time. The boundary should not be visible."""
    rng = np.random.default_rng(0)
    ag = [tuple(x) for x in rng.normal(scale=6, size=(1500, 3))]
    a = atoms(B=[(0.0, 0.0, 0.0)], C=ag)
    got = ep.contact_residues(a, ("B",), ("C",), cutoff=4.5)
    want = {i + 1 for i, x in enumerate(ag) if np.linalg.norm(x) <= 4.5}
    assert got == want


def test_jaccard_is_as_expected():
    assert ep.jaccard({1, 2, 3}, {1, 2, 3}) == 1.0
    assert ep.jaccard({1, 2}, {3, 4}) == 0.0
    assert ep.jaccard({1, 2, 3, 4}, {3, 4}) == pytest.approx(0.5)
    assert np.isnan(ep.jaccard(set(), set()))
