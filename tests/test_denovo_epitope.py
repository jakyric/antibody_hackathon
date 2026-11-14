"""Tests for the epitope measurement.

The geometry is easy to make look correct while still being wrong: wrong chain,
wrong column, wrong direction of the distance. That is why it is tested on
coordinates where the answer is known in advance.
"""

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ep = load("denovo_epitope", "scripts/denovo_epitope.py")


def atom(chain, resnum, x, y=0.0, z=0.0, name="CA"):
    return (chain, resnum, name, x, y, z)


def test_only_residues_within_the_cutoff_are_included():
    atoms = [
        atom("B", 1, 0.0),          # antibody
        atom("C", 10, 4.0),         # 4 Å away  -> included
        atom("C", 11, 5.0),         # 5 Å away  -> not included
        atom("D", 12, 4.4),         # just barely included, and in the other antigen chain
    ]
    assert ep.interface(atoms) == {("C", 10), ("D", 12)}


def test_the_cutoff_is_inclusive():
    atoms = [atom("A", 1, 0.0), atom("C", 7, ep.CUTOFF)]
    assert ep.interface(atoms) == {("C", 7)}


def test_contact_between_two_antigen_residues_does_not_count():
    atoms = [atom("C", 1, 0.0), atom("D", 2, 1.0)]
    assert ep.interface(atoms) == set()


def test_jaccard():
    assert ep.jaccard({1, 2, 3}, {2, 3, 4}) == 0.5
    assert ep.jaccard({1}, {1}) == 1.0
    assert ep.jaccard({1}, {2}) == 0.0


def test_residue_numbers_are_translated_to_her2_own_numbering(tmp_path):
    q = tmp_path / "queries.json"
    q.write_text(json.dumps({"queries": [{"query_id": "x", "chains": [
        {"chain_id": "A", "sequence": "DIQ"},
        {"chain_id": "C", "sequence": "AAA", "residue_ids": [570, 571, 572]},
        {"chain_id": "D", "sequence": "AA", "residue_ids": [591, 592]},
    ]}]}))
    mapping = ep.residue_id_map(q)
    assert ep.renumber({("C", 3), ("D", 1)}, mapping) == {572, 591}


def test_unknown_position_becomes_negative_and_does_not_disappear_silently():
    assert ep.renumber({("C", 9)}, {"C": {1: 570}}) == {-9}


def test_original_numbering_is_kept_as_it_is():
    """Some models write HER2 with its own numbers, others with 489, 490 ...
    If the number is already in residue_ids, it must not be translated again."""
    mapping = {"C": {1: 570, 2: 571, 3: 572}}
    assert ep.renumber({("C", 572)}, mapping) == {572}
    assert ep.renumber({("C", 3)}, mapping) == {572}


def test_cif_is_read_with_the_columns_from_the_header(tmp_path):
    cif = tmp_path / "m.cif"
    cif.write_text(
        "loop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n"
        "ATOM 1 N N ASP A 0.0 0.0 0.0 5 B\n"
        "ATOM 2 H H1 ASP A 1.0 0.0 0.0 5 B\n"
        "ATOM 3 C CA ALA C 3.0 0.0 0.0 77 C\n"
        "#\n"
    )
    atoms = ep.read_atoms(cif)
    assert [(a[0], a[1]) for a in atoms] == [("B", 5), ("C", 77)]  # the hydrogen is left out
    assert ep.interface(atoms) == {("C", 77)}
