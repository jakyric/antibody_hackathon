"""Precedent from measured single mutants and from the geometry of 1N8Z."""

from __future__ import annotations

import pytest

from abx.data import TZ_CDRH3_10
from abx.precedent import (
    BLOSUM45,
    HER2_CONTACT_POSITIONS,
    POSITION_LABELS,
    blosum45_similarity,
    design_evidence,
    position_precedent,
    single_mutant_table,
    wildtype_kd,
)
from abx.sequences import TRASTUZUMAB_CDRH3


def row(cdrh3: str, kd: str, result: str = "quantified_binding") -> dict[str, str]:
    return {"cdrh3": cdrh3, "kd_nM": kd, "binding_result": result}


# ── the position labels ─────────────────────────────────────────────────────
def test_position_labels_match_the_sequence() -> None:
    """The label carries the residue it describes. Otherwise the risk rules point at the wrong place."""
    assert len(POSITION_LABELS) == len(TZ_CDRH3_10)
    for label, res in zip(POSITION_LABELS, TZ_CDRH3_10, strict=True):
        assert label[0] == res, f"{label} vs {res}"
    # The 10-mer is the 11-mer CDR without the final Y.
    assert TRASTUZUMAB_CDRH3.startswith(TZ_CDRH3_10)
    assert POSITION_LABELS[8] == "M107"
    assert POSITION_LABELS[9] == "D108"


def test_contact_positions_are_inside_the_loop() -> None:
    assert HER2_CONTACT_POSITIONS <= set(range(len(TZ_CDRH3_10)))


def test_no_risk_note_contains_the_separator() -> None:
    """Semicolons separate notes in the CSV cell. A note with one in its text
    would be counted as two."""
    from abx.precedent import STRUCTURAL_RISK

    for i, note in STRUCTURAL_RISK.items():
        assert ";" not in note, f"position {i}: {note}"


def test_structural_notes_are_counted_one_by_one() -> None:
    """Two notes must count as two, no matter how long the texts are."""
    ev = design_evidence("WGPDGFYIMG", {})  # G101P, A106I and D108G
    assert ev["structural_risk"].count(";") + 1 == 3
    assert ev["n_risk_flags"] == 3


# ── BLOSUM45 ────────────────────────────────────────────────────────────────
def test_blosum45_is_symmetric_and_has_known_values() -> None:
    for (a, b), v in BLOSUM45.items():
        assert v == BLOSUM45[(b, a)], f"{a}{b}"
    assert BLOSUM45[("W", "W")] == 15
    assert BLOSUM45[("C", "C")] == 12
    assert BLOSUM45[("W", "C")] == -5
    assert BLOSUM45[("I", "V")] == 3  # a conservative swap is positive


def test_blosum45_similarity_is_highest_for_the_reference() -> None:
    same = blosum45_similarity(TZ_CDRH3_10)
    assert same == sum(BLOSUM45[(a, a)] for a in TZ_CDRH3_10)
    # One conservative swap costs less than one wild one.
    conservative = blosum45_similarity("WGGDGFYAID")  # A106I
    wild = blosum45_similarity("WGGDGFYAWD")  # A106W
    assert same > conservative > wild


def test_blosum45_similarity_rejects_a_length_mismatch() -> None:
    with pytest.raises(ValueError, match="lengths differ"):
        blosum45_similarity("WGGDG")


# ── measured single mutants ─────────────────────────────────────────────────
def test_wildtype_kd_is_the_median_of_the_repeats() -> None:
    rows = [row(TZ_CDRH3_10, "0.1"), row(TZ_CDRH3_10, "0.2"), row(TZ_CDRH3_10, "0.3")]
    assert wildtype_kd(rows) == pytest.approx(0.2)


def test_single_mutant_table_keeps_only_one_mutation_away() -> None:
    rows = [
        row(TZ_CDRH3_10, "0.1"),
        row("WGGDGFYATD", "10.0"),  # M107T, one mutation
        row("WGGDGFYATE", "1.0"),   # two mutations, must not be included
    ]
    tab = single_mutant_table(rows)
    assert set(tab) == {(8, "T")}
    assert tab[(8, "T")]["fold_vs_wt"] == pytest.approx(100.0)
    assert tab[(8, "T")]["position"] == "M107"


def test_single_mutant_table_ignores_rows_without_a_quantified_kd() -> None:
    """A mutant without a reportable KD did not bind measurably; it is not a loss figure."""
    rows = [row(TZ_CDRH3_10, "0.1"), row("WGGDGFYATD", "", "no_reportable_kd")]
    assert single_mutant_table(rows) == {}


def test_single_mutant_table_keeps_the_best_measurement_of_a_repeat() -> None:
    rows = [row(TZ_CDRH3_10, "0.1"), row("WGGDGFYATD", "10.0"), row("WGGDGFYATD", "4.0")]
    assert single_mutant_table(rows)[(8, "T")]["kd_nM"] == pytest.approx(4.0)


def test_position_precedent_collects_every_measured_loss() -> None:
    rows = [row(TZ_CDRH3_10, "0.1"), row("WGGDGFYATD", "10.0"), row("WGGDGFYAKD", "5.0")]
    folds = position_precedent(single_mutant_table(rows), 8)
    assert folds == pytest.approx([50.0, 100.0])
    assert position_precedent(single_mutant_table(rows), 0) == []


# ── design ──────────────────────────────────────────────────────────────────
def test_design_evidence_on_the_reference_is_empty() -> None:
    ev = design_evidence(TZ_CDRH3_10, {})
    assert ev["n_mutations"] == 0
    assert ev["n_risk_flags"] == 0
    assert ev["her2_contacts_touched"] == ""


def test_design_evidence_flags_a_position_where_every_mutant_lost() -> None:
    rows = [
        row(TZ_CDRH3_10, "0.1"),
        row("WGGDGFYATD", "10.0"),  # 100x
        row("WGGDGFYAKD", "5.0"),   # 50x
    ]
    singles = single_mutant_table(rows)
    ev = design_evidence("WGGDGFYARD", singles)  # M107R, not itself measured
    assert ev["measured_singles"] == ""  # no measured R at that position
    assert "M107" in ev["against_measured_precedent"]
    assert ev["n_risk_flags"] >= 1


def test_design_evidence_reports_an_exact_measured_hit() -> None:
    rows = [row(TZ_CDRH3_10, "0.1"), row("WGGDGFYATD", "10.0")]
    ev = design_evidence("WGGDGFYATD", single_mutant_table(rows))
    assert ev["measured_singles"] == "M107T:100x"
    assert ev["worst_measured_fold"] == pytest.approx(100.0)


def test_design_evidence_counts_her2_contacts_and_structural_notes() -> None:
    # G101P touches no HER2 contact, but is a structural risk.
    ev = design_evidence("WGPDGFYAMD", {})
    assert ev["n_her2_contacts_touched"] == 0
    assert "G101" in ev["structural_risk"]
    # G103 is a contact.
    ev2 = design_evidence("WGGDPFYAMD", {})
    assert ev2["her2_contacts_touched"] == "G103"
    assert ev2["n_her2_contacts_touched"] == 1


def test_design_evidence_rejects_a_length_mismatch() -> None:
    with pytest.raises(ValueError, match="lengths differ"):
        design_evidence("WGG", {})
