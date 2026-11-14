"""What the measured data and the crystal structure already say about each CDR-H3 position.

Two sources of precedent, both independent of any model we have trained:

1. **Measured single mutants in the SPR table.** 15 of the 768 rows are one mutation
   from trastuzumab. They cover 7 of 10 positions and tell us what it *costs*
   to touch the position — measured, not predicted. All three measured M107 mutants
   lose 40–70-fold affinity; a design that swaps M107 therefore does so against
   known counter-evidence.
2. **The geometry in 1N8Z.** Four positions touch HER2. Three others hold the loop
   in place without touching the antigen, and are therefore "safe" only if one looks at
   contact distance and ignores why the loop has the shape it has.

Both are borrowed from the teammates' branches (jakob: structure and
single-mutant precedent; moritz: the contact positions measured from 1N8Z) and have been
recomputed against our own data files here.

BLOSUM45 similarity to trastuzumab is included for the same reason as
`edit_distance`: Ucar & Sormanni (bioRxiv 2025.10.26.684652) show that similarity
to the starting point alone ranks variants as well as generative models
on HER2 data. A model that does not beat it has not shown that it has learned anything.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from statistics import median

from abx.data import TZ_CDRH3_10

# ── What each position in the 10-mer is called in 1N8Z, chain B ─────────────
# WGGDGFYAMD = B99..B108. The order is verified in the tests against
# TRASTUZUMAB_CDRH3, so a change there shows up here.
POSITION_LABELS: tuple[str, ...] = (
    "W99", "G100", "G101", "D102", "G103", "F104", "Y105", "A106", "M107", "D108",
)

# Positions with a heavy atom closer than 4.5 Å to HER2 in 1N8Z.
HER2_CONTACT_POSITIONS = frozenset({0, 4, 5, 6})  # W99, G103, F104, Y105

# Positions that do not touch HER2 but that hold the loop. Swapping them is not
# free even though the contact distance says they are free.
# Note: none of the texts may contain a semicolon. It is the separator when
# several notes are joined into one CSV cell, and a text with a semicolon in
# it would be read as two notes.
# The texts are in English because they end up directly in shared result files.
STRUCTURAL_RISK: Mapping[int, str] = {
    2: "G101 has a left-handed backbone angle, so a side chain here reshapes the loop",
    7: "A106 is almost fully buried against VL, so a large or charged residue strains it",
    9: "D108 holds a salt bridge to the conserved ArgH94, and losing the carboxylate removes the anchor",
}

# ── BLOSUM45, NCBI standard, the 20 common residues ─────────────────────────
_BLOSUM45_ORDER = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM45_ROWS = """
 5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0
-2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2
-1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3
-2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3
-1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1
-1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3
-1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3
 0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3
-2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3
-1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3
-1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1
-1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2
-1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1
-2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0
-1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3
 1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0
-2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3
-2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1
 0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5
"""


def _build_blosum45() -> dict[tuple[str, str], int]:
    rows = [r.split() for r in _BLOSUM45_ROWS.strip().splitlines()]
    return {
        (a, b): int(rows[i][j])
        for i, a in enumerate(_BLOSUM45_ORDER)
        for j, b in enumerate(_BLOSUM45_ORDER)
    }


BLOSUM45 = _build_blosum45()


def blosum45_similarity(seq: str, reference: str = TZ_CDRH3_10) -> int:
    """Sum of BLOSUM45 scores between `seq` and `reference`, position by position.

    High score = similar to the starting point. Used as a baseline, not as a prediction.
    """
    if len(seq) != len(reference):
        msg = f"lengths differ: {len(seq)} vs {len(reference)}"
        raise ValueError(msg)
    return sum(BLOSUM45[(a, b)] for a, b in zip(reference, seq, strict=True))


# ── Measured single mutants ─────────────────────────────────────────────────
def wildtype_kd(rows: Iterable[Mapping[str, str]], reference: str = TZ_CDRH3_10) -> float:
    """Median measured KD for trastuzumab itself. The table has it several times."""
    vals = [
        float(r["kd_nM"])
        for r in rows
        if r["cdrh3"] == reference and r["kd_nM"] not in ("", "NA")
    ]
    if not vals:
        msg = "no measured KD found for the reference"
        raise ValueError(msg)
    return median(vals)


def single_mutant_table(
    rows: Iterable[Mapping[str, str]], reference: str = TZ_CDRH3_10
) -> dict[tuple[int, str], dict[str, float | int | str]]:
    """{(position, new residue): measured KD and loss versus wild type} for each measured single mutant.

    Only rows with a quantified KD count. A mutant without a reportable KD
    says something else (it did not bind measurably), and is not mixed in here.
    """
    rows = list(rows)
    wt = wildtype_kd(rows, reference)
    out: dict[tuple[int, str], dict[str, float | int | str]] = {}
    for r in rows:
        seq = r["cdrh3"]
        if len(seq) != len(reference) or r["binding_result"] != "quantified_binding":
            continue
        diff = [i for i, (a, b) in enumerate(zip(reference, seq, strict=True)) if a != b]
        if len(diff) != 1:
            continue
        i = diff[0]
        kd = float(r["kd_nM"])
        key = (i, seq[i])
        prev = out.get(key)
        if prev is None or kd < float(prev["kd_nM"]):
            out[key] = {"kd_nM": kd, "fold_vs_wt": kd / wt, "position": POSITION_LABELS[i]}
    return out


def position_precedent(
    singles: Mapping[tuple[int, str], Mapping[str, object]], index: int
) -> list[float]:
    """All measured losses for a position, regardless of which residue was introduced."""
    return sorted(
        float(v["fold_vs_wt"]) for (i, _), v in singles.items() if i == index
    )


def design_evidence(
    cdrh3: str,
    singles: Mapping[tuple[int, str], Mapping[str, object]],
    reference: str = TZ_CDRH3_10,
) -> dict[str, object]:
    """Precedent and risk flags for one design.

    Returns the mutations, which of them have a measured single mutant behind them,
    the worst measured loss among them, how many HER2 contacts the design touches,
    and the structural notes that apply.
    """
    if len(cdrh3) != len(reference):
        msg = f"lengths differ: {len(cdrh3)} vs {len(reference)}"
        raise ValueError(msg)
    muts = [i for i, (a, b) in enumerate(zip(reference, cdrh3, strict=True)) if a != b]

    exact: list[str] = []
    worst_exact = 0.0
    for i in muts:
        hit = singles.get((i, cdrh3[i]))
        if hit is not None:
            fold = float(hit["fold_vs_wt"])
            exact.append(f"{POSITION_LABELS[i]}{cdrh3[i]}:{fold:.0f}x")
            worst_exact = max(worst_exact, fold)

    # Positions where *any* measured single mutant lost at least tenfold.
    against_precedent: list[str] = []
    for i in muts:
        folds = position_precedent(singles, i)
        if folds and min(folds) >= 10.0:
            against_precedent.append(
                f"{POSITION_LABELS[i]}: all {len(folds)} measured mutants lost "
                f"{min(folds):.0f}-{max(folds):.0f}x"
            )

    notes = [STRUCTURAL_RISK[i] for i in muts if i in STRUCTURAL_RISK]
    contacts = sorted(POSITION_LABELS[i] for i in muts if i in HER2_CONTACT_POSITIONS)
    return {
        "cdrh3": cdrh3,
        "n_mutations": len(muts),
        "mutations": ";".join(f"{POSITION_LABELS[i]}{cdrh3[i]}" for i in muts),
        "blosum45": blosum45_similarity(cdrh3, reference),
        "her2_contacts_touched": ";".join(contacts),
        "n_her2_contacts_touched": len(contacts),
        "measured_singles": ";".join(exact),
        "worst_measured_fold": round(worst_exact, 1) if exact else "",
        "against_measured_precedent": ";".join(against_precedent),
        "structural_risk": ";".join(notes),
        "n_risk_flags": len(against_precedent) + len(notes),
    }
