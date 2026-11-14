"""Enumerate CDR-H3 loops that drop trastuzumab's own liability motifs.

WHY THIS EXISTS
  Scoring a published comparison set with our own tools showed one real gap,
  and it was not in our favour: nine of their ten designs reduced the
  high-severity liability count against trastuzumab, where four of our twenty
  did, and nine of ours introduced the acid-labile `FragH:DP` motif that none
  of theirs carried.

  That gap is not a modelling problem. It is a search problem. Our twenty were
  drawn at random from a position-weight matrix and then inspected for
  liabilities afterwards; nothing in the draw knew what a liability was. The
  loops that clear the bar exist, they were simply never enumerated.

WHAT A CLEAN LOOP MEANS HERE
  Trastuzumab's CDR-H3 is WGGDGFYAMD and carries three motifs of its own:

      Isom:DG   the D4-G5 pair, isomerisation, HIGH severity
      MetOx:M   the M9, methionine oxidation
      TrpOx:W   the W1, tryptophan oxidation

  Removing all three costs three substitutions (W1, M9, and one of D4/G5) and
  the novelty rule requires four, so exactly one position is free. The space
  is therefore small enough to enumerate completely rather than sample: every
  legal loop that drops all three motifs and adds none can be written down and
  scored, with no draw and no seed.

WHAT THIS MODULE DOES NOT DO
  It does not decide anything. The motif list here is our reading of what
  the official checker reported on our own designs, used to narrow a search;
  the counts that rank candidates come from the official checker itself, run
  on the survivors with `--all-liabilities`. A motif is a sequence pattern,
  not a measurement.
"""

from __future__ import annotations

from itertools import combinations, product

from .data import TZ_CDRH3_10

# Residues never introduced into the loop.
#   C  an unpaired cysteine is a hard gate failure
#   M  methionine oxidation, the motif we are removing
#   W  tryptophan oxidation, the motif we are removing
FORBIDDEN_RESIDUES = frozenset("CMW")

# Dipeptides avoided inside the loop and across its two joins with the
# framework (...YYCSR-loop-YWGQGT...).
#
# THIS LIST IS A SEARCH FILTER, NOT AN AUTHORITY. It started as DG, DP and the
# obvious deamidation pairs, and a probe of eight enumerated loops through
# the official checker (`--all-liabilities`) showed three of them still
# carrying an Isom motif in CDR-H3: the tool reads `DT` as isomerisation too,
# which this list did not. It now covers the classic isomerisation (D followed
# by a small or polar residue) and deamidation (N followed by the same)
# series, which is wider than what the official checker flagged rather than
# narrower.
#
# The finalists are still validated by the official checker itself. A motif
# list we wrote is a way to narrow a search, never a way to declare a design
# clean.
#   DP      acidic fragmentation [high]
#   D[GSTDHN]  isomerisation (DG and DT confirmed [high])
#   N[GSTADHN] deamidation (NG [high], NT [medium] confirmed)
#   TN      deamidation [low]
FORBIDDEN_PAIRS = (
    "DP",
    "DG", "DS", "DT", "DD", "DH", "DN",
    "NG", "NS", "NT", "NA", "ND", "NH", "NN",
    "TN",
)

LEFT_FLANK = "R"   # the residue immediately before the loop
RIGHT_FLANK = "Y"  # the residue immediately after it

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


# CDR-H3 residues whose heavy atoms come within 4.5 A of HER2 in 1N8Z, as
# positions in the 10-mer (0-indexed): W99, G103, F104, Y105. Measured from
# the crystal structure, not assumed, and identical at a 5.0 A cutoff -- so
# the set is not an artefact of where the cutoff was put. The other six
# positions do not touch the antigen at either cutoff.
#
# This matters because trastuzumab's TrpOx liability IS W99, a contact. You
# cannot remove it without changing a residue that touches HER2. Isomerisation
# (D102) and Met oxidation (M107) sit at 5.1 A and 9.4 A and are free.
CONTACT_POSITIONS = frozenset({0, 4, 5, 6})


def motif_positions(loop: str, pair: str) -> list[int]:
    """Where `pair` occurs in the loop, counting the framework joins."""
    padded = LEFT_FLANK + loop + RIGHT_FLANK
    return [i - 1 for i in range(len(padded) - 1) if padded[i:i + 2] == pair]


def is_clean(loop: str) -> bool:
    """No forbidden residue, and no forbidden dipeptide including the joins."""
    if any(a in FORBIDDEN_RESIDUES for a in loop):
        return False
    padded = LEFT_FLANK + loop + RIGHT_FLANK
    return not any(p in padded for p in FORBIDDEN_PAIRS)


def is_clean_outside(loop: str, frozen, parent: str = TZ_CDRH3_10) -> bool:
    """Clean everywhere the design is allowed to touch.

    A frozen position keeps the parent's residue by construction, so a
    liability sitting on one is not something the design introduced and not
    something it can remove. This asks the question that can be acted on: did
    we add a motif, or leave one we could have taken out?

    Dipeptides are still checked across the whole loop including the framework
    joins, because a pair formed between a frozen residue and a mutated one IS
    the design's doing.
    """
    for i, a in enumerate(loop):
        if a in FORBIDDEN_RESIDUES and i not in frozen:
            return False
        if i in frozen and a != parent[i]:
            return False
    padded = LEFT_FLANK + loop + RIGHT_FLANK
    return not any(p in padded for p in FORBIDDEN_PAIRS)


def enumerate_contact_safe(
    n_mut: int = 4, frozen=CONTACT_POSITIONS, parent: str = TZ_CDRH3_10
) -> list[str]:
    """Every loop at `n_mut` substitutions that touches no contact residue.

    The trade this encodes, stated plainly: keeping W99 means keeping the Trp
    oxidation liability (medium severity). It buys the whole measured
    paratope. The judged axes are binding, epitope and thermostability;
    liabilities are judged by nobody and matter only through expression. So
    a contact is worth more than a medium-severity motif, and this enumeration
    is the version of the search that believes that.
    """
    allowed = [a for a in ALPHABET if a not in FORBIDDEN_RESIDUES]
    free = [p for p in range(len(parent)) if p not in frozen]

    # Positions that must change for a removable motif to go. A forbidden pair
    # whose every position is frozen cannot be removed and is not required.
    forced = {p for p in free if parent[p] in FORBIDDEN_RESIDUES}
    choices: list[list[int]] = [[p] for p in sorted(forced)]
    for pair in FORBIDDEN_PAIRS:
        for start in motif_positions(parent, pair):
            inside = [p for p in (start, start + 1) if p in free]
            if not inside:
                continue  # wholly on frozen ground; cannot be removed
            if any(p in forced for p in inside):
                continue
            choices.append(inside)

    seen: set[str] = set()
    for combo in product(*choices) if choices else [()]:
        base = set(combo)
        if len(base) > n_mut:
            continue
        rest = [p for p in free if p not in base]
        for extra in combinations(rest, n_mut - len(base)):
            sites = sorted(base | set(extra))
            options = [[a for a in allowed if a != parent[p]] for p in sites]
            for residues in product(*options):
                loop = list(parent)
                for p, a in zip(sites, residues):
                    loop[p] = a
                s = "".join(loop)
                if is_clean_outside(s, frozen, parent):
                    seen.add(s)
    return sorted(seen)


def required_positions(parent: str = TZ_CDRH3_10) -> set[int]:
    """Positions that must change for the parent's own motifs to disappear.

    Every forbidden residue in the parent must go, and every forbidden pair
    must be broken -- but a pair can be broken at either of its two positions,
    so those are returned as alternatives by `mutation_sites`, not here.
    """
    return {i for i, a in enumerate(parent) if a in FORBIDDEN_RESIDUES}


def mutation_sites(n_mut: int, parent: str = TZ_CDRH3_10) -> list[tuple[int, ...]]:
    """Every set of exactly `n_mut` positions that can clean the parent.

    Each forbidden pair in the parent contributes one position from that pair,
    each forbidden residue contributes its own position, and the remainder are
    free. Positions outside the loop (a pair straddling a framework join)
    cannot be mutated and are skipped -- no such pair exists in trastuzumab,
    and if one appeared the loop could not be cleaned at all, which is a thing
    to find out rather than silently ignore.
    """
    forced = required_positions(parent)
    choices: list[list[int]] = [[p] for p in sorted(forced)]
    for pair in FORBIDDEN_PAIRS:
        for start in motif_positions(parent, pair):
            inside = [p for p in (start, start + 1) if 0 <= p < len(parent)]
            if not inside:
                raise ValueError(f"{pair} at the framework join cannot be removed")
            if any(p in forced for p in inside):
                continue
            choices.append(inside)

    out: set[tuple[int, ...]] = set()
    for combo in product(*choices) if choices else [()]:
        base = set(combo)
        if len(base) > n_mut:
            continue
        free = [p for p in range(len(parent)) if p not in base]
        for extra in combinations(free, n_mut - len(base)):
            out.add(tuple(sorted(base | set(extra))))
    return sorted(out)


def enumerate_clean(n_mut: int = 4, parent: str = TZ_CDRH3_10) -> list[str]:
    """Every loop at exactly `n_mut` substitutions that carries no motif.

    Complete, not sampled: the same list every time, with no seed to report.
    """
    allowed = [a for a in ALPHABET if a not in FORBIDDEN_RESIDUES]
    seen: set[str] = set()
    for sites in mutation_sites(n_mut, parent):
        options = [[a for a in allowed if a != parent[p]] for p in sites]
        for combo in product(*options):
            loop = list(parent)
            for p, a in zip(sites, combo):
                loop[p] = a
            s = "".join(loop)
            if is_clean(s):
                seen.add(s)
    return sorted(seen)

