"""Thermostability: the axis the novelty filter does not touch.

WHY THIS IS THE CHEAPEST AXIS
  Both novelty requirements concern CDR-H3:
     requirement A  CDR-H3 identity below 70% against trastuzumab
     requirement B  more than 2 point mutations from each provided data point
  The framework is free. A stabilizing framework mutation therefore costs
  nothing in novelty, and counts on its own axis in the scoring.
  A design that does not beat affinity can still earn a point here.

WHAT THIS FILE DOES
  Finds places in the real structure where a disulfide bridge between the two
  variable domains is geometrically possible. No positions are quoted
  from memory: everything is measured in 1N8Z.

  A disulfide bridge has an Sγ–Sγ distance of about 2.05 Å. With no sulfur atoms to measure on,
  we use the backbone geometry to decide whether two side chains can reach
  each other at all:
      Cβ–Cβ   3.0 – 4.8 Å
      Cα–Cα   4.0 – 7.0 Å
  Glycine has no Cβ and therefore drops out on its own.

WARNING
  This is geometry, not a stability prediction. A pair that passes
  is a CANDIDATE. It has to be folded and preferably measured before anyone believes in it.
  It also changes the cysteine count: two new cysteines keep the parity
  even, one would break the hard filter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[2] / "data"

# Disulfide geometry, from the backbone
#
# The limits are set from a pair we KNOW works, not from a gut feeling.
# The most widely used engineered bridge between VH and VL (Kabat VH44 + VL100; in
# trastuzumab the sequential positions heavy 44 and light 100) lies in 1N8Z at
# Cβ–Cβ 5.18 Å and Cα–Cα 5.46 Å. A window that stops at 4.8 Å thus throws out
# exactly that pair. The side chains rotate when both become cysteine, so the Cβ distance
# in the wild type is a coarser measure than it appears. The Cα window is the same as
# in the literature's own geometric scans (roughly 4.2–6.6 Å).
CB_MIN, CB_MAX = 3.0, 5.5
CA_MIN, CA_MAX = 4.0, 7.0

# A pair that has been built and measured in the lab by others, on 4D5 — the antibody
# trastuzumab is humanized from. We use it as the ground truth that the scan looks
# in the right place. (light, heavy) in the 1N8Z numbering.
KNOWN_ENGINEERED_PAIR = (100, 44)

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Where Cβ would have been, computed from the backbone.

    WHY THIS HAS TO BE HERE
      Glycine has no Cβ. A scan that requires a Cβ atom therefore skips
      every single glycine — and the most widely used engineered bridge between VH and
      VL in the literature sits exactly on a glycine (Kabat VH44, which is Gly in
      trastuzumab). We did not see that pair until this function existed.
      When the glycine is mutated to cysteine, it GETS a Cβ, roughly here.

    The standard construction from the three backbone atoms. It hits a real Cβ
    to within a few tenths of an ångström; this is tested in test_thermo.py.
    """
    b = ca - n
    cc = c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


@dataclass(frozen=True)
class Residue:
    chain: str
    resi: int
    icode: str
    resn: str
    ca: np.ndarray
    cb: np.ndarray | None
    cb_is_virtual: bool = False

    @property
    def one(self) -> str:
        return THREE2ONE.get(self.resn, "X")

    def __str__(self) -> str:
        return f"{self.chain}{self.resi}{self.icode}{self.one}"


def load_residues(pdb: Path | None = None, chains=("A", "B")) -> list[Residue]:
    pdb = pdb or (DATA / "1n8z.pdb")
    atoms: dict[tuple, dict] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        ch = line[21]
        if ch not in chains:
            continue
        name = line[12:16].strip()
        if name not in ("N", "CA", "C", "CB"):
            continue
        key = (ch, int(line[22:26]), line[26].strip())
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        atoms.setdefault(key, {"resn": line[17:20].strip()})[name] = xyz

    out = []
    for (ch, resi, icode), d in atoms.items():
        if "CA" not in d:
            continue
        cb, virtual = d.get("CB"), False
        if cb is None and all(k in d for k in ("N", "C")):
            # glycine (and any residue without Cβ in the file): compute where it would have been
            cb, virtual = virtual_cb(d["N"], d["CA"], d["C"]), True
        out.append(Residue(ch, resi, icode, d["resn"], d["CA"], cb, virtual))
    return sorted(out, key=lambda r: (r.chain, r.resi, r.icode))


def disulfide_candidates(
    residues: list[Residue],
    chain_a: str = "A",
    chain_b: str = "B",
    exclude: set[tuple[str, int]] | None = None,
) -> list[dict]:
    """Pairs across the domain boundary where a disulfide bridge is geometrically possible."""
    exclude = exclude or set()
    A = [r for r in residues if r.chain == chain_a and r.cb is not None
         and (r.chain, r.resi) not in exclude and r.resn != "CYS"]
    B = [r for r in residues if r.chain == chain_b and r.cb is not None
         and (r.chain, r.resi) not in exclude and r.resn != "CYS"]

    hits = []
    for ra in A:
        for rb in B:
            dcb = float(np.linalg.norm(ra.cb - rb.cb))
            if not (CB_MIN <= dcb <= CB_MAX):
                continue
            dca = float(np.linalg.norm(ra.ca - rb.ca))
            if not (CA_MIN <= dca <= CA_MAX):
                continue
            hits.append({
                "a": ra, "b": rb, "cb_dist": dcb, "ca_dist": dca,
                # the closer to 3.8 Å, the more typical of a real bridge
                "deviation": abs(dcb - 3.8),
            })
    hits.sort(key=lambda h: h["deviation"])
    return hits


def existing_disulfides(residues: list[Residue], cutoff: float = 5.0) -> list[tuple]:
    """Cysteine pairs that already lie close enough to be bridged."""
    cys = [r for r in residues if r.resn == "CYS" and r.cb is not None]
    out = []
    for i, a in enumerate(cys):
        for b in cys[i + 1:]:
            d = float(np.linalg.norm(a.cb - b.cb))
            if d <= cutoff:
                out.append((a, b, d))
    return sorted(out, key=lambda t: t[2])


# The boundaries of the variable domains in 1N8Z, read from the sequence:
# VL ends at ...KVEIK (residue 107), VH at ...VTVSS (residue 113).
# Everything beyond this is CL and CH1 — they are NOT in an scFv, so a bridge
# there is useless to us no matter how nice the geometry is.
VL_LAST, VH_LAST = 107, 113

# Residues not to be touched without good reason. W in WGQGT is the J segment's
# conserved tryptophan and packs against the core of the immunoglobulin fold.
HEAVY_KEEP = {110}


def in_variable_domain(r) -> bool:
    if r.chain == "A":
        return r.resi <= VL_LAST
    if r.chain == "B":
        return r.resi <= VH_LAST
    return False


def cdr_positions_heavy() -> set[int]:
    """The CDR residues in chain B of 1N8Z. These are to be left alone: they bind."""
    # CDR-H3 is verified as 99–108 in scripts/stabddg_prep.py.
    # H1 and H2 were found with the same anchor method in tests.
    return set(range(26, 36)) | set(range(50, 66)) | set(range(99, 110))


def cdr_positions_light() -> set[int]:
    return set(range(24, 35)) | set(range(50, 57)) | set(range(89, 98))
