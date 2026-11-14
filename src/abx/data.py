"""Loading of the two trastuzumab tables, with the conventions made explicit.

Source: Chinery et al., prepared by the organisers (UiO).
    data/tz_her2_facs_3classes_524346.csv   524,346 variants, FACS class
    data/trastuzumab_her2_spr_all_768.csv   768 variants, 190 with measured KD

THE CDR-H3 CONVENTION, which is worth reading twice:
    The dataset calls CDR-H3 a **10-mer**. The heavy chain is built as
        ...YYCSR  +  <10-mer>  +  Y  +  WGQGT...
    The final Y is thus NOT included in the `cdrh3` column, even though most
    numbering schemes (Kabat, IMGT) count it. Trastuzumab's own loop is
    called `WGGDGFYAMD` here, not `WGGDGFYAMDY`.

    The novelty filter in the task is measured against these tables. It must
    therefore be computed in THEIR convention. If the two are mixed, every
    distance is off by one residue, and a design that looks new gets discarded.
"""

from __future__ import annotations

import csv
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[2] / "data"

FACS_CSV = DATA / "tz_her2_facs_3classes_524346.csv"
FACS_ZIP = DATA / "tz_her2_facs_3classes_524346.csv.zip"
SPR_CSV = DATA / "trastuzumab_her2_spr_all_768.csv"

# Trastuzumab's CDR-H3 in the dataset's convention (10 residues).
TZ_CDRH3_10 = "WGGDGFYAMD"
# The anchors that surround it in the VH.
VH_PREFIX = "YYCSR"
VH_SUFFIX = "YWGQGT"          # note the Y: it belongs to the loop in Kabat, not here

# The other antibody in the SPR table. Different VH and VL; to be dropped when
# modelling trastuzumab variants (stated outright in the data README).
OTHER_AB_CDRH3 = "ARDRIVDDYYYYGMDV"


@dataclass(frozen=True)
class Facs:
    cdrh3: np.ndarray          # (N,) str, all of length 10
    binding_class: np.ndarray  # (N,) 'high' | 'medium' | 'low'
    edit_distance: np.ndarray  # (N,) int, number of positions that differ from trastuzumab
    vh: np.ndarray             # (N,) str

    def __len__(self) -> int:
        return len(self.cdrh3)

    @property
    def is_high(self) -> np.ndarray:
        return self.binding_class == "high"

    @property
    def is_low(self) -> np.ndarray:
        return self.binding_class == "low"


@dataclass(frozen=True)
class Spr:
    cdrh3: np.ndarray
    kd_nM: np.ndarray          # (N,) float, nan where KD is not reported
    binding_result: np.ndarray
    vh: np.ndarray
    vl: np.ndarray

    def __len__(self) -> int:
        return len(self.cdrh3)

    @property
    def has_kd(self) -> np.ndarray:
        return np.isfinite(self.kd_nM)


def _facs_rows():
    if FACS_CSV.exists():
        with FACS_CSV.open(newline="") as fh:
            yield from csv.DictReader(fh)
        return
    if FACS_ZIP.exists():
        with zipfile.ZipFile(FACS_ZIP) as zf:
            name = next(n for n in zf.namelist()
                        if n.endswith(".csv") and not n.startswith("__MACOSX"))
            with zf.open(name) as raw:
                text = (line.decode() for line in raw)
                yield from csv.DictReader(text)
        return
    raise FileNotFoundError(
        f"Found neither {FACS_CSV.name} nor {FACS_ZIP.name} in {DATA}. Run: make data"
    )


def load_facs() -> Facs:
    cdrh3, cls, ed, vh = [], [], [], []
    for row in _facs_rows():
        cdrh3.append(row["cdrh3"])
        cls.append(row["binding_class"])
        ed.append(int(row["edit_distance"]))
        vh.append(row["vh_sequence"])
    return Facs(np.array(cdrh3), np.array(cls), np.array(ed, dtype=np.int16), np.array(vh))


def load_spr(drop_other_antibody: bool = True) -> Spr:
    cdrh3, kd, res, vh, vl = [], [], [], [], []
    with SPR_CSV.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if drop_other_antibody and row["cdrh3"] == OTHER_AB_CDRH3:
                continue
            cdrh3.append(row["cdrh3"])
            raw = row["kd_nM"].strip()
            kd.append(float(raw) if raw else np.nan)
            res.append(row["binding_result"])
            vh.append(row["vh_sequence"])
            vl.append(row["vl_sequence"])
    return Spr(np.array(cdrh3), np.array(kd, dtype=float),
               np.array(res), np.array(vh), np.array(vl))


def provided_cdrh3s() -> set[str]:
    """All CDR-H3s the task counts as "provided data".

    The data README says 525,063 unique. That is the union of the two tables.
    A design must lie more than 2 point mutations away from every single one of them.
    """
    s = {str(x) for x in load_facs().cdrh3}
    s |= {str(x) for x in load_spr(drop_other_antibody=False).cdrh3}
    return s


def cdrh3_from_vh(vh: str) -> str | None:
    """Extract the CDR-H3 (the dataset's 10-mer convention) from a full VH domain."""
    i = vh.find(VH_PREFIX)
    if i < 0:
        return None
    start = i + len(VH_PREFIX)
    j = vh.find(VH_SUFFIX, start)
    if j < 0:
        return None
    return vh[start:j]


def vh_from_cdrh3(cdrh3: str) -> str:
    """Build a full VH domain around a CDR-H3, on trastuzumab's framework."""
    from .sequences import TRASTUZUMAB_VH
    i = TRASTUZUMAB_VH.find(VH_PREFIX) + len(VH_PREFIX)
    j = TRASTUZUMAB_VH.find(VH_SUFFIX, i)
    return TRASTUZUMAB_VH[:i] + cdrh3 + TRASTUZUMAB_VH[j:]
