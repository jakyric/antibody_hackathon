"""Shared constants, data loading and small helpers for the HER2 CDR-H3 + Vernier pipeline."""
from __future__ import annotations

import csv
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml

# ---------------------------------------------------------------- sequences
WT_H3 = "WGGDGFYAMD"          # trastuzumab CDR-H3 (positions 1-10)
AAS = "ACDEFGHIKLMNPQRSTVWY"
AI = {a: i for i, a in enumerate(AAS)}

# Trastuzumab Fv. VH numbering below is sequential: CDR-H3 = residues 99-108, H94 = R98.
WT_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
         "SRWGGDGFYAMDYWGQGTLVTVSS")
WT_VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
# HER2 domain IV (the 1N8Z construct incl. the segment 581-590 that is unresolved in the crystal)
HER2_DOM4 = ("CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQCVACAHYKDPPFCVARCPS"
             "GVKPDLSYMPIWKFPDEEGACQPCPIN")
H3_START = 99                 # sequential VH index of CDR-H3 position 1 (1-based)

# Vernier-zone residues that can touch CDR-H3 (1N8Z analysis): name -> (chain, 1-based index, WT residue, H3 positions touched)
VERNIER = {
    "H94": ("H", 98, "R", {1, 2, 10}),
    "H93": ("H", 97, "S", {9}),
    "H47": ("H", 47, "W", {9}),
    "L36": ("L", 36, "Y", {8, 9}),
    "L46": ("L", 46, "L", {6, 8, 9, 10}),
    "L49": ("L", 49, "Y", {6, 8}),
}

# CDR-H3 amino-acid background used by the paper's BLOSUM recipe (Chinery et al., SAbDab CDR-H3 counts)
H3_COUNTS = {"A": 6710, "C": 555, "D": 6442, "E": 1636, "F": 3404, "G": 6582, "H": 1081, "I": 1342, "K": 970,
             "L": 2509, "M": 1218, "N": 1658, "P": 2090, "Q": 767, "R": 5196, "S": 4350, "T": 2835, "V": 2975,
             "W": 1686, "Y": 9028}

CONTROLS = {  # id -> CDR-H3; measured binders and non-binders used to calibrate every metric
    "ctrl_wt": "WGGDGFYAMD",
    "ctrl_base": "WGSQGFYAYD",
    "ctrl_nonbinder1": "KGGDNWYAMD",
    "ctrl_nonbinder2": "EYAYARGYFE",
}


def mutations(seq: str, ref: str) -> list[str]:
    """Substitutions of seq relative to ref as 'G2A' style strings (1-based CDR-H3 positions)."""
    return [f"{ref[i]}{i + 1}{seq[i]}" for i in range(len(ref)) if seq[i] != ref[i]]


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


# ---------------------------------------------------------------- config / paths
def load_config(path: str | Path) -> dict:
    cfg = yaml.safe_load(open(path))
    cfg["_path"] = str(path)
    out = Path(cfg["paths"]["out_dir"]).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    cfg["_out"] = out
    return cfg


def out(cfg: dict, name: str) -> Path:
    p = cfg["_out"] / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def save_summary(cfg: dict, **kv) -> None:
    """Counts written by each stage; the deck stage reads them back."""
    f = out(cfg, "summary.json")
    d = json.loads(f.read_text()) if f.exists() else {}
    d.update(kv)
    f.write_text(json.dumps(d, indent=1))


def read_summary(cfg: dict) -> dict:
    f = out(cfg, "summary.json")
    return json.loads(f.read_text()) if f.exists() else {}


# ---------------------------------------------------------------- hackathon data
def spr_table(cfg: dict) -> list[dict]:
    return [d for d in csv.DictReader(open(Path(cfg["paths"]["data_dir"]).expanduser() / "trastuzumab_her2_spr_all_768.csv"))]


def spr_kd(cfg: dict) -> dict[str, float]:
    """Mean K_D (nM) per unique quantified 10-mer CDR-H3."""
    agg = defaultdict(list)
    for d in spr_table(cfg):
        if len(d["cdrh3"]) == 10 and d["binding_result"] == "quantified_binding" and d["kd_nM"]:
            agg[d["cdrh3"]].append(float(d["kd_nM"]))
    return {s: float(np.mean(v)) for s, v in agg.items()}


def known_sequences(cfg: dict) -> set[str]:
    """Every CDR-H3 in the provided data (FACS screen + SPR designs): 525,063 sequences."""
    data = Path(cfg["paths"]["data_dir"]).expanduser()
    known = {d["cdrh3"] for d in spr_table(cfg)}
    with zipfile.ZipFile(data / "tz_her2_facs_3classes_524346.csv.zip") as z:
        name = [x for x in z.namelist() if x.endswith(".csv") and "MACOS" not in x][0]
        with z.open(name) as fh:
            known.update(d["cdrh3"] for d in csv.DictReader(io.TextIOWrapper(fh, "utf-8")))
    return known


def onehot(seqs: list[str]) -> np.ndarray:
    X = np.zeros((len(seqs), 10 * len(AAS)))
    for r, s in enumerate(seqs):
        for i, a in enumerate(s):
            X[r, i * len(AAS) + AI[a]] = 1
    return X


def fit_ridge(cfg: dict) -> np.ndarray:
    """Additive ridge model on log10 K_D of the SPR designs: a (10, 20) coefficient matrix (negative = tighter)."""
    agg = defaultdict(list)
    for d in spr_table(cfg):
        if len(d["cdrh3"]) == 10 and d["binding_result"] == "quantified_binding" and d["kd_nM"]:
            agg[d["cdrh3"]].append(np.log10(float(d["kd_nM"])))
    seqs = sorted(agg)
    y = np.array([np.mean(agg[s]) for s in seqs])     # replicates averaged on the log scale
    X = onehot(seqs)
    Xc, yc = X - X.mean(0), y - y.mean()                     # unpenalised intercept = centre both sides
    w = np.linalg.solve(Xc.T @ Xc + cfg["ridge"]["alpha"] * np.eye(X.shape[1]), Xc.T @ yc)
    return w.reshape(10, len(AAS))


def ridge_score(seq: str, coef: np.ndarray, ref: str = WT_H3) -> float:
    """Predicted change in log10 K_D of seq relative to ref (lower = tighter). Ranking aid only."""
    return float(sum(coef[i, AI[a]] - coef[i, AI[ref[i]]] for i, a in enumerate(seq)))


def blosum(name: str = "BLOSUM45") -> dict[str, dict[str, int]]:
    from Bio.Align import substitution_matrices
    B = substitution_matrices.load(name)
    return {a: {b: int(B[a][b]) for b in AAS} for a in AAS}


def build_fv(h3: str, vernier: dict[str, str] | None = None) -> tuple[str, str]:
    """VH and VL with a given CDR-H3 and optional Vernier substitutions {name: new residue}."""
    vh = list(WT_VH.replace(WT_H3, h3))
    vl = list(WT_VL)
    for name, aa in (vernier or {}).items():
        chain, idx, wt, _ = VERNIER[name]
        seq = vh if chain == "H" else vl
        assert seq[idx - 1] == wt, f"{name}: expected {wt} at {idx}"
        seq[idx - 1] = aa
    return "".join(vh), "".join(vl)


def fv_yaml(vh: str, vl: str) -> str:
    """Boltz-style YAML used by JapanFold: chain A = VH, B = VL, C = HER2 domain IV."""
    return (f"sequences:\n  - protein: {{id: A, sequence: {vh}}}\n  - protein: {{id: B, sequence: {vl}}}\n"
            f"  - protein: {{id: C, sequence: {HER2_DOM4}}}\n")
