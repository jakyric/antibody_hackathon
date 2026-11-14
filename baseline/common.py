"""Shared helpers for the trastuzumab CDRH3 binding baseline.

Everything here is local and dependency-light: numpy, pandas, h5py, scikit-learn.
The authors' Keras CNN is re-implemented in numpy and loaded from its .h5 weights, so no
TensorFlow is needed and no pickle is ever loaded.
"""
from pathlib import Path

import numpy as np
import pandas as pd

import os

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
# VARIANT=contacts reruns the whole pipeline into separate files/folders (…_contacts) without touching the first run
VARIANT = os.environ.get("VARIANT", "")
SUF = f"_{VARIANT}" if VARIANT else ""


def vname(name):
    """'blosum_candidates.csv' -> 'blosum_candidates_contacts.csv' when VARIANT=contacts."""
    stem, ext = os.path.splitext(name)
    return f"{stem}{SUF}{ext}"


REFOLD = ROOT / "results" / f"refold{SUF}"
NEB_IN = ROOT / "nebius-runs" / f"refold{SUF}"
JF_IN = Path(__file__).resolve().parent / f"refold_inputs{SUF}" / "jf"
FACS_CSV = DATA / "her2" / "facs_extract" / "tz_her2_facs_3classes_524346.csv"
SPR_CSV = DATA / "her2" / "trastuzumab_her2_spr_all_768.csv"
MEHTA = DATA / "mehta"
OUT = Path(__file__).resolve().parent / "outputs"

WT = "WGGDGFYAMD"  # trastuzumab CDRH3
L = 10
AA = "ACDEFGHIKLMNPQRSTVWY"  # alphabetical; verified against the authors' CNN on their test set
_PW = 21 ** np.arange(L - 1, -1, -1, dtype=np.int64)


def encode_idx(seqs, alphabet=AA):
    """(n, 10) int8 array of amino-acid indices. All sequences must have length 10."""
    seqs = list(seqs)
    if any(len(s) != L for s in seqs):
        raise ValueError(f"all CDRH3 sequences must have length {L}")
    lut = np.full(256, -1, np.int16)
    for i, a in enumerate(alphabet):
        lut[ord(a)] = i
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), L)
    idx = lut[arr]
    if (idx < 0).any():
        raise ValueError("non-standard amino-acid letter in input")
    return idx.astype(np.int8)


def onehot(idx):
    n = idx.shape[0]
    X = np.zeros((n, L, 20), np.float32)
    X[np.arange(n)[:, None], np.arange(L)[None, :], idx] = 1.0
    return X


def load_facs():
    df = pd.read_csv(FACS_CSV, usecols=["cdrh3", "binding_class", "edit_distance"])
    df["label"] = (df["binding_class"] == "high").astype(int)
    return df


def load_spr():
    return pd.read_csv(SPR_CSV)


def hamming_to_wt(idx):
    wt = encode_idx([WT])[0]
    return (idx != wt[None, :]).sum(1)


# ---- "differs at >= 2 positions from every provided sequence" check ---------------------------
def _codes(idx):
    return (idx.astype(np.int64) * _PW).sum(1)


def masked_keys(idx):
    """For each position p: the sequence code with position p replaced by a wildcard digit.
    Two length-10 sequences are within Hamming distance <= 1 iff they share a masked key."""
    code = _codes(idx)
    return [code + (20 - idx[:, p].astype(np.int64)) * _PW[p] for p in range(L)]


def _pair_keys(idx):
    """Codes with two positions wildcarded, for every pair (p, q). Sequences within Hamming <= 2 share one."""
    code = _codes(idx)
    out = {}
    for p in range(L):
        for q in range(p + 1, L):
            out[(p, q)] = code + (20 - idx[:, p].astype(np.int64)) * _PW[p] + (20 - idx[:, q].astype(np.int64)) * _PW[q]
    return out


class ProvidedData:
    """Index of the provided sequences (FACS + SPR) for the novelty check."""

    def __init__(self, seqs):
        idx = encode_idx(sorted(set(seqs)))
        self._idx = idx
        self.n = len(idx)
        self.keys = [np.unique(k) for k in masked_keys(idx)]

    def has_neighbour_within_1(self, idx):
        q = masked_keys(idx)
        hit = np.zeros(len(idx), bool)
        for p in range(L):
            hit |= np.isin(q[p], self.keys[p])
        return hit

    def has_neighbour_within_2(self, idx):
        """True if a provided sequence is within 2 mutations (fails a '> 2 point mutations' rule)."""
        if not hasattr(self, "keys2"):
            self.keys2 = {pq: np.unique(k) for pq, k in _pair_keys(self._idx).items()}
        q = _pair_keys(idx)
        hit = np.zeros(len(idx), bool)
        for pq, k in q.items():
            hit |= np.isin(k, self.keys2[pq])
        return hit


def provided_sequences():
    seqs = set(load_facs()["cdrh3"])
    spr = load_spr()
    seqs |= set(s for s in spr["cdrh3"] if len(s) == L and s.isalpha())
    return seqs


# ---- the authors' CNN (Keras 2.4 Sequential) in numpy -----------------------------------------
class AuthorsCNN:
    """Conv1D(400, k=5, same, relu) -> Dropout -> MaxPool1D(2, stride 1) -> Flatten ->
    Dense(300, relu) -> Dense(1, sigmoid). Input: one-hot (10, 20). Inference only."""

    def __init__(self, h5_path):
        import h5py

        with h5py.File(h5_path, "r") as f:
            g = f["model_weights"]
            self.k = g["conv1d/conv1d/kernel:0"][()]  # (5, 20, 400)
            self.b = g["conv1d/conv1d/bias:0"][()]
            self.w1 = g["dense/dense/kernel:0"][()]  # (3600, 300)
            self.b1 = g["dense/dense/bias:0"][()]
            self.w2 = g["dense_1/dense_1/kernel:0"][()]  # (300, 1)
            self.b2 = g["dense_1/dense_1/bias:0"][()]

    def predict_proba(self, idx, batch=4096):
        out = []
        for s in range(0, len(idx), batch):
            X = onehot(idx[s : s + batch])
            pad = np.pad(X, ((0, 0), (2, 2), (0, 0)))
            conv = sum(pad[:, j : j + L, :] @ self.k[j] for j in range(5)) + self.b
            h = np.maximum(conv, 0.0)
            pooled = np.maximum(h[:, :-1, :], h[:, 1:, :])  # (n, 9, 400)
            d = np.maximum(pooled.reshape(len(X), -1) @ self.w1 + self.b1, 0.0)
            z = d @ self.w2 + self.b2
            out.append((1.0 / (1.0 + np.exp(-z)))[:, 0])
        return np.concatenate(out)
