"""Add real VHH-antigen calibration complexes from the Smorodina et al. (2026) benchmark to calib_set.csv.

Positives: cognate nanobody-antigen pairs of crystal structures (diagonal of their ipTM matrix).
Negatives: the same nanobodies against antigens they do not bind (shuffled pairs), plus the HARD shuffled pairs
that AlphaFold3 gave a high ipTM (the benchmark's point: confident but wrong). Their AF3 / Boltz-2 / Chai-1 ipTM are
stored in the note so our refolders can be compared on identical pairs.

Run: baseline/.venv/bin/python denovo/s3_score/make_vhh_calib.py [--n-pos 10] [--n-rand-neg 3] [--n-hard 12]
Needs data/smorodina/x/iptm_confidence_scores/*.csv (scp'd from the VM's hackathon-data) and internet (RCSB FASTA).
"""
import argparse
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "denovo"))
import dn_common as dn  # noqa: E402

HERE = Path(__file__).resolve().parent
M = ROOT / "data" / "smorodina" / "x" / "iptm_confidence_scores"
VHH_NAME = re.compile(r"nanobod|vhh|single.domain|sybody|\bnb\b|camel|llama|vicugna|heavy.chain.antibody", re.I)
cache = {}


def fasta(pdb):
    if pdb in cache:
        return cache[pdb]
    try:
        txt = urllib.request.urlopen(f"https://www.rcsb.org/fasta/entry/{pdb}/display", timeout=30).read().decode()
    except Exception:
        cache[pdb] = []
        return []
    ents, name = [], None
    for line in txt.splitlines():
        if line.startswith(">"):
            name = line[1:]
            ents.append([name, ""])
        elif ents:
            ents[-1][1] += line.strip()
    cache[pdb] = ents
    return ents


def split_pair(pdb):
    """(vhh_seq, antigen_seq) of one crystal entry, or None when ambiguous."""
    ents = fasta(pdb)
    vh = [e for e in ents if 100 <= len(e[1]) <= 150 and ("WGQG" in e[1] or "WGKG" in e[1] or "WGRG" in e[1] or "WGLG" in e[1])
          and (VHH_NAME.search(e[0]) or re.search(r"^[QEDK]V[QK]L", e[1]))]
    ag = [e for e in ents if e not in vh and 40 <= len(e[1]) <= 400]
    if len(vh) < 1 or len(ag) < 1:
        return None
    return vh[0][1], max(ag, key=lambda e: len(e[1]))[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-pos", type=int, default=10)
    ap.add_argument("--n-rand-neg", type=int, default=3)
    ap.add_argument("--n-hard", type=int, default=12)
    a = ap.parse_args()
    af3 = pd.read_csv(M / "af3_matrix_clean.csv", index_col=0)
    bz2 = pd.read_csv(M / "boltz2_matrix_clean.csv", index_col=0)
    ch1 = pd.read_csv(M / "chai1_matrix_clean.csv", index_col=0)
    rng = np.random.RandomState(7)
    diag = [p for p in af3.index if p in af3.columns and af3.loc[p, p] >= 0.75]
    rng.shuffle(diag)
    pairs = {}
    for p in diag:
        sp = split_pair(p)
        if sp and len(sp[1]) <= 300:
            pairs[p] = sp
        if len(pairs) >= a.n_pos:
            break
    print(f"{len(pairs)} usable cognate pairs: {list(pairs)}")
    rows = []

    def add(tag, k, nb, ag, lab, note):
        rows.append({"calib_id": f"cal-{tag}-{k:02d}", "label": lab, "group": tag, "kind": "vhh",
                     "vhh_seq": pairs[nb][0], "vl_seq": "", "antigen_seq": ag, "antigen_name": f"{tag}{k:02d}", "note": note})
    for k, p in enumerate(pairs, 1):
        add("vhhpos", k, p, pairs[p][1], 1, f"cognate {p}: AF3 {af3.loc[p, p]:.2f}, Boltz-2 {bz2.loc[p, p]:.2f}, Chai-1 {ch1.loc[p, p]:.2f}")
    k = 0
    for p in pairs:
        others = [q for q in af3.columns if q != p]
        for q in rng.choice(others, size=a.n_rand_neg, replace=False):
            sp = split_pair(q)
            if sp and len(sp[1]) <= 300:
                k += 1
                rows.append({"calib_id": f"cal-vhhneg-{k:02d}", "label": 0, "group": "vhhneg", "kind": "vhh", "vhh_seq": pairs[p][0],
                             "vl_seq": "", "antigen_seq": sp[1], "antigen_name": f"vhhneg{k:02d}",
                             "note": f"shuffled {p}->{q}: AF3 {af3.loc[p, q]:.2f}, Boltz-2 {bz2.loc[p, q]:.2f}, Chai-1 {ch1.loc[p, q]:.2f}"})
    # hard shuffled pairs for the chosen nanobodies, highest AF3 ipTM first
    cand = sorted(((af3.loc[p, q], p, q) for p in pairs for q in af3.columns if q != p and af3.loc[p, q] >= 0.5), reverse=True)
    k = 0
    for v, p, q in cand:
        if k >= a.n_hard:
            break
        sp = split_pair(q)
        if sp and len(sp[1]) <= 300:
            k += 1
            rows.append({"calib_id": f"cal-vhhhard-{k:02d}", "label": 0, "group": "vhhhard", "kind": "vhh", "vhh_seq": pairs[p][0],
                         "vl_seq": "", "antigen_seq": sp[1], "antigen_name": f"vhhhard{k:02d}",
                         "note": f"HARD shuffled {p}->{q}: AF3 {v:.2f}, Boltz-2 {bz2.loc[p, q]:.2f}, Chai-1 {ch1.loc[p, q]:.2f}"})
    f = HERE / "calib_set.csv"
    old = [r for r in dn.read_csv(f) if r["kind"] != "vhh"] if f.exists() else []
    dn.write_csv(f, old + rows, ["calib_id", "label", "group", "kind", "vhh_seq", "vl_seq", "antigen_seq", "antigen_name", "note"])
    print(f"wrote {len(rows)} VHH-antigen complexes ({sum(r['label'] for r in rows)} positive) next to {len(old)} Fv rows")


if __name__ == "__main__":
    main()
