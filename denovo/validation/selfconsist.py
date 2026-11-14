"""VAL check 8: self-consistency. Predict the VHH ALONE with NanoBodyBuilder2 (ImmuneBuilder, CPU, micromamba env `tnp`) for a
stratified sample of designs and compare with the RFD4 design's VHH: framework CA RMSD and per-CDR backbone RMSD after
superposing on the framework. Also compares S3's OpenFold3 pose (VHH aligned on HER2) with the design pose where available.

    denovo/.venv/bin/python denovo/validation/selfconsist.py [--n 14] [--workers 6]
Writes results/denovo/validity_selfconsist.csv and denovo/validation/nb2/*.pdb.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import validate as V  # noqa: E402

OUT = V.OUT
NB2 = HERE / "nb2"
NB2.mkdir(exist_ok=True)

PRED = """
import sys, os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import torch
torch.set_num_threads(2)
from ImmuneBuilder import NanoBodyBuilder2
seq, out = sys.argv[1], sys.argv[2]
m = NanoBodyBuilder2()
nb = m.predict({"H": seq})
try:
    nb.save_single_unrefined(out)
except Exception:
    nb.save(out)
print("saved", out)
"""


def pick(n):
    D = pd.read_csv(OUT / "designs_rfd4jf.csv")
    D["pose"] = D.params_json.map(lambda s: json.loads(s)["pose"])
    Vd = pd.read_csv(OUT / "validity.csv")
    M = D.merge(Vd[["design_id", "clash_intra_vhh", "clash_vhh_her2_lt2p2"]], on="design_id")
    M["bad"] = M.clash_intra_vhh + M.clash_vhh_her2_lt2p2
    chosen = []
    try:
        P = pd.read_csv(OUT / "leaderboard_pareto.csv")
        P = P[P.design_id.isin(M.design_id)].sort_values(["pareto_rank", "maximin"], ascending=[True, False])
        for fw in ("bcii10", "capla"):
            chosen += list(P[P.framework == fw].design_id.head(3))
        L = pd.read_csv(OUT / "leaderboard.csv")
        L = L[L.design_id.isin(M.design_id) & (L.is_control == 0) & L.composite.notna()].sort_values("composite")
        mid = len(L) // 2
        chosen += list(L.design_id.iloc[mid - 1:mid + 1])
    except Exception as e:  # noqa: BLE001
        print("pareto/leaderboard sample failed:", e)
    chosen += list(M.sort_values("bad", ascending=False).design_id.head(3))
    for fw in ("bcii10", "capla"):
        chosen += list(M[(M.framework == fw) & (M.pose == "p1")].design_id.head(2))
    seen, out = set(), []
    for c in chosen:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return M[M.design_id.isin(out[:n])]


def run_one(row):
    did, seq = row["design_id"], row["vhh_seq"]
    pdb = NB2 / f"{did}.pdb"
    if not pdb.exists():
        p = subprocess.run(["micromamba", "run", "-n", "tnp", "python", "-c", PRED, seq, str(pdb)], capture_output=True, text=True,
                           env={**os.environ, "CUDA_VISIBLE_DEVICES": ""}, timeout=900)
        if p.returncode != 0 or not pdb.exists():
            return {"design_id": did, "nb2_error": (p.stderr or p.stdout)[-200:]}
    return analyse(row, pdb)


def read_pdb(path, seq=None):
    """Residue list in SEQUENCE order. ImmuneBuilder writes IMGT numbers with insertion codes, and IMGT lists some insertion
    runs in descending order (111, 111A, 111B, 112B, 112A, 112), so file order is not sequence order: reorder each same-number
    group (forward or reversed) to match the design sequence."""
    res = {}
    for line in open(path):
        if line.startswith("ATOM"):
            key = (int(line[22:26]), line[26])
            r = res.setdefault(key, {"aa": V.THREE.get(line[17:20], "X"), "num": key[0], "atoms": {}})
            r["atoms"][line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    items = list(res.values())
    if seq is None:
        return [r["atoms"] for r in items]
    out, i, off = [], 0, 0
    while i < len(items):
        j = i
        while j < len(items) and items[j]["num"] == items[i]["num"]:
            j += 1
        grp = items[i:j]
        for cand in (grp, grp[::-1]):
            if "".join(r["aa"] for r in cand) == seq[off:off + len(cand)]:
                grp = cand
                break
        out += grp
        off += len(grp)
        i = j
    return [r["atoms"] for r in out]


def analyse(row, pdb):
    P = json.loads(row["params_json"])
    caa = bool(P.get("caa_fixed", 0))
    ch, _, _ = V.read_cif(V.ROOT / row["design_cif"])
    B = ch["B"]
    seq = V.seq_of(B)
    fr = V.FW[row["framework"]]["caa"]["fr"] if caa else V.FW[row["framework"]]["fr"]
    kept, bare = V.parse_contig(P["contig"])
    ok, st, cdrs, why = V.derive_segments(seq, fr, bare)
    pred = read_pdb(pdb, seq)
    out = {"design_id": row["design_id"], "framework": row["framework"], "pose": P["pose"], "nb2_n": len(pred), "design_n": len(B)}
    if not ok or len(pred) != len(B):
        out["nb2_error"] = f"length mismatch {len(pred)} vs {len(B)}" if ok else why
        return out
    idx = {"fw": [], "cdr1": [], "cdr2": [], "cdr3": []}
    spans = {"cdr1": (len(fr[0]), st[1]), "cdr2": (st[1] + len(fr[1]), st[2]), "cdr3": (st[2] + len(fr[2]), st[3])}
    for i in range(len(B)):
        k = "fw"
        for n, (a, b) in spans.items():
            if a <= i < b:
                k = n
        idx[k].append(i)
    dca = np.array([r["atoms"]["CA"] for r in B])
    pca = np.array([r["CA"] for r in pred])
    f = idx["fw"]
    R, t = V.kabsch(pca[f], dca[f])
    out["nb2_fw_ca_rmsd"] = V.rmsd((R @ pca[f].T).T + t, dca[f])
    for k in ("cdr1", "cdr2", "cdr3"):
        i = idx[k]
        out[f"nb2_{k}_ca_rmsd_on_fw"] = V.rmsd((R @ pca[i].T).T + t, dca[i])
        bd = np.array([[B[j]["atoms"][a] for a in V.BB] for j in i]).reshape(-1, 3)
        bp = np.array([[pred[j][a] for a in V.BB] for j in i]).reshape(-1, 3)
        out[f"nb2_{k}_bb_rmsd_on_fw"] = V.rmsd((R @ bp.T).T + t, bd)
    # OpenFold3 pose (S3), VHH aligned on HER2, where available
    try:
        S = pd.read_csv(OUT / "score_of3.csv")
        s = S[S.design_id == row["design_id"]]
        if len(s) and isinstance(s.iloc[0].get("of3_cif"), str):
            oc, _, _ = V.read_cif(V.ROOT / s.iloc[0]["of3_cif"])
            ov, oh = oc.get("A", []), oc.get("B", [])
            if len(ov) == len(B) and oh:
                HA = ch["A"] + ch["C"]
                # S3 folds the 119 aa domain: map by residue number onto the modelled 109
                d_num = {488 + r["num"] if k < 92 else 590 + (k - 92) + 1: r for k, r in enumerate(HA)}
                pairs = [(r, d_num[488 + r["num"]]) for r in oh if (488 + r["num"]) in d_num and "CA" in r["atoms"]]
                if len(pairs) > 50:
                    X = np.array([p[0]["atoms"]["CA"] for p in pairs])
                    Y = np.array([p[1]["atoms"]["CA"] for p in pairs])
                    R2, t2 = V.kabsch(X, Y)
                    ovca = (R2 @ np.array([r["atoms"]["CA"] for r in ov]).T).T + t2
                    out["of3_her2_fit_rmsd"] = V.rmsd((R2 @ X.T).T + t2, Y)
                    out["of3_fw_ca_rmsd_on_her2"] = V.rmsd(ovca[f], dca[f])
                    out["of3_cdr3_ca_rmsd_on_her2"] = V.rmsd(ovca[idx["cdr3"]], dca[idx["cdr3"]])
                    out["of3_cdr3_centroid_dist"] = float(np.linalg.norm(ovca[idx["cdr3"]].mean(0) - dca[idx["cdr3"]].mean(0)))
    except Exception as e:  # noqa: BLE001
        out["of3_note"] = str(e)[:80]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=14)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    sel = pick(a.n)
    print("sample:", list(sel.design_id))
    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(run_one, sel.to_dict("records")))
    pd.DataFrame(res).to_csv(OUT / "validity_selfconsist.csv", index=False)
    print(pd.DataFrame(res).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
