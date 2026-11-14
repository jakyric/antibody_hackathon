"""Collect and score the cropped-HER2 templated OpenFold3 runs from the VM (workers crop/run_vm_crop.sh, outputs out/of3c_NN/<design>/seed_42/).

    denovo/.venv/bin/python denovo/vm_templated/collect_crop.py            # one pass: rsync, parse, write results/denovo/crop/of3_crop.csv
    denovo/.venv/bin/python denovo/vm_templated/collect_crop.py --loop 240 # repeat every 240 s until all queries are done (or crop/STOP_COLLECT exists)

Read only: rsyncs scores / PAE / structures, never touches VM processes. Crop: HER2 539-607 (chain C 539-580, label_seq_id 1..42; chain D 591-607,
label_seq_id 1..17; the VM CIFs keep the original HER2 numbers in auth_seq_id and use 1..N in label_seq_id, so label_seq_id is used here).
Design CIFs: chain A = HER2 489-580 (resi 1-92), B = VHH, C = HER2 591-607 (resi 1-17); HER2 residues are matched BY HER2 NUMBER (the 59
residues 539-580 and 591-607 that both models have), Kabsch fit on their CA atoms, same transform applied to the predicted VHH.
Columns (vmof3c_*): iptm, pair_iptm (VHH-HER2 chain-pair ipTM, max over the two fragments), plddt, ipsae (standard, 10 A cutoff, d0res), ipsae_chn15
(chain-length d0, 15 A cutoff, the ranking variant), min_ipae, pdockq2 (same implementation as denovo/validation/dde_pdockq2.py), pose_rmsd (VHH CA after the
HER2 fit), scrmsd (VHH-only fit), her2_fit (HER2 CA RMSD of the fit), epi_frac (share of the touched HER2 residues, heavy atoms <= 4.5 A, inside trastuzumab's epitope),
pose_ok (pose RMSD <= 10 A). Controls (ctl-*) have no design CIF, so their pose columns stay empty.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import cbor2
import numpy as np

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / "denovo" / "s3_score"))
import geom as g  # noqa: E402

CROP = ROOT / "results" / "denovo" / "crop"
RAW, MOD = CROP / "vm_raw", CROP / "vm_models"
VM = __import__("os").environ["PIPELINE_VM_HOST"]  # user@host of the shared VM
SSH = "ssh -i " + os.environ["PIPELINE_VM_KEY"] + " -o BatchMode=yes"
NC, ND = 42, 17            # HER2 residues of chain C (539-580) and chain D (591-607)
EPI = {int(r["auth_seq"]) for r in json.load(open(ROOT / "denovo/targets/epitope.json"))["residues"]}


def sync():
    for dest, pat in ((RAW, ("*_scores.json", "*_detailed_scores.cbor.gz")), (MOD, ("*_model.cif",))):
        dest.mkdir(parents=True, exist_ok=True)
        cmd = ["rsync", "-az", "-e", SSH, "--prune-empty-dirs", "--include", "of3c_*/", "--include", "of3c_*/*/", "--include", "of3c_*/*/seed_*/"]
        for p in pat:
            cmd += ["--include", p]
        cmd += ["--exclude", "*", f"{VM}:moritz/of3t/out/", str(dest) + "/"]
        subprocess.run(cmd, check=False, capture_output=True)


def arr(tag):
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 40:
        shape, data = tag.value
        return arr(data).reshape(shape)
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 85:
        return np.frombuffer(tag.value, dtype="<f4")
    return np.array(tag)


def read(path, seq="auth_seq_id"):
    hdr, atoms = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            atoms.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", seq, "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "B_iso_or_equiv")}
    return [(a[ix["label_asym_id"]], int(a[ix[seq]]), a[ix["label_atom_id"]], np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])]),
             float(a[ix["B_iso_or_equiv"]])) for a in atoms if a[ix["type_symbol"]] != "H"]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731


def pdockq2(P, pae):
    res, toks = {}, []
    for c, r, n, x, b in P:
        k = (c, r)
        if k not in res:
            res[k] = {}
            toks.append(k)
        if n == "CA":
            res[k]["ca"], res[k]["plddt"] = x, b
        if n == "CB":
            res[k]["cb"] = x
    if pae.shape[0] != len(toks) or any("ca" not in res[k] for k in toks):
        return float("nan")
    A = [i for i, k in enumerate(toks) if k[0] == "A"]
    B = [i for i, k in enumerate(toks) if k[0] != "A"]
    pos = lambda i: res[toks[i]].get("cb", res[toks[i]]["ca"])  # noqa: E731
    XA, XB = np.array([pos(i) for i in A]), np.array([pos(i) for i in B])
    d = np.linalg.norm(XA[:, None] - XB[None], axis=2)
    ii, jj = np.where(d <= 8.0)
    if not len(ii):
        return 0.005
    a, b = np.array(A)[ii], np.array(B)[jj]
    pa = (pae[a, b] + pae[b, a]) / 2.0
    plddt = np.mean([res[toks[i]]["plddt"] for i in set(a) | set(b)])
    x = plddt * np.mean(1.0 / (1.0 + (pa / 8.0) ** 2))
    return float(1.31 / (1 + np.exp(-0.075 * (x - 84.733))) + 0.005)


def design_cif(did, D):
    r = D.get(did)
    if r and r.get("design_cif") and (ROOT / r["design_cif"]).is_file():
        return ROOT / r["design_cif"]
    p = ROOT / "results/denovo/designs" / f"{did}.cif"
    return p if p.is_file() else None


def one(did, sj, D):
    out = {"design_id": did}
    s = json.load(open(sj))
    cc = s.get("chain_chain_interface_predicted_tm_score")
    out["vmof3c_iptm"] = s.get("interface_predicted_tm_score")
    out["vmof3c_plddt"] = s.get("average_predicted_local_distance_difference_test")
    if cc:
        out["vmof3c_pair_iptm"] = max(cc[0][1], cc[0][2]) if len(cc) >= 3 else cc[0][1]
    stem = sj[: -len("_scores.json")]
    cbp = stem + "_detailed_scores.cbor.gz"
    cif = str(MOD / Path(sj).relative_to(RAW).parent / (Path(stem).name + "_model.cif"))
    if not (os.path.exists(cbp) and os.path.exists(cif)):
        return None                                      # incomplete transfer: try again next pass
    pae = arr(cbor2.loads(gzip.open(cbp).read())["token_token_predicted_aligned_error"])
    n = pae.shape[0]
    nv = n - (NC + ND)
    v, h = slice(0, nv), slice(nv, n)
    cross = np.concatenate([pae[v, h].ravel(), pae[h, v].ravel()])
    out["vmof3c_min_ipae"] = float(cross.min())
    out["vmof3c_ipsae"] = g.ipsae(pae, v, h)
    out["vmof3c_ipsae_chn15"] = g.ipsae(pae, v, h, 15.0, True)
    P = read(cif, "label_seq_id")
    out["vmof3c_pdockq2"] = pdockq2(P, pae)
    ha = [((538 + r) if c == "C" else (590 + r), x) for c, r, nm, x, b in P if c in ("C", "D")]
    va = np.array([x for c, r, nm, x, b in P if c == "A"])
    touched = {r for r, x in ha if (np.linalg.norm(va - x, axis=1) <= 4.5).any()}
    out["vmof3c_epi_frac"] = round(len(touched & EPI) / max(1, len(touched)), 3)
    out["vmof3c_n_her2_contact"] = len(touched)
    dc = design_cif(did, D)
    if dc is not None:
        Dd = read(dc)
        pv = {r: x for c, r, nm, x, b in P if c == "A" and nm == "CA"}
        ph = {((538 + r) if c == "C" else (590 + r)): x for c, r, nm, x, b in P if c in ("C", "D") and nm == "CA"}
        dv = {r: x for c, r, nm, x, b in Dd if c == "B" and nm == "CA"}
        dh = {((488 + r) if c == "A" else (590 + r)): x for c, r, nm, x, b in Dd if c in ("A", "C") and nm == "CA"}
        common = sorted(set(ph) & set(dh))
        if len(common) >= 50 and len(pv) == len(dv):
            R, t = kabsch(np.array([ph[k] for k in common]), np.array([dh[k] for k in common]))
            out["vmof3c_her2_fit"] = round(rmsd(np.array([ph[k] for k in common]) @ R.T + t, np.array([dh[k] for k in common])), 2)
            vp, vd = np.array([pv[k] for k in sorted(pv)]), np.array([dv[k] for k in sorted(dv)])
            pose = rmsd(vp @ R.T + t, vd)
            R2, t2 = kabsch(vp, vd)
            out["vmof3c_pose_rmsd"] = round(pose, 1)
            out["vmof3c_scrmsd"] = round(rmsd(vp @ R2.T + t2, vd), 2)
            out["vmof3c_pose_ok"] = int(pose <= 10.0)
    return out


def collect():
    sync()
    D = {}
    for f in ("designs_rfd4jf.csv", "designs_control.csv", "designs_control_krisha.csv"):
        p = ROOT / "results/denovo" / f
        if p.exists():
            for r in csv.DictReader(open(p)):
                D[r["design_id"]] = r
    rows = {}
    for sj in sorted(glob.glob(f"{RAW}/of3c_*/*/seed_*/*_scores.json")):
        did = Path(sj).parent.parent.name
        try:
            r = one(did, sj, D)
        except Exception as e:  # noqa: BLE001
            print("skip", did, repr(e)[:100])
            continue
        if r:
            rows[did] = r
    cols = ["design_id", "vmof3c_iptm", "vmof3c_pair_iptm", "vmof3c_plddt", "vmof3c_ipsae", "vmof3c_ipsae_chn15", "vmof3c_min_ipae", "vmof3c_pdockq2",
            "vmof3c_pose_rmsd", "vmof3c_scrmsd", "vmof3c_her2_fit", "vmof3c_epi_frac", "vmof3c_pose_ok", "vmof3c_n_her2_contact"]
    with open(CROP / "of3_crop.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows.values():
            w.writerow({k: ("" if r.get(k) is None or (isinstance(r.get(k), float) and r[k] != r[k]) else r.get(k)) for k in cols})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loop", type=int, default=0, help="seconds between passes; 0 = one pass")
    a = ap.parse_args()
    total = sum(len(json.load(open(f))["queries"]) for f in glob.glob(f"{CROP}/chunks/q_*.json"))
    while True:
        rows = collect()
        he = [r["vmof3c_her2_fit"] for r in rows.values() if "vmof3c_her2_fit" in r]
        print(time.strftime("%H:%M:%S"), f"{len(rows)}/{total} designs parsed | HER2 fit median {np.median(he):.2f} A (n={len(he)})" if he else f"{len(rows)}/{total} parsed", flush=True)
        if not a.loop or len(rows) >= total or (CROP / "STOP_COLLECT").exists():
            break
        time.sleep(a.loop)


if __name__ == "__main__":
    main()
