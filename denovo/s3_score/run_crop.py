"""OpenDDE-abag on the CROPPED HER2 (539-607, 69 aa contiguous) for the designs and controls in results/denovo/crop/set.csv.

Why: on the full 119-aa domain the models dock the VHH on the far side of HER2 (median 44 A from the designed pose, 0% epitope contact).
Cropping away 489-538 (the first ~50 residues; 539 is the first residue that severs no disulfide, all epitope residues 557-605 are kept) limits
where the model can place the VHH. Standalone: does not touch run_s3.py or the shared S3 state.

    denovo/.venv/bin/python denovo/s3_score/run_crop.py            # submit (cached by job name), wait, download, parse
    denovo/.venv/bin/python denovo/s3_score/run_crop.py --parse    # only parse what is already downloaded

Target name = <design_id>__her2c. Model opendde-abag, use_msa_server, write_pae, seed 0, cif, diffusion_samples 3 (as MODELS['dde'] in run_s3.py).
Predicted complex: chain A = VHH, chain B = HER2 539-607 (label resi k = HER2 538+k). Designed complex (RFD4 CIF): chain A = HER2 489-580
(resi 1-92), chain B = VHH, chain C = HER2 591-607 (resi 1-17); HER2 residues are matched BY HER2 NUMBER (59 common CA atoms).
Output: results/denovo/crop/dde_crop.csv (columns dde_c_*), raw files results/denovo/crop/s3_raw/ddec/<job>/.
Metrics (sample 0): iptm, ptm, plddt (JapanFold row), ipsae (standard d0res, 10 A), ipsae_chn15 (chain-length d0, 15 A), pdockq2 (same implementation
as denovo/validation/dde_pdockq2.py), pose_rmsd (VHH CA after fitting the predicted HER2 on the designed HER2), scrmsd (VHH-only fit), her2_fit,
epi_frac / n_contact (share of the HER2 residues touched by the VHH, heavy atoms <= 4.5 A, inside trastuzumab's epitope; number of touched residues).
"""
import argparse
import concurrent.futures as cf
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "denovo"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dn_common as dn  # noqa: E402
import geom as g  # noqa: E402

CROP = ROOT / "results" / "denovo" / "crop"
RAW = CROP / "s3_raw" / "ddec"
CUT = 539
HER2C = (CROP / "her2_crop.txt").read_text().strip()
assert len(HER2C) == 69
PARAMS = dict(write_pae=True, seed=0, output_format="cif", diffusion_samples=3, use_msa_server=True)
MODEL = "opendde-abag"
BATCH = 40
EPI = {int(r["auth_seq"]) for r in json.load(open(ROOT / "denovo/targets/epitope.json"))["residues"]}


def sequences():
    d = {}
    for f in ("designs_rfd4jf.csv", "designs_control.csv", "designs_control_krisha.csv"):
        for r in csv.DictReader(open(ROOT / "results/denovo" / f)):
            d[r["design_id"]] = r
    return d


def tid(i):
    return f"{i}__her2c"


# ------------------------------------------------------------------ jobs
def run_job(jf, ids, seqs, n):
    names = sorted(tid(i) for i in ids)
    jname = "crop-dde-" + hashlib.sha1((json.dumps(PARAMS, sort_keys=True) + ",".join(names)).encode()).hexdigest()[:10]
    body = {"model": MODEL, "name": jname, "params": PARAMS,
            "targets": [{"name": tid(i), "content": f"sequences:\n  - protein: {{id: A, sequence: {seqs[i]['vhh_seq']}}}\n  - protein: {{id: B, sequence: {HER2C}}}\n"}
                        for i in ids]}
    for attempt in range(5):
        try:
            jid = jf.submit(jname, "/v1/predictions", body)
            break
        except RuntimeError as e:
            print(f"[job{n}] submit refused ({str(e)[:100]}), retry {attempt + 1}/5 in 30 s", flush=True)
            time.sleep(30)
    else:
        return None
    print(f"[job{n}] {jname}: {jid}, {len(ids)} targets", flush=True)
    j = jf.wait(jid, max_s=3300, log=lambda s: print(f"[job{n}] {s}", flush=True))
    if j["status"] != "succeeded":
        print(f"[job{n}] ended {j['status']}: {json.dumps(j)[:300]}", flush=True)
        return None
    dest = RAW / jid
    dest.mkdir(parents=True, exist_ok=True)
    res = jf.results(jid)
    (dest / "s3_results.json").write_text(json.dumps(res))
    if not any(dest.glob("*.cif")):
        jf.download(jid, dest)
    print(f"[job{n}] downloaded {len(list(dest.glob('*__her2c.cif')))} structures", flush=True)
    return jid


# ------------------------------------------------------------------ metrics
def read(path):
    hdr, atoms = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            atoms.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "B_iso_or_equiv")}
    return [(a[ix["label_asym_id"]], int(a[ix["auth_seq_id"]]), a[ix["label_atom_id"]], np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])]),
             float(a[ix["B_iso_or_equiv"]])) for a in atoms if a[ix["type_symbol"]] != "H"]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731


def pdockq2(P, pae, n_vhh):
    res, order = {}, []
    for c, r, n, x, b in P:
        k = (c, r)
        if k not in res:
            res[k] = {}
            order.append(k)
        if n == "CA":
            res[k]["ca"], res[k]["plddt"] = x, b
        if n == "CB":
            res[k]["cb"] = x
    if pae.shape[0] != len(order) or any("ca" not in res[k] for k in order):
        return float("nan")
    A = [i for i, k in enumerate(order) if k[0] == "A"]
    B = [i for i, k in enumerate(order) if k[0] == "B"]
    pos = lambda i: res[order[i]].get("cb", res[order[i]]["ca"])  # noqa: E731
    d = np.linalg.norm(np.array([pos(i) for i in A])[:, None] - np.array([pos(j) for j in B])[None], axis=2)
    ii, jj = np.where(d <= 8.0)
    if not len(ii):
        return 0.005
    a, b = np.array(A)[ii], np.array(B)[jj]
    pa = (pae[a, b] + pae[b, a]) / 2.0
    plddt = np.mean([res[order[i]]["plddt"] for i in set(a) | set(b)])
    x = plddt * np.mean(1.0 / (1.0 + (pa / 8.0) ** 2))
    return float(1.31 / (1 + np.exp(-0.075 * (x - 84.733))) + 0.005)


def design_cif(i, seqs):
    c = seqs[i].get("design_cif") or ""
    for p in (ROOT / c if c else None, ROOT / "results/denovo/designs" / f"{i}.cif"):
        if p is not None and p.is_file():
            return p
    return None


def parse(i, row, dest, seqs):
    t = tid(i)
    m = {"design_id": i, "dde_c_iptm": row.get("iptm"), "dde_c_ptm": row.get("ptm"), "dde_c_plddt": row.get("complex_plddt")}
    cif, npz = dest / f"{t}.cif", dest / f"{t}_pae.npz"
    P = read(cif)
    nv, nh = len(seqs[i]["vhh_seq"]), len(HER2C)
    if npz.exists():
        pae = np.load(npz)["pae"]
        if pae.shape[0] == nv + nh:
            v, h = slice(0, nv), slice(nv, nv + nh)
            m["dde_c_ipsae"] = g.ipsae(pae, v, h)
            m["dde_c_ipsae_chn15"] = g.ipsae(pae, v, h, 15.0, True)
            m["dde_c_pdockq2"] = pdockq2(P, pae, nv)
    pv = {r: x for c, r, n, x, b in P if c == "A" and n == "CA"}
    ph = {CUT - 1 + r: x for c, r, n, x, b in P if c == "B" and n == "CA"}
    va = np.array([x for c, r, n, x, b in P if c == "A"])
    touched = {CUT - 1 + r for c, r, n, x, b in P if c == "B" and (np.linalg.norm(va - x, axis=1) <= 4.5).any()}
    m["dde_c_n_contact"] = len(touched)
    m["dde_c_epi_frac"] = round(len(touched & EPI) / max(1, len(touched)), 3)
    dc = design_cif(i, seqs) if not i.startswith("ctl-") else None
    if dc is not None:
        D = read(dc)
        dv = {r: x for c, r, n, x, b in D if c == "B" and n == "CA"}
        dh = {(488 + r if c == "A" else 590 + r): x for c, r, n, x, b in D if c in ("A", "C") and n == "CA"}
        common = sorted(set(ph) & set(dh))
        if len(common) >= 50 and len(pv) == len(dv):
            R, tt = kabsch(np.array([ph[k] for k in common]), np.array([dh[k] for k in common]))
            m["dde_c_her2_fit"] = rmsd(np.array([ph[k] for k in common]) @ R.T + tt, np.array([dh[k] for k in common]))
            vp, vd = np.array([pv[k] for k in sorted(pv)]), np.array([dv[k] for k in sorted(dv)])
            m["dde_c_pose_rmsd"] = rmsd(vp @ R.T + tt, vd)
            R2, t2 = kabsch(vp, vd)
            m["dde_c_scrmsd"] = rmsd(vp @ R2.T + t2, vd)
            m["dde_c_pose_ok"] = int(m["dde_c_pose_rmsd"] <= 10.0)
    return m


def parse_all(seqs, ids):
    want = {tid(i): i for i in ids}
    rows = {}
    for f in sorted(RAW.glob("*/s3_results.json"), key=os.path.getmtime):   # later jobs overwrite earlier ones
        res = json.loads(f.read_text())
        for r in res.get("rows", []):
            i = want.get(r.get("id"))
            if i is None or r.get("status") != "ok":
                continue
            try:
                rows[i] = dict(parse(i, r, f.parent, seqs), dde_c_job=f.parent.name)
            except Exception as e:  # noqa: BLE001
                print("parse failed", i, repr(e)[:120], flush=True)
    df = pd.DataFrame(list(rows.values()))
    if len(df):
        df = df.round(4)
        df.to_csv(CROP / "dde_crop.csv", index=False)
    print(f"parsed {len(df)} of {len(ids)} -> {CROP / 'dde_crop.csv'}", flush=True)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parse", action="store_true")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    seqs = sequences()
    ids = [i for i in pd.read_csv(CROP / "set.csv").design_id if i in seqs]
    RAW.mkdir(parents=True, exist_ok=True)
    if not a.parse:
        jf = dn.JapanFold()
        c0 = jf.credits()
        batches = [ids[k:k + BATCH] for k in range(0, len(ids), BATCH)]
        print(f"{len(ids)} targets in {len(batches)} jobs; credits ${c0:.3f}", flush=True)
        with cf.ThreadPoolExecutor(max_workers=min(a.workers, len(batches))) as ex:
            for f in [ex.submit(run_job, jf, b, seqs, n) for n, b in enumerate(batches)]:
                try:
                    f.result()
                except Exception as e:  # noqa: BLE001
                    print("job error:", repr(e)[:300], flush=True)
        c1 = jf.credits()
        print(f"credits ${c1:.3f} (spent ${c0 - c1:.3f})", flush=True)
    parse_all(seqs, ids)


if __name__ == "__main__":
    main()
