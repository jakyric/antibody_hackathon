"""Why does OpenDDE recover the pose for Jakob's redesigned trastuzumab Fvs but not for our de novo VHHs?  (user, 16:5x)

Jakob's branch (jakyric/antibody_hackathon, `jakob`): trastuzumab Fv (VH + VL) with a new CDR-H3 and Vernier swaps, folded with JapanFold
opendde-abag; 201 of 517 designs have ipTM >= 0.9 and every one of them has the HER2 epitope within 3 A (median 1.2 A) of the 1N8Z crystal after fitting
the Fv framework. Prediction settings are the same model and the same HER2 construct as ours (domain IV, 119 aa incl. 581-590), so two things differ:
(1) input type: paired Fv (A = VH, B = VL, C = HER2) versus our single VHH, (2) sequence novelty: >95 % identical to the 1N8Z antibody versus our all-CDR de novo VHHs.
This script isolates (1) and the settings: Jakob's top 20 candidates are folded with OUR OpenDDE parameters (MODELS['dde'] in run_s3.py: MSA server on, seed 0,
3 diffusion samples, cif) as (a) Fv + HER2 and (b) VH only + HER2 (single chain like our VHHs).
Pose check = Kabsch on the predicted HER2 CA atoms (matched by HER2 number) onto the 1N8Z HER2, then CA RMSD of the predicted VH onto the crystal VH.
Run: baseline/.venv/bin/python denovo/validation/jakob_check.py <path to Jakob's top_candidates.json>  -> results/denovo/jakob_check.csv
"""
import hashlib
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "denovo")
import dn_common as dn  # noqa: E402

cands = json.load(open(sys.argv[1]))["candidates"]
HER2 = open("baseline/refold_inputs/her2_domIV.txt").read().strip()
PARAMS = dict(write_pae=False, seed=0, output_format="cif", diffusion_samples=3, use_msa_server=True)
yaml = lambda ch: "sequences:\n" + "".join(f"  - protein: {{id: {c}, sequence: {s}}}\n" for c, s in ch)  # noqa: E731
targets, kind = [], {}
for c in cands:
    i = c["id"]
    targets.append({"name": f"{i}__fv", "content": yaml([("A", c["vh_sequence"]), ("B", c["vl_sequence"]), ("C", HER2)])})
    targets.append({"name": f"{i}__vh", "content": yaml([("A", c["vh_sequence"]), ("B", HER2)])})
name = "jakobcheck-" + hashlib.sha1(",".join(t["name"] for t in targets).encode()).hexdigest()[:10]
jf = dn.JapanFold()
jid = jf.submit(name, "/v1/predictions", {"model": "opendde-abag", "name": name, "targets": targets, "params": PARAMS})
print("job", jid, len(targets), "targets", flush=True)
j = jf.wait(jid, max_s=3000, log=lambda s: None)
assert j["status"] == "succeeded", j
res = jf.results(jid)
dest = dn.OUT / "jakob_check" / jid
dest.mkdir(parents=True, exist_ok=True)
jf.download(jid, dest)


def read(path, seq="auth_seq_id"):
    hdr, atoms = [], []
    for l in open(path):
        if l.startswith("_atom_site."):
            hdr.append(l.strip().split(".")[1])
        elif l.startswith(("ATOM", "HETATM")):
            atoms.append(l.split())
    ix = {h: hdr.index(h) for h in ("label_atom_id", "label_asym_id", seq, "Cartn_x", "Cartn_y", "Cartn_z")}
    return [(a[ix["label_asym_id"]], int(a[ix[seq]]), np.array([float(a[ix["Cartn_x"]]), float(a[ix["Cartn_y"]]), float(a[ix["Cartn_z"]])]))
            for a in atoms if a[ix["label_atom_id"]] == "CA"]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


rmsd = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))  # noqa: E731
X = read("denovo/targets/ref_1n8z.cif")                          # crystal: A = VL, B = VH, C/D = HER2 (original numbering)
xh = {r: x for c, r, x in X if c in ("C", "D")}
xv = np.array([x for c, r, x in sorted((t for t in X if t[0] == "B"), key=lambda t: t[1])])
rows_by_id = {r["id"]: r for r in res["rows"]}
out = []
for t in targets:
    tid = t["name"]
    row = rows_by_id.get(tid, {})
    cif = dest / f"{tid}.cif"
    if not cif.exists():
        out.append({"id": tid, "status": row.get("status", "missing")})
        continue
    P = read(cif, "label_seq_id")
    her_ch = "C" if tid.endswith("__fv") else "B"
    ph = {488 + r: x for c, r, x in P if c == her_ch}
    common = sorted(set(ph) & set(xh))
    R, tt = kabsch(np.array([ph[k] for k in common]), np.array([xh[k] for k in common]))
    pv = np.array([x for c, r, x in sorted((q for q in P if q[0] == "A"), key=lambda q: q[1])])
    pose = rmsd(pv @ R.T + tt, xv) if len(pv) == len(xv) else float("nan")
    out.append({"id": tid, "status": "ok", "iptm": row.get("iptm"), "complex_plddt": row.get("complex_plddt"),
                "her2_fit": round(rmsd(np.array([ph[k] for k in common]) @ R.T + tt, np.array([xh[k] for k in common])), 2), "vh_pose_rmsd": round(pose, 1)})
D = pd.DataFrame(out)
D["input"] = D.id.str.rsplit("__", n=1).str[-1]
D.to_csv("results/denovo/jakob_check.csv", index=False)
for k, g in D.groupby("input"):
    g = g[g.status == "ok"]
    print(f"{k}: n={len(g)} | ipTM median {pd.to_numeric(g.iptm).median():.2f} | VH pose RMSD to the crystal pose: median {g.vh_pose_rmsd.median():.1f} A, "
          f"<= 3 A: {int((g.vh_pose_rmsd <= 3).sum())}, <= 10 A: {int((g.vh_pose_rmsd <= 10).sum())} | HER2 fit median {g.her2_fit.median():.2f} A")
