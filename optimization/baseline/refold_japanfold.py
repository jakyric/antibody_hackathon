"""Submit the refolding set to JapanFold as one batch job, wait, and download the results.

  optimization/baseline/.venv/bin/python optimization/baseline/refold_japanfold.py --model opendde-abag --msa on  --tag jf_abag_msa
  optimization/baseline/.venv/bin/python optimization/baseline/refold_japanfold.py --model openfold3    --msa on  --tag jf_of3_msa
  optimization/baseline/.venv/bin/python optimization/baseline/refold_japanfold.py --model opendde-abag --msa off --tag jf_abag_nomsa

The API key is read from JAPANFOLD_API_KEY and never printed. Re-running a tag reuses its saved job id instead of
submitting (and paying for) a second job.
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import uuid
import zipfile
from pathlib import Path

import pandas as pd

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True)
ap.add_argument("--msa", choices=["on", "off"], default="on")
ap.add_argument("--tag", required=True)
ap.add_argument("--samples", type=int, default=3)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--ids", nargs="*", default=None, help="only submit these complex ids")
a = ap.parse_args()

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY = os.environ["JAPANFOLD_API_KEY"]
HDR = {"X-JapanFold-Client": "skill", "Authorization": f"Bearer {KEY}"}
OUTDIR = C.REFOLD
OUTDIR.mkdir(parents=True, exist_ok=True)
idfile = OUTDIR / f"{a.tag}.jobid"


def call(method, path, body=None, headers=None, timeout=120):
    h = dict(HDR)
    if body is not None:
        h["Content-Type"] = "application/json"
    h.update(headers or {})
    req = urllib.request.Request(BASE + path, data=json.dumps(body).encode() if body is not None else None,
                                 headers=h, method=method)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read(), r.headers


if idfile.exists():
    job_id = idfile.read_text().strip()
    print(f"reusing job {job_id} for tag {a.tag}")
else:
    S = pd.read_csv(C.OUT / C.vname("refold_set.csv"))
    if a.ids:
        S = S[S.id.isin(a.ids)]
    targets = [{"name": r.id, "content": (C.JF_IN / f"{r.id}.yaml").read_text()}
               for r in S.itertuples()]
    params = {"use_msa_server": a.msa == "on", "diffusion_samples": a.samples, "write_pae": True, "seed": a.seed,
              "output_format": "cif"}
    body = {"model": a.model, "name": a.tag, "targets": targets, "params": params}
    raw, _ = call("POST", "/v1/predictions", body, {"Idempotency-Key": str(uuid.uuid4())})
    job = json.loads(raw)
    if "id" not in job:
        sys.exit(f"submit failed: {job.get('detail') or job}")
    job_id = job["id"]
    idfile.write_text(job_id)
    print(f"submitted {a.tag}: job {job_id} ({len(targets)} targets, model {a.model}, MSA {a.msa}, samples {a.samples})")
    if job.get("warnings"):
        print("warnings:", json.dumps(job["warnings"])[:600])

last = None
while True:
    raw, _ = call("GET", f"/v1/jobs/{job_id}", headers={"Prefer": "wait=60"}, timeout=150)
    j = json.loads(raw)
    st = j.get("status")
    if st != last:
        print(f"[{a.tag}] status: {st}", (f"queue: {json.dumps(j['queue'])[:200]}" if j.get("queue") else ""), flush=True)
        last = st
    if st in ("succeeded", "failed", "canceled"):
        break
    time.sleep(3)
if st != "succeeded":
    sys.exit(f"job ended {st}: {json.dumps(j)[:600]}")

raw, _ = call("GET", f"/v1/jobs/{job_id}/results")
(OUTDIR / f"{a.tag}_results.json").write_bytes(raw)
res = json.loads(raw)
if res.get("partial"):
    print("PARTIAL:", json.dumps(res["partial"]))
raw, _ = call("GET", f"/v1/jobs/{job_id}/archive", timeout=300)
zpath = OUTDIR / f"{a.tag}.zip"
zpath.write_bytes(raw)
dest = OUTDIR / a.tag
dest.mkdir(exist_ok=True)
with zipfile.ZipFile(zpath) as z:
    for m in z.infolist():
        t = (dest / m.filename).resolve()
        if not str(t).startswith(str(dest.resolve()) + os.sep) and t != dest.resolve():
            sys.exit(f"refusing unsafe path in archive: {m.filename}")
    z.extractall(dest)
print(f"[{a.tag}] done: {len(list(dest.rglob('*')))} files in {dest}; results json saved")
