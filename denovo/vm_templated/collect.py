"""Pull the VM's templated OpenFold3 (HER2-only template, no MSA) and ESMFold-2 outputs and write results/denovo/score_vmof3.csv.

    denovo/.venv/bin/python denovo/vm_templated/collect.py        # rsyncs scores only (no structures), parses, writes the CSV
Columns: vmof3_* from OpenFold3 (iptm = interface pTM of the whole complex, pair_iptm = VHH-HER2 chain-pair ipTM (max over the two HER2
fragments), ptm, plddt, min_ipae / mean_ipae over the VHH-HER2 interface, ipsae_chn15 = ipSAE with the chain-length d0 and a 15 A PAE
cutoff as in S3), vmesm2_* from ESMFold-2 (same fields where its output has them). Keyed by design_id; controls have ctl- ids.
"""
import glob
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import cbor2
import numpy as np

sys.path.insert(0, "denovo")
import dn_common as dn  # noqa: E402

RAW = dn.OUT / "vm_raw"
RAW.mkdir(exist_ok=True)
subprocess.run(["rsync", "-az", "-e", "ssh -i " + os.environ["PIPELINE_VM_KEY"] + " -o BatchMode=yes", "--include", "*/", "--include", "*scores.json",
                "--include", "*detailed_scores.cbor.gz", "--exclude", "*", os.environ["PIPELINE_VM_HOST"] + ":moritz/of3t/out/", str(RAW) + "/"],
               check=False)


def arr(tag):
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 40:
        shape, data = tag.value
        return arr(data).reshape(shape)
    if isinstance(tag, cbor2.CBORTag) and tag.tag == 85:
        return np.frombuffer(tag.value, dtype="<f4")
    return np.array(tag)


def d0(L):
    return max(1.0, 1.24 * (L - 15.0) ** (1.0 / 3.0) - 1.8)


def ipsae(pae, ia, ib, cutoff=15.0):
    d0c = d0((ia.stop - ia.start) + (ib.stop - ib.start))
    def one(M):
        best = 0.0
        for i in range(M.shape[0]):
            ok = M[i] < cutoff
            if ok.any():
                best = max(best, float((1.0 / (1.0 + (M[i][ok] / d0c) ** 2)).mean()))
        return best
    return max(one(pae[ia, ib]), one(pae[ib, ia].T))


rows = {}
for model in ("of3", "esm2"):
    for sj in glob.glob(f"{RAW}/{model}_*/*/seed_*/*_scores.json"):
        did = Path(sj).parent.parent.name
        try:
            s = json.load(open(sj))
            pre = "vmof3_" if model == "of3" else "vmesm2_"
            r = rows.setdefault(did, {"design_id": did})
            cc = s.get("chain_chain_interface_predicted_tm_score")
            r[pre + "iptm"] = s.get("interface_predicted_tm_score")
            r[pre + "ptm"] = s.get("predicted_tm_score")
            r[pre + "plddt"] = s.get("average_predicted_local_distance_difference_test")
            if cc:
                r[pre + "pair_iptm"] = max(cc[0][1], cc[0][2]) if len(cc) >= 3 else cc[0][1]
            cb = sj.replace("_scores.json", "_detailed_scores.cbor.gz")
            if os.path.exists(cb):
                d = cbor2.loads(gzip.open(cb).read())
                pae = arr(d["token_token_predicted_aligned_error"])
                n = pae.shape[0]
                nv = n - (92 + 17)  # VHH tokens first, then HER2 489-580 (92) and 591-607 (17)
                ia, ib = slice(0, nv), slice(nv, n)
                inter = np.concatenate([pae[ia, ib].ravel(), pae[ib, ia].ravel()])
                r[pre + "min_ipae"] = float(inter.min())
                r[pre + "mean_ipae"] = float(inter.mean())
                r[pre + "ipsae_chn15"] = ipsae(pae, ia, ib, 15.0)
        except Exception as e:  # noqa: BLE001
            print("skip", sj, e)
out = list(rows.values())
cols = ["design_id"] + sorted({k for r in out for k in r if k != "design_id"})
dn.write_csv(dn.OUT / "score_vmof3.csv", out, cols)
print(f"{len(out)} designs parsed ->", dn.OUT / "score_vmof3.csv",
      "| OF3:", sum("vmof3_iptm" in r for r in out), "| ESMFold-2:", sum("vmesm2_iptm" in r for r in out))
