"""Score every refolded complex against the 1N8Z crystal and compare candidates with the controls.

Reads whatever has finished under results/refold/ (JapanFold tags and the Nebius Boltz-2 runs), writes
baseline/outputs/refold_metrics.csv and prints a summary. Safe to re-run as more runs finish.
"""
import glob
import json

import numpy as np
import pandas as pd

import common as C
import structure_tools as S

BASE = C.REFOLD
SET = pd.read_csv(C.OUT / C.vname("refold_set.csv"))
rows = []


def add(service, id_, pred, conf, pae=None):
    m = S.compare(pred)
    m.update(conf)
    if pae is not None:
        m.update(S.interface_pae(pae))
    m.update(service=service, id=id_)
    rows.append(m)


# ------------------------------------------------------------------ JapanFold (per tag)
for tagdir in sorted(p for p in BASE.glob("jf_*") if p.is_dir()):
    tag = tagdir.name
    rj = BASE / f"{tag}_results.json"
    if not rj.exists():
        continue
    res = json.load(open(rj))
    for r in res["rows"]:
        id_ = r["id"]
        if r.get("status") not in (None, "succeeded", "ok"):
            continue
        pair = r.get("pair_chains_iptm") or {}
        iptm_vh = pair.get("0", {}).get("2", np.nan)
        iptm_vl = pair.get("1", {}).get("2", np.nan)
        conf = {"confidence": r.get("confidence_score"), "iptm": r.get("iptm"), "iptm_ab_her2": np.nanmean([iptm_vh, iptm_vl]),
                "plddt": r.get("complex_plddt")}
        pae = None
        npz = tagdir / f"{id_}_pae.npz"
        if npz.exists():
            pae = np.load(npz)["pae"]
        add(tag, id_, S.parse_cif((tagdir / f"{id_}.cif").read_text()), conf, pae)
        # the other samples, to see how stable the pose is
        extra = []
        for k in (1, 2):
            f = tagdir / f"{id_}_model_{k}.cif"
            if f.exists():
                extra.append(S.compare(S.parse_cif(f.read_text())))
        if extra:
            rows[-1]["contact_recall_all_samples_mean"] = float(np.mean([rows[-1]["contact_recall"]] + [e["contact_recall"] for e in extra]))

# ------------------------------------------------------------------ Nebius Boltz-2
for d in sorted(glob.glob(str(BASE / "nebius" / "*" ))):
    d = __import__("pathlib").Path(d)
    if not d.is_dir():
        continue
    files = [f for f in d.glob("*.json") if f.name != "result.json"]
    if not files:
        continue
    j = json.load(open(files[0]))
    pred = S.parse_cif(j["structures"][0]["structure"])
    add("nebius_boltz2", d.name, pred, {"confidence": j["confidence_scores"][0], "ptm": j["ptm_scores"][0]})

df = pd.DataFrame(rows).merge(SET[["id", "role", "cdrh3"]], on="id", how="left")
order = {"positive control": 0, "negative control": 1, "candidate": 2}
df["_o"] = df.role.map(order)
df = df.sort_values(["service", "_o", "id"]).drop(columns="_o")
df.to_csv(C.OUT / C.vname("refold_metrics.csv"), index=False)

show = ["id", "role", "cdrh3", "confidence", "iptm_ab_her2", "ipae_crystal_contacts", "epitope_ca_rmsd_fvframe",
        "contact_recall", "cdrh3_bb_rmsd_pose", "her2_ca_rmsd_fit", "W99_kept", "G103_kept", "F104_kept", "Y105_kept"]
for svc, g in df.groupby("service"):
    print(f"\n=== {svc} ({len(g)} complexes) ===")
    cols = [c for c in show if c in g]
    print(g[cols].round(2).to_string(index=False))
    ctrl = g[g.role != "candidate"]
    print("controls: positives vs measured non-binders (mean)")
    print(ctrl.groupby("role")[[c for c in ("confidence", "iptm_ab_her2", "ipae_crystal_contacts", "epitope_ca_rmsd_fvframe", "contact_recall") if c in ctrl]].mean().round(2).to_string())
print("\nsaved", C.OUT / C.vname("refold_metrics.csv"))
