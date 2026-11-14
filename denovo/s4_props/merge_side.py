"""Append the outputs of one `run_vm.sh` run (nanomelt.csv, tnp.csv, my_monoyield.csv) to the side files that rank.py reads
(results/denovo/props_nm.csv, props_tnpvm.csv, props_my.csv). Only adds designs that are not in the side file yet, never rewrites existing rows.

    denovo/.venv/bin/python denovo/s4_props/merge_side.py results/denovo/props_work/final_fill2/out
"""
import sys
from pathlib import Path

import pandas as pd

run = Path(sys.argv[1])
out = Path("results/denovo")


def add(side, new, key="design_id"):
    p = out / side
    old = pd.read_csv(p) if p.exists() else pd.DataFrame(columns=[key])
    new = new[~new[key].isin(set(old[key]))]
    pd.concat([old, new], ignore_index=True).to_csv(p, index=False)
    print(f"{side}: +{len(new)} rows (now {len(old) + len(new)})")


f = run / "nanomelt.csv"
if f.exists() and f.stat().st_size:
    t = pd.read_csv(f)
    tm = [c for c in t.columns if "Tm" in c][0]
    add("props_nm.csv", pd.DataFrame({"design_id": t["ID"], "nm_tm": t[tm].round(2)}).dropna())
f = run / "tnp.csv"
if f.exists() and f.stat().st_size:
    t = pd.read_csv(f)
    if "design_id" not in t.columns and "name" in t.columns:
        t = t.rename(columns={"name": "design_id"})
    keep = [c for c in pd.read_csv(out / "props_tnpvm.csv", nrows=1).columns]
    t = t[[c for c in keep if c in t.columns]]
    add("props_tnpvm.csv", t[t["tnp_n_red"].notna()])  # rows where TNP failed have no flags; they stay missing
f = run / "my_monoyield.csv"
if f.exists() and f.stat().st_size:
    t = pd.read_csv(f)
    add("props_my.csv", pd.DataFrame({"design_id": t["name"], "my_yield": t["pred_mean"], "my_yield_sd": t["pred_std"]}))
