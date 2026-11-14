"""Build the Nebius protenix-v2 input (native Protenix JSON list of jobs) for the top Pareto designs + controls.

    baseline/.venv/bin/python denovo/s3d_protenix/build_input.py [--cap 40]
Writes nebius-runs/ptxtop-input.json and denovo/s3d_protenix/manifest.csv (design_id, chains, kind).
Selection: union of results/denovo/submission_top20_pareto.csv and Pareto front-1 of leaderboard_pareto.csv (gate-passing
population by construction), ranked by pareto_rank, maximin, composite, capped. Controls: trastuzumab Fv (VH+VL+HER2) and
three framework-only decoys. Chain 1 = VHH, last chain = HER2 domain IV (119 aa, first a3m record).
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dn_common as dn  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--cap", type=int, default=40)
a = ap.parse_args()

her2 = [l.strip() for l in (dn.TARGETS / "her2_domIV.a3m").read_text().splitlines()[:2]][1]
assert len(her2) == 119, len(her2)
P = pd.read_csv(dn.OUT / "leaderboard_pareto.csv")
S = pd.read_csv(dn.OUT / "submission_top20_pareto.csv")
sel = P[(P.pareto_rank == 1) | P.design_id.isin(S.design_id)].copy()
sel = sel.sort_values(["pareto_rank", "maximin", "composite"], ascending=[True, False, False]).head(a.cap)
C = pd.read_csv(dn.OUT / "designs_control.csv").set_index("design_id")


def chain(seq):
    return {"proteinChain": {"sequence": seq, "count": 1}}


jobs, man = [], []
for r in sel.itertuples():
    jobs.append({"name": r.design_id, "sequences": [chain(r.vhh_seq), chain(her2)]})
    man.append(dict(design_id=r.design_id, kind="design", is_control=0, chains="VHH,HER2"))
tz = C.loc["ctl-tz-fv"]
jobs.append({"name": "ctl-tz-fv", "sequences": [chain(tz.vhh_seq), chain(json.loads(tz.params_json)["vl_seq"]), chain(her2)]})
man.append(dict(design_id="ctl-tz-fv", kind="fv", is_control=1, chains="VH,VL,HER2"))
for d in ("ctl-fwrand-01", "ctl-fwrand-02", "ctl-fwrand-03"):
    jobs.append({"name": d, "sequences": [chain(C.loc[d].vhh_seq), chain(her2)]})
    man.append(dict(design_id=d, kind="decoy", is_control=1, chains="VHH,HER2"))
Path("nebius-runs").mkdir(exist_ok=True)
Path(dn.ROOT / "nebius-runs" / "ptxtop-input.json").write_text(json.dumps(jobs))
dn.write_csv(Path(__file__).parent / "manifest.csv", man)
print(len(jobs), "jobs;", sum(1 for m in man if m["kind"] == "design"), "designs +", sum(m["is_control"] for m in man), "controls")
