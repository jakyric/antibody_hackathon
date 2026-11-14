"""Stage 4: fold every validated sequence (plus four controls) with JapanFold and download the structures.

The model is OpenDDE antibody-antigen: of the three models tried it was the only one that put both measured binders
above both measured non-binders. HER2 domain IV is chain C. PAE matrices are written (needed for ipSAE / iLIS) but
only downloaded later, for designs that survive the contact filter, because they are ~1.4 MB each.
"""
from __future__ import annotations

import pandas as pd

from .common import CONTROLS, build_fv, fv_yaml, out, save_summary
from .japanfold import JapanFold


def run(cfg: dict) -> None:
    f = cfg["fold"]
    jf = JapanFold(f["base_url"])
    valid = pd.read_csv(out(cfg, "panel_valid.csv"))
    targets = [{"name": r.sequence_id, "content": fv_yaml(r.heavy_sequence, r.light_sequence)} for r in valid.itertuples()]
    for name, h3 in CONTROLS.items():
        targets.append({"name": name, "content": fv_yaml(*build_fv(h3))})
    job = f.get("job_id")
    if not job:
        job = jf.submit(f["model"], "her2-cdrh3-vernier", targets, {"output_format": "pdb", "write_pae": bool(f["write_pae"])})
        print(f"  submitted {len(targets)} targets as job {job}")
    state = jf.wait(job)
    if state["status"] != "succeeded":
        raise SystemExit(f"fold job ended with status {state['status']}: {state.get('error')}")
    res = jf.results(job)
    rows = [{"id": x["id"], "iptm": x["iptm"], "complex_plddt": x.get("complex_plddt"), "confidence_score": x.get("confidence_score")}
            for x in res["rows"]]
    pd.DataFrame(rows).to_csv(out(cfg, "fold_scores.csv"), index=False)
    out(cfg, "job_id.txt").write_text(job)
    jf.fetch_many(job, [f"{r['id']}.pdb" for r in rows], out(cfg, "structures"))
    print(f"  {len(rows)} structures in {out(cfg, 'structures')}")
    save_summary(cfg, n_folded=len(rows), fold_job=job, n_docked=int(sum(r["iptm"] >= cfg["rank"]["min_iptm"] for r in rows if not r["id"].startswith("ctrl_"))))
