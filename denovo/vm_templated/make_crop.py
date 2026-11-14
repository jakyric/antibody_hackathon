"""Cropped-HER2 re-run (user, 16:4x): drop the first ~50 residues of the 119-aa domain IV (HER2 489-538), keep 539-607, so the models
cannot place the VHH on the far side of the domain. 539 is the first residue that severs no disulfide (Cys538-Cys522 are both removed) and
every trastuzumab epitope residue (557-605) is kept. Builds, from results/denovo/leaderboard.csv:
  results/denovo/crop/set.csv          design_id, priority, why     (gate_strict first, then the best-evidence eligible designs, then controls)
  results/denovo/crop/template_her2_crop.cif   template with chains C (539-580) and D (591-607) only, residue numbers unchanged
  results/denovo/crop/chunks/q_NN.json  OpenFold3 queries of 20 designs each, in priority order
Usage: denovo/.venv/bin/python denovo/vm_templated/make_crop.py [--n 150]
"""
import argparse
import csv
import json
import os

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=150, help="designs (besides controls) to re-run")
a = ap.parse_args()
OUT = "results/denovo/crop"
os.makedirs(f"{OUT}/chunks", exist_ok=True)
CUT = 539                                      # first kept HER2 residue
FULL = open("denovo/targets/her2_domIV.a3m").read().split("\n")[1].strip()   # 489-607
assert len(FULL) == 119
C_FULL, D_FULL = FULL[:92], FULL[102:]         # chain C = 489-580, chain D = 591-607 (581-590 are unresolved in 1N8Z)
assert len(D_FULL) == 17, len(D_FULL)
C = C_FULL[CUT - 489:]                         # 539-580
HER2C = FULL[CUT - 489:]                       # 539-607 incl. the 581-590 loop, contiguous (for OpenDDE)
print(f"crop: chain C {CUT}-580 ({len(C)} aa), chain D 591-607 ({len(D_FULL)} aa); contiguous crop for OpenDDE {len(HER2C)} aa")

# ---- template
lines = open("denovo/vm_templated/template_her2.cif").read().splitlines()
hdr = [l.strip().split(".")[1] for l in lines if l.startswith("_atom_site.")]
ix = {h: i for i, h in enumerate(hdr)}
out = []
for l in lines:
    p = l.split()
    if l.startswith("ATOM"):
        if p[ix["label_asym_id"]] == "C":
            k = int(p[ix["label_seq_id"]])
            if k < CUT - 488:
                continue
            p[ix["label_seq_id"]] = str(k - (CUT - 489))
            l = " ".join(p)
    elif len(p) == 4 and p[0] == "3" and p[3] == "n" and p[1].isdigit():   # entity_poly_seq rows of chain C
        k = int(p[1])
        if k < CUT - 488:
            continue
        l = f"3 {k - (CUT - 489)} {p[2]} n"
    out.append(l)
open(f"{OUT}/template_her2_crop.cif", "w").write("\n".join(out) + "\n")
print("template atoms:", sum(1 for l in out if l.startswith("ATOM")), "(was", sum(1 for l in lines if l.startswith("ATOM")), ")")

# ---- design set
L = pd.read_csv("results/denovo/leaderboard.csv")
D = {r["design_id"]: r for r in csv.DictReader(open("results/denovo/designs_rfd4jf.csv"))}
for f in ("designs_control.csv", "designs_control_krisha.csv"):
    for r in csv.DictReader(open(f"results/denovo/{f}")):
        D[r["design_id"]] = r
real = L[L.is_control == 0]
strict = real[real.gate_strict == 1].sort_values("t_bind", ascending=False)
rest = real[(real.gate_vm == 1) & (real.gate_strict != 1)].sort_values("t_bind", ascending=False, na_position="last")
chosen = list(strict.design_id) + list(rest.design_id)[: max(0, a.n - len(strict))]
ctl = [c for c in L[L.is_control == 1].design_id if c.startswith(("ctl-kri-0", "ctl-bcii10", "ctl-fwrand"))][:20] + ["ctl-kri-01"]
ctl = list(dict.fromkeys(c for c in ctl if c in D and c != "ctl-tz-fv"))
rows = [(i, k, "gate_strict" if i in set(strict.design_id) else "eligible_top_binding") for k, i in enumerate(chosen)] + [(i, len(chosen) + k, "control") for k, i in enumerate(ctl)]
pd.DataFrame(rows, columns=["design_id", "priority", "why"]).to_csv(f"{OUT}/set.csv", index=False)
print(f"set: {len(chosen)} designs ({len(strict)} strict) + {len(ctl)} controls")

# ---- OpenFold3 chunks, 20 designs per chunk
ids = [r[0] for r in rows]
rid = lambda a, b: [str(x) for x in range(a, b + 1)]  # noqa: E731
for n, k in enumerate(range(0, len(ids), 20)):
    qs = []
    for i in ids[k:k + 20]:
        qs.append({"query_id": i, "chains": [
            {"chain_id": "A", "sequence": D[i]["vhh_seq"], "msa_path": None},
            {"chain_id": "C", "sequence": C, "residue_ids": rid(CUT, 580), "msa_path": None},
            {"chain_id": "D", "sequence": D_FULL, "residue_ids": rid(591, 607), "msa_path": None}]})
    json.dump({"queries": qs}, open(f"{OUT}/chunks/q_{n:02d}.json", "w"))
print(f"{n + 1} chunks")
open(f"{OUT}/her2_crop.txt", "w").write(HER2C)
