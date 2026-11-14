"""Build the refolding set (top shortlist designs + controls) and the input files for each service.

Chains in every complex: A = trastuzumab VH carrying the CDRH3 under test, B = trastuzumab VL,
C = HER2 domain IV (residues 489-607, UniProt P04626).

Writes
  baseline/outputs/refold_set.csv           one row per complex
  baseline/refold_inputs/jf/<id>.yaml       JapanFold / Boltz-style YAML, one per complex
  nebius-runs/refold/<id>.json              Boltz-2 args for the Nebius platform (needs HER2 a3m; see --her2-a3m)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--top", type=int, default=12)
ap.add_argument("--her2-a3m", default=str(HERE / "refold_inputs" / "her2_domIV.a3m"),
                help="a3m for HER2 domain IV (default: the PDB70 alignment from Nebius msa-search-pdb70, 43 sequences)")
a = ap.parse_args()

HERE = Path(__file__).resolve().parent
VH_T = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
        "SR{H3}YWGQGTLVTVSS")
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
HER2 = (HERE / "refold_inputs" / "her2_domIV.txt").read_text().strip()  # UniProt P04626 residues 489-607
assert len(HER2) == 119 and len(VL) == 107 and len(VH_T.format(H3=C.WT)) == 120

wt = C.encode_idx([C.WT])[0]
rows = []
cand = pd.read_csv(C.OUT / C.vname("blosum_candidates.csv")).head(a.top)
for i, r in enumerate(cand.itertuples(), 1):
    rows.append((f"cand{i:02d}", "candidate", r.cdrh3, f"BLOSUM rank {i}, changed {r.changed_positions}, gate {r.p_mean:.3f}"))
rows.append(("ctrl_wt", "positive control", C.WT, "trastuzumab"))
rows.append(("ctrl_D102E", "positive control", "WGGEGFYAMD", "D102E single, 0.117 nM in SPR (about trastuzumab)"))

# two measured non-binders from the paper's SPR set that also change HER2-contact positions
sub = pd.read_csv(C.DATA / "paper" / "spr_submitted.csv")
spr = C.load_spr()
m = sub.merge(spr.drop_duplicates("cdrh3"), left_on="H3", right_on="cdrh3", how="left")
m = m[m.H3.str.len().eq(10) & m.binding_result.eq("no_reportable_kd") & m.Label.isin(["blosum", "ablang_one", "esm_one"])]
idx = C.encode_idx(m.H3)
chg = idx != wt[None, :]
k = chg.sum(1)
pick = m[(k >= 4) & (k <= 6) & chg[:, [0, 4, 5, 6]].any(1)].sort_values("H3").drop_duplicates("H3")
for j, r in enumerate(pick.head(2).itertuples(), 1):
    rows.append((f"ctrl_nonbinder{j}", "negative control", r.H3, f"measured non-binder in SPR ({r.Label}), no reportable KD"))

df = pd.DataFrame(rows, columns=["id", "role", "cdrh3", "note"])
df["vh"] = [VH_T.format(H3=s) for s in df.cdrh3]
df["n_changed_vs_wt"] = C.hamming_to_wt(C.encode_idx(df.cdrh3))
C.OUT.mkdir(exist_ok=True)
df.to_csv(C.OUT / C.vname("refold_set.csv"), index=False)

jf = C.JF_IN
jf.mkdir(parents=True, exist_ok=True)
neb = C.NEB_IN
neb.mkdir(parents=True, exist_ok=True)
her2_a3m = Path(a.her2_a3m).read_text() if a.her2_a3m else f">query\n{HER2}\n"
for r in df.itertuples():
    (jf / f"{r.id}.yaml").write_text(
        "sequences:\n"
        f"  - protein: {{id: A, sequence: {r.vh}}}\n"
        f"  - protein: {{id: B, sequence: {VL}}}\n"
        f"  - protein: {{id: C, sequence: {HER2}}}\n")

    def poly(i, seq, a3m):
        return {"id": i, "molecule_type": "protein", "sequence": seq,
                "msa": {"msa_search": {"a3m": {"alignment": a3m}}}}

    args = {"polymers": [poly("A", r.vh, f">query\n{r.vh}\n"), poly("B", VL, f">query\n{VL}\n"), poly("C", HER2, her2_a3m)]}
    (neb / f"{r.id}.json").write_text(json.dumps(args))
print(df[["id", "role", "cdrh3", "n_changed_vs_wt", "note"]].to_string(index=False))
print(f"\nHER2 domain IV: {len(HER2)} aa; VH 120 aa; VL 107 aa; complex 346 residues")
print("wrote", C.OUT / "refold_set.csv", "|", jf, "|", neb)
