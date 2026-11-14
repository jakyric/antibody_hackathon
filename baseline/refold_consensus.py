"""Combine the refolding services into one verdict per design.

A complex "recovers the crystal binding mode" if, after aligning the Fv, the HER2 epitope residues sit within
EPITOPE_MAX angstrom of the 1N8Z crystal AND at least CONTACT_MIN of the crystal's 42 contact pairs are reproduced.
Each service is first judged on the four controls (2 measured binders, 2 measured non-binders); services that do
not separate them are shown but not counted towards the consensus.
"""
import pandas as pd

import common as C

EPITOPE_MAX = 4.0
CONTACT_MIN = 0.5

m = pd.read_csv(C.OUT / C.vname("refold_metrics.csv"))
m = m[~m.service.str.contains(r"_s\d+$")]  # primary (seed 0 / default) runs only; seed runs are checked separately
m["recovers"] = (m.epitope_ca_rmsd_fvframe <= EPITOPE_MAX) & (m.contact_recall >= CONTACT_MIN)

print("Control check (does the service separate measured binders from measured non-binders?)")
trusted = []
for svc, g in m.groupby("service"):
    pos = g[g.role == "positive control"]
    neg = g[g.role == "negative control"]
    ok_pos, ok_neg = int(pos.recovers.sum()), int((~neg.recovers).sum())
    good = ok_pos == len(pos) and ok_neg == len(neg) and len(pos) > 0
    if good:
        trusted.append(svc)
    print(f"  {svc:16s} binders recovered {ok_pos}/{len(pos)} | non-binders rejected {ok_neg}/{len(neg)} "
          f"-> {'TRUSTED' if good else 'not counted'}")

cand = m[m.role == "candidate"]
t = cand.pivot(index="id", columns="service", values="recovers")
cols = list(t.columns)
t["n_trusted_pass"] = t[[c for c in cols if c in trusted]].sum(axis=1) if trusted else 0
t["n_all_pass"] = t[cols].sum(axis=1)
meta = pd.read_csv(C.OUT / C.vname("refold_set.csv")).set_index("id")[["cdrh3"]]
paratope = cand.groupby("id")[["W99_kept", "G103_kept", "F104_kept", "Y105_kept"]].mean().mean(axis=1).rename("paratope_kept_mean")
out = meta.join(t, how="inner").join(paratope).sort_values(["n_trusted_pass", "n_all_pass", "paratope_kept_mean"], ascending=False)
print(f"\nCandidates (pose recovered = epitope within {EPITOPE_MAX} A of crystal and >= {CONTACT_MIN:.0%} of crystal contacts):")
print(out.round(2).to_string())
out.to_csv(C.OUT / C.vname("refold_consensus.csv"))
print("\nsaved", C.OUT / C.vname("refold_consensus.csv"))
