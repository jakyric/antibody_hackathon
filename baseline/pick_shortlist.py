"""Pick and rank the final designs from the refolding evidence.

Rules (all printed, nothing hidden):
  1. Vetted = OpenFold3 recovers the crystal binding mode in EVERY seed run (the only service that separated
     the measured binders from the measured non-binders).
  2. Rank vetted designs by: number of services recovering the pose (OpenFold3, Boltz-2, OpenDDE) > developability
     flag (fewer aromatic residues first) > OpenFold3 interface PAE on the crystal contacts (lower first)
     > predicted paratope contacts kept > BLOSUM rank.
"""
import numpy as np
import pandas as pd

import common as C

EPITOPE_MAX, CONTACT_MIN = 4.0, 0.5
m = pd.read_csv(C.OUT / C.vname("refold_metrics.csv"))
m["recovers"] = (m.epitope_ca_rmsd_fvframe <= EPITOPE_MAX) & (m.contact_recall >= CONTACT_MIN)
base = pd.read_csv(C.OUT / C.vname("blosum_candidates.csv"))
base["blosum_rank"] = range(1, len(base) + 1)
rset = pd.read_csv(C.OUT / C.vname("refold_set.csv"))

# ---- control check across all OpenFold3 seed runs
of3 = m[m.service.str.startswith("jf_of3_msa")]
print("OpenFold3 control check across seeds (pose recovered / runs):")
for role in ("positive control", "negative control"):
    for id_, g in of3[of3.role == role].groupby("id"):
        print(f"  {role:17s} {id_:16s} {int(g.recovers.sum())}/{len(g)}")

# ---- per candidate evidence
cand = of3[of3.role == "candidate"].groupby("id").agg(
    of3_seeds=("recovers", "size"), of3_pass_frac=("recovers", "mean"),
    of3_ipae=("ipae_crystal_contacts", "median"), of3_iptm=("iptm_ab_her2", "median"))
prim = m[m.service.isin(["jf_abag_msa", "jf_of3_msa", "nebius_boltz2"]) & (m.role == "candidate")]
cand["models_recovering"] = prim.groupby("id").recovers.sum()
cand["which"] = prim[prim.recovers].groupby("id").service.apply(
    lambda s: ",".join({"jf_abag_msa": "OpenDDE", "jf_of3_msa": "OpenFold3", "nebius_boltz2": "Boltz-2"}[x] for x in sorted(s)))
cand["paratope_kept"] = prim.groupby("id")[["W99_kept", "G103_kept", "F104_kept", "Y105_kept"]].mean().mean(axis=1)
cand = cand.join(rset.set_index("id")[["cdrh3"]]).reset_index()
cand["cand_no"] = cand.id.str[-2:].astype(int)
cand = cand.merge(base[["cdrh3", "blosum_rank", "changed_positions", "p_mean", "liabilities"]], on="cdrh3", how="left")


def dev(s):
    arom = sum(s.count(x) for x in "WFY")
    charge = sum(s.count(x) for x in "KR") - sum(s.count(x) for x in "DE")
    return arom, charge


cand["aromatics"], cand["net_charge"] = zip(*[dev(s) for s in cand.cdrh3])
cand["keeps_M107"] = cand.cdrh3.str[8] == "M"
cand["keeps_DG"] = cand.cdrh3.str[3:5] == "DG"
cand["arom_flag"] = cand.aromatics >= 6

vet = cand[cand.of3_pass_frac == 1.0].copy()
vet = vet.sort_values(["models_recovering", "arom_flag", "of3_ipae", "paratope_kept", "blosum_rank"],
                      ascending=[False, True, True, False, True]).reset_index(drop=True)
vet.insert(0, "final_rank", range(1, len(vet) + 1))
vet["tier"] = np.where(vet.models_recovering == 3, "A: all three services", "B: OpenFold3 + Boltz-2")
cols = ["final_rank", "id", "cdrh3", "tier", "of3_seeds", "of3_ipae", "of3_iptm", "paratope_kept", "aromatics",
        "net_charge", "keeps_M107", "keeps_DG", "liabilities", "p_mean", "blosum_rank", "changed_positions", "which"]
vet[cols].round(2).to_csv(C.OUT / C.vname("final_shortlist.csv"), index=False)
print("\nVetted designs (OpenFold3 recovers the pose in every seed):")
print(vet[cols].round(2).to_string(index=False))
out = cand[~cand.id.isin(vet.id)].sort_values("cand_no")
print("\nNot selected:")
print(out[["id", "cdrh3", "of3_seeds", "of3_pass_frac", "models_recovering", "which"]].round(2).to_string(index=False))
print("\nsaved", C.OUT / C.vname("final_shortlist.csv"))
