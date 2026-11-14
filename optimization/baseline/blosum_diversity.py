"""Binding-preserving CDRH3 diversification: data-chosen positions + BLOSUM-45 substitutions.

Pipeline
  1. Position evidence: HER2 contact distance (1N8Z), FACS and SPR "does changing position p hurt?" models.
     Contact positions are excluded; the best 5 of the rest (by the worse of the two data estimates) are used.
  2. Substitutions ranked with the paper's own BLOSUM recipe (Chinery et al.: BLOSUM-45, lambda 0.25,
     SAbDab CDRH3 background frequencies), restricted to real substitutions.
  3. Candidates change 4 or 5 of those positions (4 = the minimum for CDRH3 identity < 70% under either the
     IMGT 13-residue or Kabat 11-residue CDR3 that immunum reports).
  4. Filters: not within 2 mutations of any FACS/SPR sequence (so already-tested designs are out), no new
     liability motifs, gate score from the classifiers, and an "abguard-lite" check with immunum.
  5. Greedy diverse shortlist ranked by BLOSUM joint probability.

Run: optimization/baseline/.venv/bin/python optimization/baseline/blosum_diversity.py [--m 6] [--top 24] [--threshold 0.9]
"""
import argparse
import itertools
import math
import re

import blosum as bl
import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--m", type=int, default=6, help="top substitutions per position")
ap.add_argument("--top", type=int, default=24)
ap.add_argument("--threshold", type=float, default=0.9, help="min gate score p_mean")
ap.add_argument("--positions", type=int, default=5, help="number of positions to allow changes at")
ap.add_argument("--allow-contacts", action="store_true",
                help="do not exclude HER2-contact positions (1, 5, 6, 7) from the candidate positions")
ap.add_argument("--allow-seen", action="store_true",
                help="keep substitutions already measured as single mutants in the paper's SPR set (default: drop them)")
a = ap.parse_args()

wt = C.encode_idx([C.WT])[0]
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- 1. position evidence
facs = C.load_facs()
chg = (C.encode_idx(facs.cdrh3) != wt[None, :]).astype(float)
facs_coef = LogisticRegression(C=100, max_iter=500).fit(chg, facs.label).coef_[0]

sub = pd.read_csv(C.DATA / "paper" / "spr_submitted.csv")
spr = C.load_spr()

# single substitutions already measured in the paper's SPR set ("checked in the paper")
seen = set()
for s_ in spr.cdrh3[spr.cdrh3.str.len() == C.L]:
    d = [(i, c) for i, (c, w) in enumerate(zip(s_, C.WT)) if c != w]
    if len(d) == 1:
        seen.add((d[0][0] + 1, d[0][1]))
m = sub.merge(spr.drop_duplicates("cdrh3"), left_on="H3", right_on="cdrh3", how="left")
m = m[m.Label.isin(["blosum", "ablang_all", "ablang_one", "esm_one", "protein_mpnn"]) & m.H3.str.len().eq(10)]
m = m[m.binding_result.isin(["quantified_binding", "no_reportable_kd"])]
X = (C.encode_idx(m.H3) != wt[None, :]).astype(float)
y = (m.binding_result == "quantified_binding").astype(int).values
boots = []
for _ in range(200):
    i = rng.integers(0, len(y), len(y))
    boots.append(LogisticRegression(C=1.0, max_iter=300).fit(X[i], y[i]).coef_[0])
spr_coef = np.mean(boots, 0)

# HER2 contact distance from the 1N8Z Fv-HER2 structure used by the Apheris HER2 example
hdr, rows = [], []
for line in open(C.DATA / "her2" / "her2_fv_1n8z.cif"):
    if line.startswith("_atom_site."):
        hdr.append(line.strip().split(".")[1])
    elif line.startswith(("ATOM", "HETATM")):
        rows.append(line.split())
cif = pd.DataFrame(rows, columns=hdr[: len(rows[0])])
xyz = cif[["Cartn_x", "Cartn_y", "Cartn_z"]].astype(float).to_numpy()
resn = cif.label_seq_id.astype(int).to_numpy()
asym = cif.label_asym_id.to_numpy()
her2 = xyz[np.isin(asym, ["C", "D"])]
dist = []
for s in range(99, 109):
    atoms = xyz[(asym == "B") & (resn == s)]
    dist.append(float(np.sqrt(((atoms[:, None] - her2[None]) ** 2).sum(-1)).min()))

ev = pd.DataFrame({"pos": range(1, 11), "wt": list(C.WT), "her2_dist_A": np.round(dist, 1),
                   "facs_logit_if_mutated": facs_coef.round(2), "spr_logit_if_mutated": spr_coef.round(2)})
ev["contact"] = ev.her2_dist_A < 4.5
ev["worst_case"] = np.minimum(ev.facs_logit_if_mutated, ev.spr_logit_if_mutated)
elig = (ev if a.allow_contacts else ev[~ev.contact]).sort_values("worst_case", ascending=False)
chosen = sorted(elig.pos.head(a.positions).tolist())
ev["chosen"] = ev.pos.isin(chosen)
print("Position evidence (more negative = changing it hurts more):")
print(ev.to_string(index=False))
print(f"\nchosen positions (best {a.positions} {'of all 10' if a.allow_contacts else 'non-contact'} by worst-case of FACS and SPR): {chosen}")

# ---------------------------------------------------------------- 2. BLOSUM-45 substitutions (paper recipe)
H3_COUNTS = {"A": 6710, "C": 555, "D": 6442, "E": 1636, "F": 3404, "G": 6582, "H": 1081, "I": 1342, "K": 970,
             "L": 2509, "M": 1218, "N": 1658, "P": 2090, "Q": 767, "R": 5196, "S": 4350, "T": 2835, "V": 2975,
             "W": 1686, "Y": 9028}
tot = sum(H3_COUNTS.values())
bg = {k: v / tot for k, v in H3_COUNTS.items()}
B45 = bl.BLOSUM(45)
LAMBDA = 0.25  # the paper's value (default lowered from ~0.35 to down-weight the original residue)


def subst_probs(a_res):
    """P(b | a, b != a) proportional to f_b * exp(lambda * s_ab): the paper's reverse-BLOSUM, minus 'stay'."""
    w = {b: bg[b] * math.exp(LAMBDA * B45[a_res][b]) for b in C.AA if b != a_res}
    z = sum(w.values())
    return {b: v / z for b, v in w.items()}


subs = {}
print("\nTop BLOSUM-45 substitutions per chosen position (paper recipe):")
for p in chosen:
    q = subst_probs(C.WT[p - 1])
    ranked = sorted(q.items(), key=lambda kv: -kv[1])
    dropped = [b for b, _ in ranked if (p, b) in seen] if not a.allow_seen else []
    top = [(b, v) for b, v in ranked if b not in dropped][: a.m]
    subs[p] = top
    print(f"  pos {p:2d} {C.WT[p - 1]} -> " + "  ".join(f"{b}:{v:.2f}" for b, v in top)
          + (f"   [dropped, already measured as singles: {''.join(dropped)}]" if dropped else ""))

# ---------------------------------------------------------------- 3. enumerate candidates
cands = []
for k in (4, 5):
    if k > len(chosen):
        continue
    for pos_set in itertools.combinations(chosen, k):
        for combo in itertools.product(*[subs[p] for p in pos_set]):
            seq = list(C.WT)
            lp = 0.0
            for p, (b, pr) in zip(pos_set, combo):
                seq[p - 1] = b
                lp += math.log(pr)
            cands.append(("".join(seq), k, lp, "".join(str(p) if p < 10 else "X" for p in pos_set)))
df = pd.DataFrame(cands, columns=["cdrh3", "n_changed", "blosum_logp", "changed_positions"])
df["changed_positions"] = [",".join(str(p) for p in itertools.compress(range(1, 11), [c != w for c, w in zip(s, C.WT)]))
                           for s in df.cdrh3]
print(f"\ncandidates enumerated: {len(df):,}")

# ---------------------------------------------------------------- 4. filters
idx = C.encode_idx(df.cdrh3)
prov = C.ProvidedData(C.provided_sequences())
df["within_1_of_provided"] = prov.has_neighbour_within_1(idx)
df["within_2_of_provided"] = prov.has_neighbour_within_2(idx)

LIAB = ["NG", "NS", "NT", "NN", "DG", "DS", "DD", "DT", "DH", "DN", "DP", "DQ", "KE", "KD", "EK", "ED", "M", "W"]


def liabilities(s):
    n = sum(len(re.findall(f"(?={r})", s)) for r in LIAB)
    n += len(re.findall(r"N[^P][ST]", s)) + sum(s.count(x) for x in ("RGD", "RYD", "LDV", "GPR"))
    return n


wt_l = liabilities(C.WT)
df["liabilities"] = [liabilities(s) for s in df.cdrh3]
df["has_cys"] = df.cdrh3.str.contains("C")
hgb = joblib.load(C.OUT / "hgb_full.joblib")["model"]
cnn = C.AuthorsCNN(C.MEHTA / "saved_weights.h5")
df["p_hgb"] = hgb.predict_proba(idx.astype(np.float32))[:, 1]
df["p_cnn"] = cnn.predict_proba(idx)
df["p_mean"] = (df.p_hgb + df.p_cnn) / 2
print(f"WT CDRH3 liabilities (paper's regex set): {wt_l}")
steps = [("all enumerated", np.ones(len(df), bool)),
         ("not within 2 mutations of any provided sequence", ~df.within_2_of_provided),
         ("... and no more liabilities than trastuzumab", ~df.within_2_of_provided & (df.liabilities <= wt_l)),
         ("... and no cysteine", ~df.within_2_of_provided & (df.liabilities <= wt_l) & ~df.has_cys),
         (f"... and gate score >= {a.threshold}", ~df.within_2_of_provided & (df.liabilities <= wt_l) & ~df.has_cys & (df.p_mean >= a.threshold))]
for name, msk in steps:
    print(f"  {name:55s} {int(msk.sum()):7,d}")
ok = steps[-1][1]
pool = df[ok].sort_values("blosum_logp", ascending=False)

# ---------------------------------------------------------------- 5. diverse shortlist
chosen_rows, chosen_idx = [], []
for _, r in pool.iterrows():
    v = C.encode_idx([r.cdrh3])[0]
    if all((v != u).sum() >= 3 for u in chosen_idx):  # shortlist members differ from each other by >= 3
        chosen_rows.append(r)
        chosen_idx.append(v)
    if len(chosen_rows) >= a.top:
        break
short = pd.DataFrame(chosen_rows).reset_index(drop=True)

# ---------------------------------------------------------------- abguard-lite with immunum
import immunum

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC" \
     "SR{H3}YWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
assert VH.format(H3=C.WT).startswith("EVQLVESGGG") and "WGGDGFYAMDYWGQ" in VH.format(H3=C.WT)
annH = {s: immunum.Annotator(chains=["H"], scheme=s, min_confidence=0.0) for s in ("imgt", "kabat")}
annL = immunum.Annotator(chains=["K", "L"], scheme="imgt", min_confidence=0.0)
wt_seg = {s: annH[s].segment(VH.format(H3=C.WT)).cdr3 for s in annH}


def ident(a_, b_):
    return sum(x == y for x, y in zip(a_, b_)) / max(len(a_), len(b_)) if len(a_) == len(b_) else float("nan")


rows = []
for s in short.cdrh3:
    vh = VH.format(H3=s)
    nh, nl = annH["imgt"].number(vh), annL.number(VL)
    out = {"cdrh3": s, "numbering_conf_H": round(nh.confidence, 2), "numbering_conf_L": round(nl.confidence, 2),
           "fv_length": len(vh) + len(VL), "cys_H": vh.count("C"), "cys_L": VL.count("C")}
    for sc in ("imgt", "kabat"):
        cdr3 = annH[sc].segment(vh).cdr3
        out[f"cdr3_{sc}"] = cdr3
        out[f"cdr3_{sc}_identity"] = round(ident(cdr3, wt_seg[sc]), 3)
    out["abguard_lite_pass"] = (out["numbering_conf_H"] >= 0.1 and out["fv_length"] <= 250 and out["cys_H"] % 2 == 0
                                and out["cdr3_imgt_identity"] < 0.70 and out["cdr3_kabat_identity"] < 0.70)
    rows.append(out)
guard = pd.DataFrame(rows)
short = short.merge(guard, on="cdrh3")

# single substitutions of each candidate that were already measured in the paper's SPR set
short["singles_seen_in_spr"] = [
    " ".join(f"{C.WT[p - 1]}{p}{s_[p - 1]}" for p in [int(x) for x in cp.split(",")] if (p, s_[p - 1]) in seen) or "-"
    for s_, cp in zip(short.cdrh3, short.changed_positions)]

C.OUT.mkdir(exist_ok=True)
cols = ["cdrh3", "changed_positions", "n_changed", "blosum_logp", "p_mean", "p_hgb", "p_cnn", "liabilities",
        "within_1_of_provided", "within_2_of_provided", "cdr3_imgt", "cdr3_imgt_identity", "cdr3_kabat_identity",
        "numbering_conf_H", "abguard_lite_pass", "singles_seen_in_spr"]
short[cols].round(4).to_csv(C.OUT / C.vname("blosum_candidates.csv"), index=False)
df.sort_values("blosum_logp", ascending=False).round(4).to_csv(C.OUT / C.vname("blosum_all_candidates.csv"), index=False)
ev.to_csv(C.OUT / C.vname("position_evidence.csv"), index=False)

print(f"\nshortlist of {len(short)} (greedy, ranked by BLOSUM joint probability, members >= 3 mutations apart):")
show = short[["cdrh3", "changed_positions", "blosum_logp", "p_mean", "liabilities", "abguard_lite_pass", "singles_seen_in_spr"]].round(3)
print(show.to_string(index=False))
print(f"\nimmunum: trastuzumab CDR3 = {wt_seg['imgt']} (IMGT, {len(wt_seg['imgt'])} aa) / {wt_seg['kabat']} (Kabat, {len(wt_seg['kabat'])} aa)")
print("saved baseline/outputs/blosum_candidates.csv, blosum_all_candidates.csv, position_evidence.csv")
