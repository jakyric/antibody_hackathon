"""S6, the last step: Pareto scoring on top of the leaderboard.

    denovo/.venv/bin/python denovo/pareto.py                 # writes leaderboard_pareto.csv, submission_top20_pareto.csv
    denovo/.venv/bin/python denovo/pareto.py --objectives bind,comp,therm --top 20

Why Pareto, not just the weighted composite: the organisers tier a design by HOW MANY of the three judged properties
(affinity, epitope preservation, thermostability) beat trastuzumab, and S-tier needs all of them. A weighted sum lets a
very strong binder with a poor Tm outrank a balanced design; the weakest objective is what loses the tier. So:

  1. Population = designs passing the hard gate (`gate_strict`: abguard de novo + 2-point margin + canonical Cys etc.;
     falls back to `gate_vm` if fewer than 20 pass). Controls are never ranked.
  2. Objectives, all maximised, each a within-run percentile (monotone, so dominance is unchanged by the transform):
       bind  = OpenFold3 consensus blended with the orthogonal models (OpenDDE-abag, ESMFold-2)   -> affinity proxy
       comp  = epitope overlap with trastuzumab + CDR/Fv contact fraction                          -> epitope proxy
       therm = NanoMelt predicted Tm                                                              -> thermostability
       yield = ESMC-Nb monoyield (falls back to the nine-rule score if monoyield is not yet available)
       dev   = TNP flags and sequence liabilities
     An objective with data for fewer than 30% of the population is dropped (printed), so the step works from the first
     minute and sharpens as stages fill in. A design missing one used objective gets 0 there (pessimistic): an unmeasured
     design can never dominate a measured one, and `pareto_complete` says whether every used objective was measured.
  3. Non-dominated sorting (NSGA-II style): `pareto_rank` 1 = the Pareto front, 2 = the front after removing it, ...
     `crowding` = distance to neighbours on the same front (spread), `maximin` = the design's WORST objective percentile,
     `weakest` = which objective that is.
  4. Submission list: fronts in order; within a front by maximin (balanced designs first, aligned with S-tier), then the
     weighted composite; skip near-duplicate CDR3s (edit distance < 3) and keep framework x hotspot diversity (<= 50% of
     the list from one group); prefer `pareto_complete` designs, fill with incomplete ones only if needed (flagged).

Limits: percentiles are within this run, not against trastuzumab; proxies are filters, not measurements of affinity.
"""
from __future__ import annotations

import argparse
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dn_common as dn  # noqa: E402
import rank as R  # noqa: E402

OBJ = {"bind": "t_bind", "comp": "t_comp", "therm": "t_therm", "yield": "t_yield", "dev": "t_dev"}


def dominated_by(P: np.ndarray) -> np.ndarray:
    n = len(P)
    out = np.zeros(n, int)
    for i in range(n):
        out[i] = int(((P >= P[i]).all(1) & (P > P[i]).any(1)).sum())
    return out


def fronts(P: np.ndarray) -> np.ndarray:
    rank = np.zeros(len(P), int)
    remaining = np.arange(len(P))
    k = 1
    while len(remaining):
        sub = P[remaining]
        nd = dominated_by(sub) == 0
        rank[remaining[nd]] = k
        remaining = remaining[~nd]
        k += 1
    return rank


def crowding(P: np.ndarray, rank: np.ndarray) -> np.ndarray:
    d = np.zeros(len(P))
    for k in np.unique(rank):
        idx = np.where(rank == k)[0]
        if len(idx) <= 2:
            d[idx] = np.inf
            continue
        for m in range(P.shape[1]):
            o = idx[np.argsort(P[idx, m])]
            span = P[o[-1], m] - P[o[0], m]
            d[o[0]] = d[o[-1]] = np.inf
            if span > 0:
                d[o[1:-1]] += (P[o[2:], m] - P[o[:-2], m]) / span
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--objectives", default="bind,dev,yield,therm")  # iPTM/ipSAE per scorer, TNP, yield + nine rules, NanoMelt Tm (user 15:44)
    ap.add_argument("--require", default="bind,dev", help="objectives that must be MEASURED for a design to enter the front")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--min-coverage", type=float, default=0.30)
    a = ap.parse_args()
    ns = types.SimpleNamespace(w_bind=0.50, w_comp=0.0, w_therm=0.20, w_yield=0.15, w_dev=0.15, w_rules=0.0)
    D = R.build(ns)
    D = D[D.is_control == 0].copy()
    pop = D[D.gate_strict == 1]
    gate = "gate_strict"
    if len(pop) < 20:
        pop, gate = D[D.gate_vm == 1], "gate_vm"
    pop = pop.copy()
    if pop.empty:
        sys.exit("no designs pass the gate yet")
    # TNP is run only on the most promising designs, so instead of dropping the TNP objective for low coverage, restrict the
    # population to designs where the required objectives are really measured (binding scores AND TNP).
    for o in [x for x in a.require.split(",") if x]:
        c = OBJ[o]
        if pop[c].notna().sum() >= 10:
            n0 = len(pop)
            pop = pop[pop[c].notna()].copy()
            print(f"required objective '{o}': population {n0} -> {len(pop)} (measured only)")
    # monoyield may not exist yet: nine-rule score is the stand-in yield prior
    if pop["t_yield"].notna().mean() < a.min_coverage and pop["t_rules"].notna().mean() >= a.min_coverage:
        pop["t_yield"] = pop["t_rules"]
        print("note: monoyield coverage low, using the nine-rule score as the yield objective")
    used, dropped = [], []
    for o in a.objectives.split(","):
        c = OBJ[o]
        (used if pop[c].notna().mean() >= a.min_coverage else dropped).append(o)
    if dropped:
        print("objectives dropped (coverage < %d%%): %s" % (a.min_coverage * 100, ", ".join(dropped)))
    if not used:
        sys.exit("no objective has enough data yet")
    cols = [OBJ[o] for o in used]
    M = pop[cols].to_numpy(float)
    P = np.nan_to_num(M, nan=0.0)
    pop["pareto_rank"] = fronts(P)
    pop["n_dominating"] = dominated_by(P)
    pop["crowding"] = crowding(P, pop["pareto_rank"].to_numpy())
    pop["maximin"] = P.min(1)
    pop["weakest"] = [used[i] for i in P.argmin(1)]
    pop["pareto_objectives"] = ",".join(used)
    pop["pareto_complete"] = (~np.isnan(M)).all(1).astype(int)
    pop["pareto_reason"] = [
        f"Pareto front {int(r.pareto_rank)} over {','.join(used)}; weakest objective {r.weakest} (pctl {r.maximin:.2f})"
        + ("" if r.pareto_complete else "; some objectives not measured yet") for r in pop.itertuples()]
    keep = ["pareto_rank", "design_id", "generator", "framework", "hotspot_set", "maximin", "weakest", "crowding",
            "n_dominating", "pareto_complete", "n_flags", "flags", "missing_evidence", "composite", "pareto_reason"] + cols + ["rank_reason", "cdr3", "vhh_seq"]
    out = pop.sort_values(["pareto_rank", "maximin", "composite"], ascending=[True, False, False])
    out[[c for c in keep if c in out] + [c for c in out.columns if c not in keep]].to_csv(
        dn.OUT / "leaderboard_pareto.csv", index=False)

    # diversity-aware selection
    # binding first: within a front the best iPTM/ipSAE designs come first (user: 'definitely has the best iPTM')
    # order: eligible (gate) -> complete evidence -> Pareto tier -> fewest review flags -> best binding
    cand = out.sort_values(["pareto_complete", "pareto_rank", "n_flags", "t_bind", "maximin"], ascending=[False, True, True, False, False])
    chosen, cdrs, count = [], [], {}
    cap = max(1, int(np.ceil(0.5 * a.top)))
    for _, r in cand.iterrows():
        c3 = str(r.get("cdr3", ""))
        grp = f"{r.get('framework', '')}|{r.get('hotspot_set', '')}"
        if any(R.edit(c3, o) < 3 for o in cdrs) or count.get(grp, 0) >= cap:
            continue
        chosen.append(r)
        cdrs.append(c3)
        count[grp] = count.get(grp, 0) + 1
        if len(chosen) == a.top:
            break
    S = pd.DataFrame(chosen)
    if not S.empty:
        S.insert(0, "slot", np.arange(1, len(S) + 1))
        S[["slot"] + [c for c in keep if c in S]].to_csv(dn.OUT / "submission_top20_pareto.csv", index=False)
    sizes = out.pareto_rank.value_counts().sort_index().head(5).to_dict()
    print(f"pareto: population {len(out)} ({gate}), objectives {used}, front sizes {sizes}, "
          f"complete {int(out.pareto_complete.sum())}, submission rows {len(S)}")
    if not S.empty:
        print(S[["slot", "design_id", "pareto_rank", "maximin", "weakest", "pareto_complete"]].head(10).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
