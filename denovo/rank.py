"""S5: merge every stage into one leaderboard + a diversity-aware top-20 submission list.

    denovo/.venv/bin/python denovo/rank.py            # writes results/denovo/leaderboard.csv, submission_top20.csv
    denovo/.venv/bin/python denovo/rank.py --top 20 --w-bind 0.4 ...

Reads (all optional except designs): designs_*.csv, rules.csv, score_of3.csv, score_orth.csv, props.csv,
designs_boltzgen_metrics.csv. Missing columns/files are skipped and the weights renormalised, so it is usable from the
first minute and gets sharper as stages fill in.

Composite (all terms are within-run percentile ranks in 0..1, controls excluded from the reference population):
  bind  : S3 `of3_consensus_rank_score` (else mean pct of of3_ipsae, of3_iptm, -of3_min_ipae), blended 60/40 with
          the orthogonal `consensus_rank_score` (OpenDDE / ESMFold-2 / ...) when that exists
  comp  : epitope/competition proxy = mean pct of of3_epi_overlap and of3_cdr_fv_frac4
  therm : nm_tm        yield : my_yield        dev : 1 - 0.5*tnp_n_red - 0.15*n_liab (floored at 0)
Weights default 0.40 / 0.15 / 0.20 / 0.15 / 0.10. These are priors, not fitted: the point is a transparent, tunable
filter-then-rank, and every row carries a plain-language `rank_reason`.

Gates (de novo track): `gate_strict` = abguard (VM) pass AND local identity to nearest known < 0.70 AND hard rules;
`gate_vm` = what the VM abguard alone says. Designs that fail only the identity gate are flagged `needs_fwdiv=1`
(framework diversification candidates). If nothing passes strict yet, the submission list falls back to the
best-by-composite designs and says so in `gate_note`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dn_common as dn  # noqa: E402

OUT = dn.OUT


def rd(name: str) -> pd.DataFrame | None:
    p = OUT / name
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p)
    except Exception as e:  # half-written file or bad row: skip this cycle
        print(f"warn: cannot read {name}: {e}", file=sys.stderr)
        return None


def pct(s: pd.Series, ref: pd.Series, higher_better: bool = True) -> pd.Series:
    """Percentile of each value of s inside the reference population ref (NaN stays NaN)."""
    r = np.sort(ref.dropna().to_numpy(float))
    if len(r) < 2:
        return pd.Series(np.where(s.notna(), 0.5, np.nan), index=s.index)
    v = s.to_numpy(float)
    p = np.searchsorted(r, v, side="right") / len(r)
    p = np.where(np.isnan(v), np.nan, p)
    return pd.Series(p if higher_better else 1 - p, index=s.index)


def edit(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[-1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def build(args) -> pd.DataFrame:
    D = pd.DataFrame(dn.load_designs())
    if D.empty:
        sys.exit("no designs yet")
    # framework-diversification variants and BoltzGen designs were cancelled/parked by the user: never ranked
    D = D[~D.generator.isin(["fwdiv", "boltzgen", "boltzgen-cys"])].copy()
    D["is_control"] = (D.generator == "control").astype(int)
    for name, keep in (("rules.csv", None), ("score_of3.csv", None), ("score_orth.csv", None), ("score_ptx.csv", None),
                       ("props.csv", None), ("designs_boltzgen_metrics.csv", None)):
        t = rd(name)
        if t is None or "design_id" not in t:
            continue
        t = t.drop_duplicates("design_id", keep="last")
        drop = [c for c in t.columns if c in D.columns and c != "design_id"]
        if name == "score_orth.csv":
            t = t.rename(columns={c: "orth_" + c for c in t.columns if c != "design_id"})
            drop = []
        D = D.merge(t.drop(columns=drop), on="design_id", how="left")
    D = D.copy()  # defragment after the many merges
    ctl = D.is_control == 1
    pop = D[~ctl]

    def col(name):
        return D[name] if name in D else pd.Series(np.nan, index=D.index)

    # ---- binding evidence
    if "of3_consensus_rank_score" in D and pop["of3_consensus_rank_score"].notna().any():
        bind = col("of3_consensus_rank_score")
    else:
        parts = [pct(col("of3_ipsae"), pop.get("of3_ipsae", pd.Series(dtype=float))),
                 pct(col("of3_iptm"), pop.get("of3_iptm", pd.Series(dtype=float))),
                 pct(col("of3_min_ipae"), pop.get("of3_min_ipae", pd.Series(dtype=float)), higher_better=False)]
        bind = pd.concat(parts, axis=1).mean(axis=1, skipna=True)
    # orthogonal evidence = mean of the available orthogonal scores: S3's OpenDDE/ESMFold-2 consensus and Nebius Protenix-v2
    # (`ptx_score`, worker PX, results/denovo/score_ptx.csv), both within-run percentile-style scores in 0..1
    orth = pd.concat([col("orth_consensus_rank_score"), col("ptx_score")], axis=1).mean(axis=1, skipna=True)
    D["t_bind"] = np.where(orth.notna() & bind.notna(), 0.6 * bind + 0.4 * orth, bind)
    comp = pd.concat([pct(col("of3_epi_overlap"), pop.get("of3_epi_overlap", pd.Series(dtype=float))),
                      pct(col("of3_cdr_fv_frac4"), pop.get("of3_cdr_fv_frac4", pd.Series(dtype=float)))], axis=1)
    D["t_comp"] = comp.mean(axis=1, skipna=True)
    D["t_therm"] = pct(col("nm_tm"), pop.get("nm_tm", pd.Series(dtype=float)))
    D["t_yield"] = pct(col("my_yield"), pop.get("my_yield", pd.Series(dtype=float)))
    # developability: sequence liabilities exist for EVERY design (S4 `liab_n_flags` = number of flagged liabilities: free Cys,
    # N-glyc, NG/DG, Met/Trp in CDRs, hydrophobic run, pI, CDR charge), so the objective never depends on TNP being run.
    # TNP red flags, when available, are an extra penalty (0.5 each). (An earlier version counted every numeric liab_* column,
    # including measurements such as length and pI, and only existed for designs with TNP: both were wrong.)
    nflag = col("liab_n_flags")
    D["t_dev"] = np.where(nflag.notna(), (1 - 0.15 * nflag.fillna(0) - 0.5 * col("tnp_n_red").fillna(0)).clip(lower=0), np.nan)
    # nine-rule CDR2/CDR3 sequence score from the monoyield paper (Fig. 5h), as a small yield prior: nb9_score / nb9_applicable
    D["t_rules"] = (col("nb9_score") / col("nb9_applicable").replace(0, np.nan)).clip(0, 1)
    W = {"t_bind": args.w_bind, "t_comp": args.w_comp, "t_therm": args.w_therm, "t_yield": args.w_yield, "t_dev": args.w_dev,
         "t_rules": args.w_rules}
    # Missing evidence counts as a below-median prior, not as "perfect": a design measured on one term only (e.g. just the
    # nine-rule score) must not top the list. `coverage` = share of the weight actually backed by data.
    PRIOR = 0.3
    num = sum(D[k].fillna(PRIOR) * w for k, w in W.items())
    den = sum(D[k].notna().astype(float) * w for k, w in W.items())
    D["composite"] = np.where(den > 0, num / sum(W.values()), np.nan)
    D["coverage"] = den / sum(W.values())

    # ---- gates
    one = lambda n: col(n).fillna(0).astype(float) > 0  # noqa: E731
    # hard rules: length, no unpaired Cys (odd Cys count, abguard's definition), numbering, AND the canonical VHH disulfide
    # (Cys at IMGT 23 and 104): a design with zero Cys passes abguard trivially but loses a lot of Tm (BoltzGen's
    # nanobody protocol designs with Cys excluded, so its raw designs have none).
    hard = (one("len_ok") & (col("unpaired_cys").fillna(0) == 0) & (one("numbering_ok") | col("numbering_ok").isna())
            & (one("canon_cys_ok") | col("canon_cys_ok").isna()))
    ident_ok = one("ident_lt70") if "ident_lt70" in D else pd.Series(False, index=D.index)
    abg = one("abguard_pass")
    # S2 columns win when present: `denovo_vm_pass` = VM abguard says de_novo (whole-Fv identity to trastuzumab VH < 70%),
    # `denovo_strict_pass` = that AND < 70% to every known VHH/heavy chain in S2's panel. `abguard_pass` alone only
    # means "valid antibody", not "de novo".
    # Organisers (via the user, 14:15): fully new CDRs on an intact framework should land under 70%, so the framework
    # is NOT diversified. Gate = VM abguard says de_novo (whole-Fv identity to trastuzumab VH < 70%); `gate_strict` adds a
    # 2-point safety margin (<= 0.68) because the jury's exact identity definition is unconfirmed. Identity to the nearest
    # known VHH in S2's panel (~0.80 for RFD4 designs on the intact framework) is reported but is not a gate.
    vm = one("denovo_vm_pass") if "denovo_vm_pass" in D else (hard & abg)
    idt = col("identity_to_trastuzumab")
    margin_ok = (idt <= 0.68) | idt.isna()
    D["gate_vm"] = (hard & vm).astype(int)
    D["gate_strict"] = (hard & vm & margin_ok).astype(int)
    D["panel_ident_lt70"] = ident_ok.astype(int)
    D["needs_fwdiv"] = 0

    def reason(r):
        bits = []
        if pd.notna(r.get("of3_iptm")):
            bits.append(f"OF3 ipTM {r['of3_iptm']:.2f}, ipSAE {r.get('of3_ipsae', np.nan):.2f}")
        if pd.notna(r.get("of3_epi_overlap")):
            bits.append(f"epitope overlap {r['of3_epi_overlap']:.0%}")
        if pd.notna(r.get("orth_n_models_pass")):
            bits.append(f"{int(r['orth_n_models_pass'])}/{int(r.get('orth_n_models_scored', 0))} orthogonal models agree")
        if pd.notna(r.get("nm_tm")):
            bits.append(f"NanoMelt Tm {r['nm_tm']:.1f} C")
        if pd.notna(r.get("my_yield")):
            bits.append(f"yield {r['my_yield']:.2f}")
        if pd.notna(r.get("tnp_n_red")):
            bits.append(f"TNP red flags {int(r['tnp_n_red'])}")
        if pd.notna(r.get("identity_nearest_known")):
            bits.append(f"nearest known VHH {r['identity_nearest_known']:.0%}")
        if pd.notna(r.get("identity_to_trastuzumab")):
            bits.append(f"{r['identity_to_trastuzumab']:.0%} to trastuzumab VH")
        if pd.notna(r.get("nb9_score")):
            bits.append(f"nine-rule score {int(r['nb9_score'])}/{int(r.get('nb9_applicable', 9))}")
        bits.append("abguard de novo" if r["gate_vm"] else "abguard NOT de novo")
        return "; ".join(bits)

    D["rank_reason"] = D.apply(reason, axis=1)
    D = D.sort_values(["is_control", "gate_strict", "composite"], ascending=[True, False, False], na_position="last")
    D["rank_all"] = np.where(D.is_control == 1, np.nan, np.arange(1, len(D) + 1))
    return D


def pick_top(D: pd.DataFrame, n: int, min_cdr3_dist: int = 3) -> pd.DataFrame:
    """Greedy: best composite first, skip near-duplicate CDR3s (edit distance < 3), and keep diversity across
    framework x hotspot-set groups (<= 50% of the list from one group), since RFD4 is now the only generator."""
    pool = D[(D.is_control == 0) & D.composite.notna()]
    strict = pool[pool.gate_strict == 1]
    note = "strict gate (abguard + identity < 70%)"
    if len(strict) < n:
        extra = pool[(pool.gate_strict == 0) & (pool.gate_vm == 1)]
        pool = pd.concat([strict, extra.sort_values("composite", ascending=False)])
        note = f"only {len(strict)} pass the strict identity gate; remainder are VM-abguard-only (need framework diversification)"
    chosen, cdrs = [], []
    cap = max(1, int(np.ceil(0.5 * n)))
    count = {}
    for _, r in pool.iterrows():
        c3 = str(r.get("cdr3", ""))
        if any(edit(c3, o) < min_cdr3_dist for o in cdrs):
            continue
        grp = f"{r.get('framework', '')}|{r.get('hotspot_set', '')}"
        if count.get(grp, 0) >= cap:
            continue
        chosen.append(r)
        cdrs.append(c3)
        count[grp] = count.get(grp, 0) + 1
        if len(chosen) == n:
            break
    S = pd.DataFrame(chosen)
    if not S.empty:
        S.insert(0, "slot", np.arange(1, len(S) + 1))
        S["gate_note"] = note
    return S


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--w-bind", type=float, default=0.40)
    ap.add_argument("--w-comp", type=float, default=0.15)
    ap.add_argument("--w-therm", type=float, default=0.20)
    ap.add_argument("--w-yield", type=float, default=0.15)
    ap.add_argument("--w-dev", type=float, default=0.10)
    ap.add_argument("--w-rules", type=float, default=0.05)
    a = ap.parse_args()
    D = build(a)
    keep = ["rank_all", "design_id", "parent_id", "generator", "framework", "hotspot_set", "gate_strict", "gate_vm",
            "needs_fwdiv", "composite", "coverage", "t_bind", "t_comp", "t_therm", "t_yield", "t_dev", "cdr3", "vhh_seq",
            "rank_reason"]
    cols = [c for c in keep if c in D] + [c for c in D.columns if c not in keep]
    D[cols].to_csv(OUT / "leaderboard.csv", index=False)
    S = pick_top(D, a.top)
    if not S.empty:
        S[["slot"] + [c for c in cols if c in S]].to_csv(OUT / "submission_top20.csv", index=False)
    n = int((D.is_control == 0).sum())
    print(f"leaderboard: {n} designs, strict-gate {int(D.gate_strict.sum())}, vm-gate {int(D.gate_vm.sum())}, "
          f"scored composite {int(D.composite.notna().sum())}; submission rows {len(S)}")
    if not S.empty:
        print(S[["slot", "design_id", "generator", "composite", "gate_strict"]].head(10).to_string(index=False))
