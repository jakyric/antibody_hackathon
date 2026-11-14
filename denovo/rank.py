"""S5: merge every stage into one leaderboard + a diversity-aware top-20 submission list.

    baseline/.venv/bin/python denovo/rank.py            # writes results/denovo/leaderboard.csv, submission_top20.csv
    baseline/.venv/bin/python denovo/rank.py --top 20 --w-bind 0.4 ...

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
LEN_MIN, LEN_MAX = 12, 18  # allowed CDR3 length, IMGT block incl. the 'AA' after Cys104 (Kabat-like loop = this minus 2); tightened 20 -> 18 by the user at 15:20


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
    for name, keep in (("rules.csv", None), ("score_of3.csv", None), ("score_orth.csv", None), ("score_ptx.csv", None), ("score_vmof3.csv", None), ("score_dde_extra.csv", None),
                       ("props.csv", None), ("props_nm.csv", None), ("props_tnpvm.csv", None), ("validity.csv", None), ("gapcheck.csv", None), ("designs_boltzgen_metrics.csv", None)):
        t = rd(name)
        if t is None or "design_id" not in t:
            continue
        t = t.drop_duplicates("design_id", keep="last")
        if name in ("props_nm.csv", "props_tnpvm.csv"):  # NanoMelt (local CPU) and TNP (VM CPU) run in their own processes/files
            tt = t.set_index("design_id")                # so they cannot slow or race the main TNP/monoyield loop; merged by fill
            pre = "nm_" if name == "props_nm.csv" else "tnp_"
            for c in [c for c in tt.columns if c.startswith(pre)]:
                m = D.design_id.map(tt[c])
                D[c] = D[c].fillna(m) if c in D else m
            continue
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

    # ---- binding evidence = iPTM and ipSAE of EACH scorer (user, 15:20: keep it simple, binding first).
    # For every scorer that has data for >= 5 designs: within-run percentile of its iPTM and of its ipSAE (chain-length-scaled,
    # 15 A cutoff: `ipsae_chn15`, the variant that separates the trastuzumab control from the decoys); t_bind = mean of all those
    # percentiles. Scorers: OpenFold3 (2 seeds), OpenDDE-abag, ESMFold-2, AF2-IG, Protenix-v2.
    SCORERS = {"of3": ("of3_iptm", "of3_ipsae_chn15"), "of3s1": ("orth_of3s1_iptm", "orth_of3s1_ipsae_chn15"),
               "dde": ("orth_dde_iptm", "orth_dde_ipsae_chn15"), "esm2": ("orth_esm2_iptm", "orth_esm2_ipsae_chn15"),
               "af2ig": ("orth_af2ig_iptm", "orth_af2ig_ipsae_chn15"), "ptx": ("ptx_iptm", None),
               # VM GPU runs (user, 15:30): OpenFold3 OpenBind with the HER2 template and no MSA (the upstream-style recipe) and
               # untemplated ESMFold-2; pair_iptm = VHH-HER2 chain-pair ipTM (the whole-complex iptm also counts HER2 fragment-fragment)
               "vmof3": ("vmof3_pair_iptm", "vmof3_ipsae_chn15"), "vmesm2": ("vmesm2_pair_iptm", "vmesm2_ipsae_chn15")}
    parts, iptm_cols, ipsae_cols = [], [], []
    for _name, (ci, cs) in SCORERS.items():
        for c, lst in ((ci, iptm_cols), (cs, ipsae_cols)):
            if c and c in D and pop[c].notna().sum() >= 5:
                parts.append(pct(col(c), pop[c]).rename("p_" + c))
                lst.append(c)
    D["t_bind"] = pd.concat(parts, axis=1).mean(axis=1, skipna=True) if parts else np.nan
    D["iptm_mean"] = D[iptm_cols].mean(axis=1) if iptm_cols else np.nan      # raw mean iPTM over the scorers that ran
    D["ipsae_mean"] = D[ipsae_cols].mean(axis=1) if ipsae_cols else np.nan
    D["n_scorers"] = D[iptm_cols].notna().sum(axis=1) if iptm_cols else 0
    comp = pd.concat([pct(col("of3_epi_overlap"), pop.get("of3_epi_overlap", pd.Series(dtype=float))),
                      pct(col("of3_cdr_fv_frac4"), pop.get("of3_cdr_fv_frac4", pd.Series(dtype=float)))], axis=1)
    D["t_comp"] = comp.mean(axis=1, skipna=True)
    D["t_therm"] = pct(col("nm_tm"), pop.get("nm_tm", pd.Series(dtype=float)))
    # TNP axis (user, 15:20): developability = TNP only, 1 - 0.5 per red flag - 0.2 per amber flag; NaN where TNP has not run.
    # (Sequence liabilities and NanoMelt Tm stay as plain columns, they are not scored.)
    D["t_dev"] = np.where(col("tnp_n_red").notna(), (1 - 0.5 * col("tnp_n_red").fillna(0) - 0.2 * col("tnp_n_amber").fillna(0)).clip(lower=0), np.nan)
    # yield/rules axis: mean of the within-run percentile of monoyield (ESMC-Nb) and the nine-rule CDR2/CDR3 score of the
    # monoyield paper (nb9_score / applicable rules); the nine-rule fraction alone where monoyield has not run.
    D["t_rules"] = (col("nb9_score") / col("nb9_applicable").replace(0, np.nan)).clip(0, 1)
    D["t_yield"] = pd.concat([pct(col("my_yield"), pop.get("my_yield", pd.Series(dtype=float))), D["t_rules"]], axis=1).mean(axis=1, skipna=True)
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
    # structural validity (worker VAL, results/denovo/validity.csv): severe flags invalidate a design outright. The generator
    # sometimes fails to hold its conditioning (kept HER2+framework group moved 13 A, `motif_bad`), or table and structure
    # disagree. Mild flags (clashes, stretched disulfide, noisy bond lengths) are NOT gates: raw decoder output is not idealised.
    sev_cols = ["val_parse_fail", "val_chain_layout", "val_seq_mismatch", "val_fr_not_intact", "val_cdr_table_mismatch",
                "val_her2_seq_changed", "val_fw_seq_changed", "val_motif_bad", "val_incomplete", "val_no_contact",
                "val_clash_severe", "val_cys_odd"]
    sev = sum((col(c).fillna(0).astype(float) for c in sev_cols if c in D), pd.Series(0.0, index=D.index))
    D["val_severe"] = (sev > 0).astype(int)
    hard = hard & (sev == 0)
    # reasonable CDR3 length (user, 15:12): the IMGT block incl. the 'AA' after Cys104 must be LEN_MIN..LEN_MAX aa (a Kabat-like
    # loop of 10-18). References: BoltzGen designs realised a median of 15 (IMGT block, range 7-23); the monoyield paper's
    # BoltzGen/RFantibody libraries on the same scaffold used 5-13 (Kabat-like, i.e. 7-15 IMGT block).
    c3len = col("cdr3").astype(str).str.len()
    D["cdr3_block_len"] = c3len
    D["len_reasonable"] = ((c3len >= LEN_MIN) & (c3len <= LEN_MAX)).astype(int)
    hard = hard & (D.len_reasonable == 1)
    # no homopolymer run of 6+ residues in the CDRs (VAL: e.g. SSSSSSSSSS), and no conflict with the unresolved HER2 loop
    # 581-590 (gapcheck.csv: >= 1 VHH heavy atom within 3 A of the loop-anchor atoms, closer than trastuzumab's Fv ever gets, 4.7 A)
    cdrs = col("cdr1").astype(str) + "|" + col("cdr2").astype(str) + "|" + col("cdr3").astype(str)
    D["cdr_run6"] = cdrs.str.contains(r"(.)\1{5,}", regex=True).astype(int)
    D["gap_conflict"] = (col("gap_anchor_atoms_lt3").fillna(0) > 0).astype(int)
    hard = hard & (D.cdr_run6 == 0) & (D.gap_conflict == 0)

    # ---- review FLAGS: soft items that never exclude a design but order and annotate it. Typed rules: a `gate` = eligibility (a missing value FAILS), a `flag` = review item (a missing value does not
    # trigger). Flag thresholds: the tail of the decoy distribution (95th percentile of the framework-only decoys, per scorer metric).
    decoys = D[D.design_id.str.startswith("ctl-fwrand") | (D.design_id == "ctl-bcii10-wt")]
    below = []
    for c in iptm_cols + ipsae_cols:
        ref = decoys[c].dropna() if c in decoys else pd.Series(dtype=float)
        if len(ref) >= 3:
            below.append((col(c) <= ref.quantile(0.95)).where(col(c).notna()))
    if below:
        B = pd.concat(below, axis=1)
        n_meas = B.notna().sum(axis=1)
        D["n_below_decoy"] = B.sum(axis=1, min_count=1)
        D["flag_below_decoy"] = ((n_meas > 0) & (D["n_below_decoy"] / n_meas.replace(0, np.nan) >= 0.5)).astype(int)
    else:
        D["n_below_decoy"], D["flag_below_decoy"] = np.nan, 0
    FLAGS = {"below_decoy": D["flag_below_decoy"] == 1, "nglyc_cdr": col("liab_nglyc_cdr").fillna(0) > 0,
             "low_complexity": one("val_low_complexity"), "cys_in_cdr": one("val_cys_unexpected"), "clash": one("val_clash_bad"),
             "tnp_red": col("tnp_n_red").fillna(0) > 0, "tnp_amber2": col("tnp_n_amber").fillna(0) >= 2, "rules_low": D["t_rules"] < 0.5,
             # OpenDDE binder self-consistency RMSD (design VHH vs predicted VHH, scRMSD) > 1.5 A: the self-consistency flag threshold
             "dde_fold_rmsd": col("orth_dde_own_vhh_rmsd_fold") > 1.5}
    D["n_flags"] = sum(v.astype(int) for v in FLAGS.values())
    D["flags"] = [";".join(k for k, v in FLAGS.items() if v.iloc[i]) for i in range(len(D))]
    # evidence completeness (order: eligible -> complete evidence -> Pareto tier -> fewest flags -> rank)
    me = pd.DataFrame({"binding(<2 scorers)": D["n_scorers"].fillna(0) < 2, "tnp": col("tnp_n_red").isna(), "monoyield": col("my_yield").isna()})
    D["missing_evidence"] = me.apply(lambda r: ",".join(k for k, v in r.items() if v), axis=1)
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
    # OpenDDE binding GATE (user, 15:44): OpenDDE-abag iPTM AND standard ipSAE (10 A cutoff) must both reach the 95th percentile of the
    # framework-only decoys (adverse-tail logic; thresholds follow the controls, fallbacks 0.49 / 0.12). A design
    # without an OpenDDE score FAILS the gate (a gate with a missing value fails).
    _dec = D[D.design_id.str.startswith("ctl-fwrand") | (D.design_id == "ctl-bcii10-wt")]
    def _q95(c, fallback):
        s = _dec[c].dropna() if c in _dec else pd.Series(dtype=float)
        return float(s.quantile(0.95)) if len(s) >= 3 else fallback
    D["dde_thr_iptm"], D["dde_thr_ipsae"] = _q95("orth_dde_iptm", 0.49), _q95("orth_dde_ipsae", 0.12)
    D["gate_dde"] = ((col("orth_dde_iptm") >= D["dde_thr_iptm"]) & (col("orth_dde_ipsae") >= D["dde_thr_ipsae"])).astype(int)
    D["gate_vm"] = (hard & vm).astype(int)
    D["gate_strict"] = (hard & vm & margin_ok & (D.gate_dde == 1)).astype(int)
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
        if isinstance(r.get("flags"), str) and r["flags"]:
            bits.append("flags: " + r["flags"].replace(";", ", "))
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
    ap.add_argument("--w-bind", type=float, default=0.50)   # iPTM / ipSAE of each scorer
    ap.add_argument("--w-comp", type=float, default=0.0)    # epitope/competition proxy: shown, not scored
    ap.add_argument("--w-therm", type=float, default=0.20)  # NanoMelt Tm (added back by the user at 15:44)
    ap.add_argument("--w-yield", type=float, default=0.15)  # monoyield + nine rules
    ap.add_argument("--w-dev", type=float, default=0.15)    # TNP flags
    ap.add_argument("--w-rules", type=float, default=0.0)   # folded into the yield axis
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
