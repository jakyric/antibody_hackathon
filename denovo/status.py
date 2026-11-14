"""Funnel view of the continuously running pipeline: how many designs reached each stage, per generator,
how the controls (known positives / decoys) sit relative to the designs, and the JapanFold balance.

    baseline/.venv/bin/python denovo/status.py        # prints and writes results/denovo/status.md
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dn_common as dn  # noqa: E402

OUT = dn.OUT


def md(df, index=False):
    """Tiny markdown table (no tabulate dependency)."""
    if index:
        df = df.reset_index()
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for row in df.itertuples(index=False):
        out.append("| " + " | ".join("" if (isinstance(v, float) and v != v) else str(v) for v in row) + " |")
    return "\n".join(out)


SKIP = {"rank_all", "is_control", "seed", "numbering_ok", "imm_chain"}


def md_top(top):
    """Top designs with ALL numeric metrics: metrics as rows, one column per design (readable when ~90 metrics)."""
    cols = [c for c in top.columns if c not in SKIP and pd.api.types.is_numeric_dtype(top[c])]
    # keep the families together: scores first, then per-stage metrics
    order = ["composite", "coverage"] + [c for c in cols if c.startswith("t_")]
    for pre in ("of3_", "orth_", "nm_", "my_", "tnp_", "liab_", "nb"):
        order += [c for c in cols if c.startswith(pre) and c not in order]
    order += [c for c in cols if c not in order]
    short = [str(r.design_id).replace("-rfd4jf", "") for r in top.itertuples()]
    rows = ["| metric | " + " | ".join(short) + " |", "|---|" + "---|" * len(short)]
    for c in order:
        vals = ["" if pd.isna(v) else (f"{v:.3g}" if isinstance(v, float) else str(v)) for v in top[c]]
        rows.append(f"| {c} | " + " | ".join(vals) + " |")
    return "\n".join(rows)


def md_focus(top):
    """Top designs, one row each, only the three axes the user asked for: iPTM/ipSAE of each scorer, TNP, yield/rules."""
    cols = ["design_id", "framework", "cdr3_block_len", "composite", "t_bind", "iptm_mean", "ipsae_mean", "of3_iptm", "orth_dde_iptm",
            "orth_esm2_iptm", "orth_af2ig_iptm", "ptx_iptm", "of3_ipsae_chn15", "orth_dde_ipsae_chn15", "orth_esm2_ipsae_chn15",
            "orth_af2ig_ipsae_chn15", "vmof3_pair_iptm", "vmof3_ipsae_chn15", "vmesm2_pair_iptm", "vmesm2_ipsae_chn15", "ddex_pdockq2", "orth_dde_own_vhh_rmsd_fold", "orth_dde_own_vhh_rmsd_her2frame", "tnp_n_red", "tnp_n_amber", "my_yield", "nm_tm", "nb9_score", "nb9_applicable", "gate_dde", "gate_strict", "n_flags", "flags"]
    t = top[[c for c in cols if c in top.columns]].copy()
    t["design_id"] = t.design_id.str.replace("-rfd4jf", "", regex=False)
    t = t.rename(columns={"cdr3_block_len": "cdr3", "orth_dde_iptm": "dde_iptm", "orth_esm2_iptm": "esm2_iptm", "orth_af2ig_iptm": "af2ig_iptm",
                          "of3_ipsae_chn15": "of3_ipsae", "orth_dde_ipsae_chn15": "dde_ipsae", "orth_esm2_ipsae_chn15": "esm2_ipsae",
                          "orth_af2ig_ipsae_chn15": "af2ig_ipsae", "nb9_applicable": "nb9_of"})
    return md(t.round(2))


def rd(name):
    p = OUT / name
    try:
        return pd.read_csv(p) if p.exists() and p.stat().st_size else None
    except Exception:
        return None


def main():
    D = pd.DataFrame(dn.load_designs())
    if not D.empty:  # parked/cancelled generators never appear in the overview
        D = D[~D.generator.isin(["fwdiv", "boltzgen", "boltzgen-cys"])]
    lines =[f"# Pipeline status, {time.strftime('%H:%M:%S %Z')}", ""]
    if D.empty:
        lines.append("no designs yet")
    else:
        L = rd("leaderboard.csv")
        R, S3, SO, P = rd("rules.csv"), rd("score_of3.csv"), rd("score_orth.csv"), rd("props.csv")
        real = D[D.generator != "control"]
        ids = lambda t: set(t.design_id) if t is not None and "design_id" in t else set()  # noqa: E731
        rows = []
        for g, sub in real.groupby("generator"):
            s = set(sub.design_id)
            row = {"generator": g, "designs": len(s), "rules": len(s & ids(R)), "of3": len(s & ids(S3)),
                   "orth": len(s & ids(SO)), "props": len(s & ids(P))}
            if L is not None:
                l = L[L.design_id.isin(s)]
                row["gate_vm"] = int(l.gate_vm.sum()) if "gate_vm" in l else 0
                row["gate_strict"] = int(l.gate_strict.sum()) if "gate_strict" in l else 0
                row["composite"] = int(l.composite.notna().sum()) if "composite" in l else 0
            rows.append(row)
        T = pd.DataFrame(rows)
        lines += ["## Funnel (designs that reached each stage)", "", md(T) if len(T) else "", ""]
        if S3 is not None and "of3_iptm" in S3:
            m = S3.merge(D[["design_id", "generator"]], on="design_id", how="inner")  # only designs/controls still in play
            m["kind"] = m.design_id.str.startswith("ctl-").map({True: "control", False: "design"})
            cols = [c for c in ("of3_iptm", "of3_ipsae", "of3_min_ipae", "of3_epi_overlap") if c in m]
            lines += ["## OF3 metrics: controls vs designs (median)", "", md(m.groupby("kind")[cols].median().round(3), index=True), ""]
            c = m[m.kind == "control"][["design_id"] + cols].round(3)
            lines += ["### Control rows", "", md(c), ""]
        if L is not None and "composite" in L:
            top = L[(L.is_control == 0) & L.composite.notna()].head(10)
            if len(top):
                lines += ["## Current top 10 by composite", "",
                          md_focus(top), ""]
            P = rd("leaderboard_pareto.csv")
            if P is not None and len(P):
                fr = P.pareto_rank.value_counts().sort_index().head(4)
                lines += ["## Pareto (last step, `denovo/pareto.py`)", "",
                          f"Population {len(P)} gate-passing designs; objectives {P.pareto_objectives.iloc[0]}; "
                          f"front sizes { {int(k): int(v) for k, v in fr.items()} }; fully measured {int(P.pareto_complete.sum())}.", ""]
    try:
        import json
        import urllib.request
        b = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
        r = urllib.request.Request(b + "/v1/credits", headers={"Authorization": "Bearer " + os.environ["JAPANFOLD_API_KEY"],
                                                               "X-JapanFold-Client": "skill"})
        lines.append(f"JapanFold balance: ${json.load(urllib.request.urlopen(r, timeout=20))['balance']:.2f}")
    except Exception as e:
        lines.append(f"(balance unavailable: {type(e).__name__})")
    txt = "\n".join(lines)
    (OUT / "status.md").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
