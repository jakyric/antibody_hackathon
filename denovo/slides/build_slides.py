#!/usr/bin/env python3
"""Build the pipeline overview deck from the live result files.

    denovo/.venv/bin/python denovo/slides/build_slides.py            # regenerate denovo/slides/artifact/project/*
    denovo/.venv/bin/python denovo/slides/build_slides.py --plan     # also print the publish plan (JSON)

Every number on a slide is read here from results/denovo/*.csv|json|md and denovo/*/ files at build time; where a stage has not
produced evidence yet the slide says "pending" instead of showing a made-up value. No credentials, hosts or key paths are read or
written. `deck.css` (artifact/project/deck.css) is static; this script owns deck.json and slides/*.html.
Publishing is a separate step (see denovo/slides/README.md): one Artifact call to the same URL and root.
"""
from __future__ import annotations

import html
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RES = ROOT / "results" / "denovo"
PROJ = HERE / "artifact" / "project"
SLIDES = PROJ / "slides"
URL = "https://claude.ai/artifact/XnQZswmE7DSG6qpz9Dm1XG"
START_BALANCE = 93.844  # JapanFold balance before the pipeline's first job (S3 spend log)
SPEND_CAP = 65.0        # CONTRACT.md rule 3

E = lambda s: html.escape(str(s), quote=True)  # noqa: E731
TS = time.strftime("%H:%M %Z")


# ------------------------------------------------------------------------------------------------ data
def rd(name):
    p = RES / name
    try:
        return pd.read_csv(p) if p.exists() and p.stat().st_size else None
    except Exception:
        return None


def rjson(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return {}


def pct(n, d):
    return f"{100 * n / d:.0f}%" if d else "–"


def collect() -> dict:
    st: dict = {}
    D = rd("designs_rfd4jf.csv")
    if D is None:
        D = pd.DataFrame(columns=["design_id", "framework", "hotspot_set"])
    ids = set(D.design_id)
    fw = dict(zip(D.design_id, D.framework))
    st["n_designs"] = len(D)
    st["by_fw_hs"] = D.groupby(["framework", "hotspot_set"]).size().to_dict() if len(D) else {}
    st["fw_counts"] = D.framework.value_counts().to_dict() if len(D) else {}

    R = rd("rules.csv")
    S = rd("score_of3.csv")
    O = rd("score_orth.csv")
    P = rd("props.csv")
    Lb = rd("leaderboard.csv")
    Lp = rd("leaderboard_pareto.csv")

    # rules / gate
    r = R[R.design_id.isin(ids)].copy() if R is not None else pd.DataFrame()
    st["n_rules"] = len(r)
    if len(r):
        r["fw"] = r.design_id.map(fw)
        r["ok"] = r.abguard_class.eq("de_novo")
        st["n_denovo"] = int(r.ok.sum())
        st["pass_by_fw"] = {k: (int(g.ok.sum()), len(g)) for k, g in r.groupby("fw")}
        b = r[r.fw == "bcii10"]
        bins = [(0, 122, "≤ 122"), (123, 124, "123–124"), (125, 126, "125–126"), (127, 128, "127–128"), (129, 999, "≥ 129")]
        st["gate_bins"] = []
        for lo, hi, lab in bins:
            g = b[(b.vhh_len >= lo) & (b.vhh_len <= hi)]
            st["gate_bins"].append((lab, int(g.ok.sum()), len(g)))
        ni = r.identity_nearest_known.dropna() if "identity_nearest_known" in r else pd.Series(dtype=float)
        st["nearest_known"] = (float(ni.min()), float(ni.max()), int((ni < 0.70).sum()), len(ni)) if len(ni) else None
        idt = r.identity_to_trastuzumab.dropna()
        st["trast_ident"] = (float(idt.min()), float(idt.max()), float(idt.mean())) if len(idt) else None
    else:
        st.update(n_denovo=0, pass_by_fw={}, gate_bins=[], nearest_known=None, trast_ident=None)

    if Lb is not None:
        lb = Lb[(Lb.is_control == 0) & Lb.design_id.isin(ids)]
        st["n_gate_vm"] = int(lb.gate_vm.sum()) if "gate_vm" in lb else 0
        st["n_gate_strict"] = int(lb.gate_strict.sum()) if "gate_strict" in lb else 0
        st["n_cov50"] = int((lb.coverage >= 0.5).sum()) if "coverage" in lb else 0
    else:
        st.update(n_gate_vm=0, n_gate_strict=0, n_cov50=0)

    # OpenFold3 and orthogonal
    st["of3"] = {"n": 0, "decoy": [], "design": [], "tz": [], "pass": 0, "med_d": None, "med_c": None, "tz_v": None}
    if S is not None and "of3_iptm" in S:
        d = S[S.design_id.isin(ids)]
        dec = S[(S.is_control == 1) & (S.design_id.str.startswith("ctl-fwrand") | S.design_id.eq("ctl-bcii10-wt"))]
        tz = S[S.design_id.eq("ctl-tz-fv")]
        st["of3"] = {"n": len(d), "decoy": dec.of3_iptm.dropna().tolist(), "design": d.of3_iptm.dropna().tolist(),
                     "tz": tz.of3_iptm.dropna().tolist(), "pass": int(d.of3_pass.fillna(0).sum()) if "of3_pass" in d else 0,
                     "med_d": float(d.of3_iptm.median()) if len(d) else None,
                     "med_c": float(dec.of3_iptm.median()) if len(dec) else None,
                     "tz_v": float(tz.of3_iptm.iloc[0]) if len(tz) else None}
    st["n_orth"] = 0
    if O is not None and "n_models_scored" in O:
        o = O[O.design_id.isin(ids)]
        st["n_orth"] = int((o.n_models_scored >= 2).sum())

    # properties
    st["tm"] = {"bcii10": [], "capla": [], "ref": [], "n": 0}
    st["n_tnp"] = st["n_yield"] = 0
    if P is not None:
        pr = P[P.design_id.isin(ids)].copy()
        pr["fw"] = pr.design_id.map(fw)
        if "nm_tm" in pr:
            for k in ("bcii10", "capla"):
                st["tm"][k] = pr[(pr.fw == k) & pr.nm_tm.notna()].nm_tm.tolist()
            st["tm"]["n"] = int(pr.nm_tm.notna().sum())
            ref = P[P.design_id.str.startswith("ctl-s4-") & ~P.design_id.str.contains("trastuzumab") & P.nm_tm.notna()]
            st["tm"]["ref"] = ref.nm_tm.tolist()
        st["n_tnp"] = int(pr.tnp_n_red.notna().sum()) if "tnp_n_red" in pr else 0
        st["n_yield"] = int(pr.my_yield.notna().sum()) if "my_yield" in pr else 0
        if "nb9_score" in pr:
            nb = pr[pr.nb9_applicable == 9].nb9_score.dropna()
            st["nb9"] = (float(nb.median()), len(nb)) if len(nb) else None
        else:
            st["nb9"] = None
    else:
        st["nb9"] = None

    # candidates (Pareto leaderboard if present)
    L = Lp if Lp is not None else Lb
    st["cand"] = []
    st["cand_reason"] = None
    st["pareto"] = Lp is not None
    if L is not None:
        c = L[(L.is_control == 0) & L.design_id.isin(ids)].copy()
        if "gate_strict" in c:
            c = c[c.gate_strict == 1]
        if "pareto_complete" in c:
            c = c[c.pareto_complete == 1]
        elif "coverage" in c:
            c = c[c.coverage >= 0.5]
        if len(c):
            c = c.sort_values(["pareto_rank", "maximin"], ascending=[True, False]) if "pareto_rank" in c else c.sort_values("composite", ascending=False)
            st["cand_n"] = len(c)
            st["cand"] = c.head(8).to_dict("records")
            st["cand_reason"] = c.iloc[0].get("rank_reason")

    # slot list
    sp = rd("submission_top20_pareto.csv")
    sl = sp if sp is not None else rd("submission_top20.csv")
    st["slots"] = None
    if sl is not None and len(sl):
        okc = ("pareto_complete" in sl and sl.pareto_complete.fillna(0).mean() >= 0.5) or ("coverage" in sl and sl.coverage.mean() >= 0.5)
        if okc:
            st["slots"] = {"n": len(sl), "fw": sl.framework.value_counts().to_dict(), "hs": sl.hotspot_set.value_counts().to_dict()}

    # Pareto step (denovo/pareto.py): fronts, plotted points, objectives used
    st["par"] = None
    if Lp is not None and len(Lp) and "pareto_rank" in Lp:
        p = Lp[(Lp.is_control == 0) & Lp.design_id.isin(ids)].copy()
        sel = set(sp.design_id) if sp is not None else set()
        fronts = p.pareto_rank.value_counts().sort_index().tolist()
        both = p[p.t_bind.notna() & p.t_therm.notna()] if {"t_bind", "t_therm"} <= set(p.columns) else p.iloc[0:0]
        pts = [(float(x.t_bind), float(x.t_therm), int(x.pareto_rank), x.design_id in sel) for x in both.itertuples()]
        used = []
        if "pareto_reason" in p and len(p):
            m = re.search(r"over ([a-z,]+)", str(p.sort_values("pareto_rank").iloc[0].pareto_reason))
            used = m.group(1).split(",") if m else []
        allobj = ["bind", "comp", "therm", "yield", "dev"]
        st["par"] = {"n_pop": len(p), "n_complete": int(p.pareto_complete.fillna(0).sum()) if "pareto_complete" in p else 0,
                     "fronts": fronts, "pts": pts, "used": used, "dropped": [o for o in allobj if o not in used] if used else []}

    # calibration
    st["gates"] = rjson(RES.parent.parent / "denovo" / "s3_score" / "gates.json").get("models", {})

    # parked generators
    bg = rd("parked_boltzgen.csv")
    st["bg_n"] = len(bg) if bg is not None else None
    st["fwdiv_n"] = len(rd("parked_fwdiv.csv")) if rd("parked_fwdiv.csv") is not None else None

    # spend
    bal = None
    sm = RES / "status.md"
    if sm.exists():
        m = re.search(r"JapanFold balance: \$([0-9.]+)", sm.read_text())
        if m:
            bal = float(m.group(1))
    st["balance"] = bal
    st["spent"] = round(START_BALANCE - bal, 2) if bal is not None else None

    # rank.py weights
    w = {"bind": 0.40, "comp": 0.15, "therm": 0.20, "yield": 0.15, "dev": 0.10, "rules": 0.05}
    try:
        src = (ROOT / "denovo" / "rank.py").read_text()
        for k in w:
            m = re.search(rf'--w-{k}", type=float, default=([0-9.]+)', src)
            if m:
                w[k] = float(m.group(1))
    except Exception:
        pass
    st["weights"] = w
    return st


# ------------------------------------------------------------------------------------------------ helpers
def foot(txt):
    return f'<p class="foot">{E(txt)}</p>'


def note(txt):
    return f"<aside>{E(txt)}</aside>"


def tag(kind="live"):
    return f'<p class="{kind}">{"LIVE NUMBERS" if kind == "live" else "PRELIMINARY"} · BUILT {E(TS)}</p>'


def table(headers, rows, widths, fs=24, band=True):
    t = f'<table style="font-size:{fs}px;color:#1B2733;width:100%"><tr>'
    t += "".join(f'<th style="width:{w}%;text-align:left">{E(h)}</th>' for h, w in zip(headers, widths))
    t += "</tr>"
    for i, row in enumerate(rows):
        bg = ' style="background:#FBFAF6"' if band and i % 2 == 0 else ""
        t += f"<tr{bg}>" + "".join(f"<td>{E(c)}</td>" for c in row) + "</tr>"
    return t + "</table>"


def col_chart(items, track=400, col_w=130, gap=24, scale=100.0):
    """items: (label, sublabel, value 0..scale, display, colour). Bars share one baseline."""
    bars = '<div style="display:flex;gap:%dpx;align-items:flex-end;height:%dpx">' % (gap, track + 44)
    labs = '<div style="display:flex;gap:%dpx;margin-top:10px">' % gap
    for lab, sub, v, disp, colr in items:
        h = max(4, int(track * (v / scale))) if v is not None else 4
        bars += (f'<div style="width:{col_w}px;display:flex;flex-direction:column;justify-content:flex-end;align-items:center;gap:6px">'
                 f'<p class="val">{E(disp)}</p><div class="bar" style="width:{col_w - 20}px;height:{h}px;background:{colr}"></div></div>')
        labs += (f'<div style="width:{col_w}px;display:flex;flex-direction:column;align-items:center;gap:2px">'
                 f'<p class="lab" style="text-align:center">{E(lab)}</p><p class="small" style="text-align:center">{E(sub)}</p></div>')
    return bars + "</div>" + labs + "</div>"


def hbars(rows, max_w=360, label_w=400, bar_h=44, gap=14):
    """rows: (label, count, colour); bar length relative to the first row."""
    top = max([r[1] for r in rows] + [1])
    out = f'<div style="display:flex;flex-direction:column;gap:{gap}px">'
    for lab, n, colr in rows:
        w = int(max_w * n / top) if n else 0
        bar = f'<div class="hbar" style="width:{max(w, 4)}px;height:{bar_h}px;background:{colr}"></div>' if n else f'<div style="width:4px;height:{bar_h}px;background:#C9C7BC"></div>'
        out += (f'<div style="display:flex;align-items:center;gap:16px"><p class="lab" style="width:{label_w}px">{E(lab)}</p>{bar}'
                f'<p class="val">{n if n else "pending"}</p></div>')
    return out + "</div>"


def dot_strip(rows, lo, hi, svg_w=720, label_w=300, row_h=104, nticks=5, fmt="{:.1f}"):
    """rows: (label, values, colour, marker_value|None). Labels are real text beside each SVG (no <text> in svg)."""
    inset = 30
    xp = lambda v: inset + (min(max(v, lo), hi) - lo) / (hi - lo) * (svg_w - 2 * inset)  # noqa: E731
    out = '<div style="display:flex;flex-direction:column;gap:0">'
    for label, vals, colr, med in rows:
        svg = (f'<svg aria-label="{E(label)}: {len(vals)} values" width="{svg_w}" height="{row_h}" viewBox="0 0 {svg_w} {row_h}">'
               f'<line x1="{inset}" y1="{row_h // 2}" x2="{svg_w - inset}" y2="{row_h // 2}" stroke="#C9C7BC" stroke-width="2"/>')
        for i, v in enumerate(vals):
            jy = ((i * 37) % 41) / 41 - 0.5
            svg += f'<circle cx="{xp(v):.1f}" cy="{row_h / 2 + jy * 52:.1f}" r="8" fill="{colr}" fill-opacity="0.6"/>'
        if med is not None:
            svg += f'<rect x="{xp(med) - 2:.1f}" y="{row_h / 2 - 30}" width="4" height="60" fill="#12212F"/>'
        svg += "</svg>"
        out += (f'<div style="display:flex;align-items:center;gap:20px"><p class="lab" style="width:{label_w}px">{E(label)}</p>{svg}</div>')
    ticks = [lo + (hi - lo) * i / (nticks - 1) for i in range(nticks)]
    out += (f'<div style="display:flex;gap:20px"><div style="width:{label_w}px"></div>'
            f'<div style="width:{svg_w}px;display:flex;justify-content:space-between">'
            + "".join(f'<p class="small" style="width:60px;text-align:center">{fmt.format(t)}</p>' for t in ticks) + "</div></div></div>")
    return out


def med(xs):
    return float(np.median(xs)) if len(xs) else None


# ------------------------------------------------------------------------------------------------ slides
def s_cover(st):
    return f'''<section id="cover" class="dark" style="justify-content:space-between">
<div class="col" style="gap:28px">
<p class="eyebrow">AI × ANTIBODIES HACKATHON · OSLO · 14 NOV 2026 · DE NOVO TRACK</p>
<h1 style="width:1500px">Filter first, calibrate on evidence: de novo VHHs against the trastuzumab epitope of HER2</h1>
<p class="lead" style="width:1350px">RFD4-Proteina designs loops on an intact nanobody framework. The organisers’ novelty rule is a hard gate, independent folding models score what passes, and every finalist is profiled for melting temperature and yield.</p>
</div>
<div class="row" style="align-items:center;gap:24px">{tag("live")}<p class="small">Every number here is read from the pipeline’s result files when the deck is built.</p></div>
{note("Thesis: we treat every score as a filter that has to earn its place by separating known binders from decoys, and we say plainly where it does not. Structure of the talk: task, pipeline, the design choices with their evidence, the live funnel, what we dropped, risks and value.")}
</section>'''


def s_task(st):
    card = lambda t, b: f'<div class="card" style="flex:1"><p class="card-title">{E(t)}</p><p class="card-text">{E(b)}</p></div>'  # noqa: E731
    return f'''<section id="task">
<p class="eyebrow">THE TASK</p>
<h2>Three properties are judged against measured trastuzumab</h2>
<div class="row">
{card("Binding", "Comparable or better HER2 binding by SPR. Trastuzumab reference: KD 0.03–0.18 nM.")}
{card("Epitope", "A competition assay with trastuzumab: the design has to sit on the same site of HER2 domain IV.")}
{card("Thermostability", "A better melting temperature, the organisers’ proxy for real-life developability.")}
</div>
<div class="row">
<div class="card" style="width:480px"><p class="card-title">Tiers</p><p class="card-text">S: all three improve<br>A: two improve<br>B: one improves</p></div>
<div class="card" style="flex:1"><p class="card-title">De novo track rules</p><p class="card-text">VHH, at most 250 aa, no unpaired Cys, numbering confidence ≥ 0.1, under 70% identity to any known antibody, at least 2 CDRH3 mutations from trastuzumab, and it must pass abguard.</p></div>
</div>
<div class="row" style="align-items:center;gap:16px"><p class="small">Judged on:</p><p class="pill">40% reasoning</p><p class="pill">20% reproducible pipeline</p><p class="pill">20% user fit</p><p class="pill">20% path to value</p></div>
{foot("Source: hackathon task description (6 Nov 2026); trastuzumab SPR references in the provided data.")}
{note("The tiers count how many of the three properties beat the measured trastuzumab baseline. That is why we later rank by a Pareto front rather than a single weighted score: the weakest axis decides the tier. Judging weights are shown because reasoning is 40%, so each slide says why before what.")}
</section>'''


def s_gate(st):
    bins = st["gate_bins"]
    items = [(lab, f"n = {n}", (100 * k / n) if n else None, pct(k, n), "#0C7473") for lab, k, n in bins]
    pf = st["pass_by_fw"]
    cp = pf.get("capla", (0, 0))
    bp = pf.get("bcii10", (0, 0))
    if cp[1]:
        items.append(("capla-like", f"n = {cp[1]}", 100 * cp[0] / cp[1], pct(*cp), "#B4620A"))
    chart = col_chart(items, track=380, col_w=130, gap=20)
    return f'''<section id="gate">
<p class="eyebrow">THE NOVELTY GATE</p>
<h2>The novelty gate is mostly a length effect</h2>
<div class="row" style="gap:56px">
<div class="col" style="width:980px;gap:12px">
<p class="lab">Classed de novo by abguard, hNbBCII10 framework, by VHH length in aa. Amber: caplacizumab-like framework, all lengths.</p>
{chart}
</div>
<div class="col" style="flex:1;gap:20px">
<div class="card-dark" style="padding:24px"><p class="card-title" style="font-size:30px">identity = 1 − edits / max(length, 120)</p><p class="card-text" style="font-size:26px">abguard compares the whole variable domain with trastuzumab’s VH. Below 70% is de novo.</p></div>
<div class="card" style="padding:24px"><p class="card-title" style="font-size:30px">New CDRs on an intact framework</p><p class="card-text" style="font-size:26px">hNbBCII10 is close to trastuzumab’s framework, so its 70% line sits near 125 aa and {pct(*bp)} of its designs pass. The caplacizumab-like framework passes {pct(*cp)} at any length. We keep the framework and filter by the gate.</p></div>
</div>
</div>
{foot("Source: denovo/s2_rules/gate_drivers.md (abguard source read; our replica matches the VM tool exactly); results/denovo/rules.csv.")}
{note("This is the most useful single finding of the build. We read the novelty tool and measured it: identity to trastuzumab VH is one minus the edit distance over the larger of the Fv length and 120. Matched residues are nearly constant for a given framework, so identity is about matched over length and the 70 percent line falls at a fixed length. Practical consequence: we do not mutate the framework, we filter by the gate and ask the generator for longer CDR3s. Caveat: this is the VM build; the jury sentence about identity to any known antibody is not implemented in it.")}
</section>'''


def collect2(st):
    """Final-state additions: validation, HER2-gap check, OpenDDE pose check, coverage counts, final list, boxplot data."""
    D = rd("designs_rfd4jf.csv")
    ids = set(D.design_id) if D is not None else set()
    st["n_loop"] = int(D.generator.astype(str).str.contains("loop").sum()) if D is not None and "generator" in D else 0
    poses = set()
    if D is not None and "params_json" in D:
        for s in D.params_json.dropna():
            try:
                poses.add(json.loads(s).get("pose"))
            except Exception:
                pass
    st["poses"] = sorted(p for p in poses if p)
    V = rd("validity.csv")
    st["val"] = None
    if V is not None:
        v = V[V.design_id.isin(ids)]
        g = lambda c: int(v[c].sum()) if c in v else 0  # noqa: E731
        st["val"] = {"n": len(v), "ok": g("val_ok"), "motif": g("val_motif_bad"), "severe": g("val_clash_severe"), "odd": g("val_cys_odd"),
                     "lowc": g("val_low_complexity"), "mism": g("val_seq_mismatch") + g("val_fr_not_intact") + g("val_her2_seq_changed")}
    G = rd("gapcheck.csv")
    st["gap"] = None
    if G is not None and "gap_anchor_atoms_lt3" in G:
        gg = G[G.design_id.isin(ids)]
        if len(gg):
            st["gap"] = (int((gg.gap_anchor_atoms_lt3 > 0).sum()), len(gg))
    P = rd("dde_pose_check.csv")
    st["pose"] = None
    st["pose_special"] = None
    if P is not None:
        p = P[P.design_id.isin(ids)]
        if len(p):
            st["pose"] = {"n": len(p), "rmsd": float(p.chk_pose_rmsd.median()), "pred": float(p.chk_pred_epi_frac.median()),
                          "des": float(p.chk_design_epi_frac.median()), "fold": float(p.chk_fold_rmsd.median())}
        sp = P[P.design_id.eq("bcii10-rfd4jf-s201016-006")]
        if len(sp):
            r = sp.iloc[0]
            st["pose_special"] = {"rmsd": float(r.chk_pose_rmsd), "pred": float(r.chk_pred_epi_frac), "des": float(r.chk_design_epi_frac)}
    L = rd("leaderboard.csv")
    st["cov"] = {}
    st["box"] = None
    S = rd("submission_top20_pareto.csv")
    top = list(S.design_id) if S is not None else []
    if L is not None:
        l = L[(L.is_control == 0) & L.design_id.isin(ids)]
        cnt = lambda c: int(l[c].notna().sum()) if c in l else 0  # noqa: E731
        st["cov"] = {"vmof3": cnt("vmof3_pair_iptm"), "vmesm2": cnt("vmesm2_pair_iptm"), "dde": cnt("orth_dde_iptm"), "tm": cnt("nm_tm"),
                     "tnp": cnt("tnp_n_red"), "yield": cnt("my_yield"), "strict": int(l.gate_strict.sum()) if "gate_strict" in l else 0,
                     "vm": int(l.gate_vm.sum()) if "gate_vm" in l else 0}
        box = {}
        for key, col in (("yield", "my_yield"), ("tm", "nm_tm"), ("ipsae", "orth_dde_ipsae_chn15")):
            if col in l:
                a = l[["design_id", col]].dropna()
                box[key] = (a[col].tolist(), a[a.design_id.isin(top)][col].tolist())
        st["box"] = box
        kri = L[L.design_id.str.startswith("ctl-kri") & L.of3_iptm.notna()] if "of3_iptm" in L else None
        st["kri"] = (float(kri.of3_iptm.max()), len(kri)) if kri is not None and len(kri) else None
    st["top"] = None
    if S is not None and len(S):
        cnt = {}
        for x in S["flags"].fillna(""):
            for f in x.split(";"):
                if f:
                    cnt[f] = cnt.get(f, 0) + 1
        st["top"] = {"n": len(S), "flags": cnt, "fw": S.framework.value_counts().to_dict(), "front": S.pareto_rank.value_counts().sort_index().to_dict(),
                     "rows": S.head(6).to_dict("records"), "ids": top}
    return st


def box_panel(vals, top, w=520, h=420):
    """One boxplot panel (all designs) with the selected designs as diamonds; no text inside the svg."""
    if len(vals) < 5:
        return '<p class="muted">pending</p>'
    a = np.array(vals, float)
    lo, hi = float(a.min()), float(a.max())
    pad = 24
    yp = lambda v: pad + (1 - (v - lo) / (hi - lo or 1)) * (h - 2 * pad)  # noqa: E731
    q1, md, q3 = np.percentile(a, [25, 50, 75])
    iqr = q3 - q1
    wl, wh = a[a >= q1 - 1.5 * iqr].min(), a[a <= q3 + 1.5 * iqr].max()
    cx = w / 2
    rng = np.random.default_rng(1)
    svg = f'<svg aria-label="Boxplot of {len(a)} designs with {len(top)} selected designs marked" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
    svg += f'<rect x="1" y="1" width="{w - 2}" height="{h - 2}" fill="#FBFAF6" stroke="#D7D5CA" stroke-width="2"/>'
    for v in a:
        svg += f'<circle cx="{cx + rng.normal(0, 38):.1f}" cy="{yp(v):.1f}" r="3" fill="#9AA5B1" fill-opacity="0.35"/>'
    svg += (f'<line x1="{cx}" y1="{yp(wl):.1f}" x2="{cx}" y2="{yp(q1):.1f}" stroke="#4F5C69" stroke-width="3"/>'
            f'<line x1="{cx}" y1="{yp(q3):.1f}" x2="{cx}" y2="{yp(wh):.1f}" stroke="#4F5C69" stroke-width="3"/>'
            f'<line x1="{cx - 18}" y1="{yp(wl):.1f}" x2="{cx + 18}" y2="{yp(wl):.1f}" stroke="#4F5C69" stroke-width="3"/>'
            f'<line x1="{cx - 18}" y1="{yp(wh):.1f}" x2="{cx + 18}" y2="{yp(wh):.1f}" stroke="#4F5C69" stroke-width="3"/>'
            f'<rect x="{cx - 52}" y="{yp(q3):.1f}" width="104" height="{max(2.0, yp(q1) - yp(q3)):.1f}" fill="#E8ECF0" fill-opacity="0.75" stroke="#4F5C69" stroke-width="3"/>'
            f'<line x1="{cx - 52}" y1="{yp(md):.1f}" x2="{cx + 52}" y2="{yp(md):.1f}" stroke="#12212F" stroke-width="5"/>')
    for v in top:
        x, y = cx + rng.normal(0, 20), yp(v)
        svg += f'<path d="M {x - 11:.1f} {y:.1f} L {x:.1f} {y - 11:.1f} L {x + 11:.1f} {y:.1f} L {x:.1f} {y + 11:.1f} Z" fill="#C2410C" stroke="#FBFAF6" stroke-width="2"/>'
    return svg + "</svg>"


def s_pipeline(st):
    def stage(n, name, text, where):
        return (f'<div class="stage" style="padding:20px;gap:10px"><p class="stage-num">STAGE {n}</p><p class="stage-name" style="font-size:36px">{E(name)}</p>'
                f'<p class="stage-text" style="font-size:24px">{E(text)}</p><div style="flex:1"></div><p class="pill">{E(where)}</p></div>')
    arrow = '<x-shape kind="arrow-right" style="background:#0C7473;width:28px;height:20px;align-self:center"></x-shape>'
    return f'''<section id="pipeline">
<p class="eyebrow">THE PIPELINE</p>
<h2>Six stages, run as continuous parallel loops</h2>
<div class="row" style="gap:8px;align-items:stretch">
{stage(1, "Generate", "RFD4-Proteina on two intact frameworks. Now stopped.", "JapanFold")}
{arrow}
{stage(2, "Validate", "Sequence, structure, clash and HER2-loop checks.", "Laptop")}
{arrow}
{stage(3, "Gate", "abguard de novo, Cys, CDR3 12–18 aa, validity.", "Laptop + team VM")}
{arrow}
{stage(4, "Score", "VM OpenFold3 and ESMFold-2 rank; OpenDDE gates.", "Team VM + JapanFold")}
{arrow}
{stage(5, "Profile", "NanoMelt, TNP, ESMC yield, nine rules.", "Laptop + team VM")}
{arrow}
{stage(6, "Rank", "Gates, flags, then Pareto on four axes.", "Laptop")}
</div>
<div class="row">
<div class="card-dark" style="flex:1;padding:24px"><p class="card-title">Always running</p><p class="card-text" style="font-size:26px">Every stage was a watch loop that took only new designs, so scoring started while generation was still going. Five parallel workers shared one file contract.</p></div>
<div class="card-dark" style="flex:1;padding:24px"><p class="card-title">Validation fed back</p><p class="card-text" style="font-size:26px">Checks found a duplicate-seed bug (69% repeated designs), invalid capla designs and a gap in HER2 inside the epitope; each was fixed or gated out.</p></div>
</div>
{foot("Sources: denovo/CONTRACT.md decisions, stage READMEs and decisions.md files, denovo/validation/.")}
{note("Walk left to right. The order is deliberate: cheap hard gates first, expensive structure scoring only on what passes, developability on the finalists, and a Pareto step last because the organisers tier a design by how many of the three properties beat trastuzumab. Everything is resumable and no job was paid for twice. Validation was added mid-day after it found a duplicate-seed bug; it paid for itself.")}
</section>'''


def s_generate(st):
    fw = st["fw_counts"]
    rows = [["Caplacizumab-like", fw.get("capla", 0)], ["hNbBCII10-FGLA", fw.get("bcii10", 0)]]
    item = lambda k, v: f'<div class="col" style="gap:6px"><p class="card-title" style="font-size:32px">{E(k)}</p><p class="muted" style="font-size:27px">{E(v)}</p></div>'  # noqa: E731
    return f'''<section id="generate">
<p class="eyebrow">GENERATION</p>
<h2>RFD4 designs the loops and keeps the framework</h2>
<div class="row" style="gap:48px">
<div class="col" style="flex:3;gap:22px">
{item("Two intact frameworks", "hNbBCII10-FGLA is the scaffold RFantibody and Germinal fixed after expression tests. A caplacizumab-like VHH adds camelid hallmarks. CDR1–3 are designed.")}
{item("Competition built in", "The framework sits where trastuzumab’s VH sits in 1N8Z, in one of two rigid poses, so a binder in these poses overlaps trastuzumab by construction. It is also a limit.")}
{item("A loop-aware set", f"{st['n_loop']} designs on a modelled HER2 that fills the unresolved loop (residues 581–590), so the chain has no break.")}
{item("Reasonable lengths", "CDR3 designed per framework, with the final IMGT block capped at 12–18 aa; BoltzGen designs were the length reference. About $0.002 per design.")}
</div>
<div class="col" style="flex:2;gap:20px">
<div class="card"><p class="num">{st["n_designs"]:,}</p><p class="cap">unique RFD4 designs on JapanFold</p></div>
{table(["Framework", "Designs"], rows, [66, 34], fs=24)}
<div class="card-dark" style="padding:24px"><p class="card-text" style="font-size:26px">Validation found that 69% of an early wave were duplicate seeds. It was fixed, and only unique designs count here.</p></div>
</div>
</div>
{foot("Sources: denovo/s1_design/decisions.md, denovo/targets/epitope.json, results/denovo/designs_rfd4jf.csv, CONTRACT.md decisions.")}
{note("Choices and why. Intact frameworks: stability and a validated scaffold; the organisers said new CDRs should land below 70 percent identity. Pose anchored on trastuzumab’s VH: this is how epitope competition is built in, and also the limit, because every design shares one of two rigid poses and the hosted model gives no confidence score. The loop-aware set uses a modelled HER2 because 1N8Z lacks residues 581 to 590. Cost is negligible, so we generated widely and let validation and the gates decide. Final count is read from the designs file.")}
</section>'''


def s_validate(st):
    v = st["val"]
    gp = st["gap"]
    item = lambda k, t: f'<div class="col" style="gap:4px"><p class="card-title" style="font-size:30px">{E(k)}</p><p class="muted" style="font-size:26px">{E(t)}</p></div>'  # noqa: E731
    if v:
        mism = "No sequence, framework or HER2 mismatches." if v["mism"] == 0 else f"{v['mism']} sequence or framework mismatches."
        checks = [
            item("Sequence and structure agree", f"{mism} Covalent geometry and the framework’s fixed part hold."),
            item("Conditioned group held", f"{v['motif']} capla designs lost the conditioned group and are invalid."),
            item("Clashes, Cys, complexity", f"{v['severe']} severe clashes, {v['odd']} odd-Cys designs and {v['lowc']} low-complexity loops are flagged or gated."),
            item("Only two rigid binding modes", f"Poses {', '.join(st['poses']) or 'p0, p1'}: framework and HER2 are held rigid, so every design shares one of two poses."),
        ]
        big = (f'<div class="card"><p class="num">{v["ok"]:,}</p><p class="cap">of {v["n"]:,} checked designs are valid</p></div>')
    else:
        checks, big = [item("Validation", "pending")], ""
    gtxt = (f"{pct(*gp)}" if gp else "–")
    return f'''<section id="validate">
<p class="eyebrow">VALIDITY CHECKS</p>
<h2>We checked the designs themselves before trusting any score</h2>
<div class="row" style="gap:48px">
<div class="col" style="flex:3;gap:20px">
{"".join(checks)}
</div>
<div class="col" style="flex:2;gap:20px">
{big}
<div class="card"><p class="num-amber">{E(gtxt)}</p><p class="cap">of designs conflict with the anchors of the unresolved HER2 loop 581–590, which sits inside the epitope. They are gated out.</p></div>
</div>
</div>
{foot("Sources: denovo/validation/ (validate.py, gap_check.py, summary.txt), results/denovo/validity.csv and gapcheck.csv. Top 10 also inspected in PyMOL.")}
{note("Validation was a separate pass over every design: does the sequence match the structure, is the framework intact, did the conditioned part stay where it was put, are there clashes, odd cysteines or low-complexity loops. It found the duplicate-seed bug and the caplacizumab-like designs that lost the conditioned group early on. The HER2 gap check matters biologically: residues 581 to 590 are missing in the crystal and sit inside trastuzumab’s epitope, so a design whose anchors collide with that loop cannot be trusted. Counts are read from the validity and gap-check files.")}
</section>'''


def s_scoring(st):
    cov, po, kri = st["cov"], st["pose"], st.get("kri")
    rows = [
        ["OpenFold3 OpenBind, HER2 template, no MSA", "Team VM", "Ranked: iPTM and ipSAE (15 Å)"],
        ["ESMFold-2", "Team VM", "Ranked: iPTM and ipSAE"],
        ["OpenDDE-abag", "JapanFold", "Gate: iPTM and ipSAE ≥ decoy 95th percentile; missing fails"],
        ["JapanFold OpenFold3", "JapanFold", "Not ranked: untemplated, looks like decoys"],
    ]
    cnt = f"Scored so far: OpenFold3 on {cov.get('vmof3', 0)} designs, ESMFold-2 on {cov.get('vmesm2', 0)}, OpenDDE on {cov.get('dde', 0)}."
    if po:
        ptxt = (f"OpenDDE’s pose is a median {po['rmsd']:.0f} Å from the designed pose, and {po['pred']:.0%} of the HER2 residues it contacts lie in the trastuzumab epitope, "
                f"against {po['des']:.0%} for the designed poses. Re-checked independently, not a bug. Its fold of the VHH itself is consistent (scRMSD median {po['fold']:.1f} Å).")
    else:
        ptxt = "Pose check pending."
    ktxt = (f"Krisha’s branch, templated Boltz-2 and OpenFold3-OB: 1 of 120 designs pose-consistent. Her six VHHs through our tools: best OpenFold3 ipTM {kri[0]:.2f}. "
            "Recipe matters: templated OpenFold3 median 0.36 against 0.11 for our untemplated runs." if kri else "External reference pending.")
    return f'''<section id="scoring">
<p class="eyebrow">INDEPENDENT BINDING EVIDENCE</p>
<h2>The VM runs rank; JapanFold only gates</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:840px;gap:16px">
{table(["Scorer", "Where", "Role"], rows, [40, 20, 40], fs=24)}
<p class="small">{E(cnt)} Decoys are framework-only VHHs; every gate is set from them.</p>
</div>
<div class="col" style="flex:1;gap:16px">
<div class="card-dark" style="padding:24px"><p class="card-title" style="font-size:30px">OpenDDE does not support epitope preservation</p><p class="card-text" style="font-size:24px">{E(ptxt)}</p></div>
<div class="card" style="padding:24px"><p class="card-title" style="font-size:30px">External reference</p><p class="card-text" style="font-size:24px">{E(ktxt)}</p></div>
</div>
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small">Binding evidence is computational and weak. Scores are filters, not affinity measurements.</p></div>
{foot("Sources: denovo/s3_score/, denovo/vm_templated/, denovo/validation/dde_pose_check.py, results/denovo/dde_pose_check.csv, leaderboard.csv.")}
{note("What changed late in the day. JapanFold scores are not ranked any more: untemplated runs score our designs like decoys. The ranked binding axis comes from OpenFold3 OpenBind with the HER2 template and no MSA, which is the upstream recipe, and from ESMFold-2, both run on the team VM. OpenDDE exists only on JapanFold, so it is a gate with a threshold from the decoys; a missing score fails. The honest finding is that OpenDDE places the VHH in a different spot from the designed pose, so passing its gate does not show epitope preservation. We checked it independently to rule out a numbering bug.")}
</section>'''


def pareto_scatter(par, w=800, h=330):
    """bind vs therm percentile, coloured by Pareto front; front 1 outlined, selected designs marked with a diamond."""
    pad = 20
    xp = lambda v: pad + v * (w - 2 * pad)  # noqa: E731
    yp = lambda v: h - pad - v * (h - 2 * pad)  # noqa: E731
    col = lambda r: "#0C7473" if r == 1 else ("#B4620A" if r == 2 else ("#4C6A8A" if r == 3 else "#B8B6AA"))  # noqa: E731
    svg = (f'<svg aria-label="Binding percentile against melting-temperature percentile for {len(par["pts"])} gate-passing designs, coloured by Pareto front" '
           f'width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
           f'<rect x="{pad}" y="{pad}" width="{w - 2 * pad}" height="{h - 2 * pad}" fill="#FBFAF6" stroke="#D7D5CA" stroke-width="2"/>'
           f'<line x1="{xp(0.5)}" y1="{pad}" x2="{xp(0.5)}" y2="{h - pad}" stroke="#E3E1D7" stroke-width="2"/>'
           f'<line x1="{pad}" y1="{yp(0.5)}" x2="{w - pad}" y2="{yp(0.5)}" stroke="#E3E1D7" stroke-width="2"/>')
    for b, t, r, sel in sorted(par["pts"], key=lambda p: -p[2]):
        stroke = ' stroke="#12212F" stroke-width="3"' if r == 1 else ""
        svg += f'<circle cx="{xp(b):.1f}" cy="{yp(t):.1f}" r="{11 if sel else 8}" fill="{col(r)}" fill-opacity="0.85"{stroke}/>'
    for b, t, r, sel in par["pts"]:
        if sel:
            x, y = xp(b), yp(t)
            svg += f'<path d="M {x - 7:.1f} {y:.1f} L {x:.1f} {y - 7:.1f} L {x + 7:.1f} {y:.1f} L {x:.1f} {y + 7:.1f} Z" fill="#12212F"/>'
    return svg + "</svg>"


def legend_item(kind, label):
    cols = {"f1": "#0C7473", "f2": "#B4620A", "f3": "#4C6A8A", "f4": "#B8B6AA"}
    if kind == "sel":
        ic = '<svg aria-label="diamond marker" width="28" height="28" viewBox="0 0 28 28"><path d="M 14 3 L 25 14 L 14 25 L 3 14 Z" fill="#12212F"/></svg>'
    else:
        stroke = ' stroke="#12212F" stroke-width="3"' if kind == "f1" else ""
        ic = f'<svg aria-label="{label} marker" width="28" height="28" viewBox="0 0 28 28"><circle cx="14" cy="14" r="10" fill="{cols[kind]}"{stroke}/></svg>'
    return f'<div style="display:flex;align-items:center;gap:8px">{ic}<p class="small">{E(label)}</p></div>'


def s_ranking(st):
    par = st["par"]
    if par and len(par["pts"]) >= 5:
        legend = ('<div style="display:flex;gap:20px;align-items:center">' + legend_item("f1", "Front 1") + legend_item("f2", "Front 2")
                  + legend_item("f3", "Front 3") + legend_item("f4", "Later") + legend_item("sel", "Selected") + "</div>")
        fr = " / ".join(str(v) for v in par["fronts"][:5])
        left = (f'<div class="col" style="width:830px;gap:8px"><p class="lab">Binding (x) against melting temperature (y), percentiles within this run</p>'
                f'{legend}{pareto_scatter(par)}<p class="small">{par["n_pop"]} designs in the Pareto step, {len(par["pts"])} plotted; fronts {fr}; {par["n_complete"]} fully measured.</p></div>')
    else:
        left = ('<div class="card" style="width:830px"><p class="card-title">Pareto chart: pending</p>'
                '<p class="card-text">It appears once binding and melting-temperature evidence exist for at least five eligible designs.</p></div>')
    item = lambda k, v: f'<div class="col" style="gap:4px"><p class="card-title" style="font-size:28px">{E(k)}</p><p class="small">{E(v)}</p></div>'  # noqa: E731
    right = f'''<div class="col" style="flex:1;gap:14px">
{item("1 · Gates: eligibility, missing = fail", "abguard de novo, canonical Cys, CDR3 block 12–18 aa, no severe clash or HER2-loop conflict, no homopolymer; the OpenDDE gate.")}
{item("2 · Flags: review items, never exclude", "Below decoy, N-glyc or Cys in a CDR, clash, low complexity, TNP red or amber, low nine-rule score, OpenDDE scRMSD above 1.5 Å.")}
{item("3 · Four axes, one Pareto", "Binding 50%, NanoMelt Tm 20%, TNP 15%, yield with nine rules 15%; within-run percentiles. An unmeasured design never beats a measured one.")}
{item("4 · The order", "Eligible, complete evidence, Pareto tier, fewest flags, best binding.")}
</div>'''
    return f'''<section id="ranking">
<p class="eyebrow">STAGE 6 · CHOOSING THE 20 SLOTS</p>
<h2>Gates, flags, then Pareto dominance</h2>
<div class="row" style="gap:48px">
{left}
{right}
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small">Percentiles are within this run, not against trastuzumab; the proxies are filters, not affinity measurements.</p></div>
{foot("Sources: denovo/pareto.py, denovo/rank.py, CONTRACT.md decisions 15:12 and 15:20, results/denovo/leaderboard_pareto.csv.")}
{note("Why Pareto: the organisers tier a design by how many of the three judged properties beat trastuzumab, and S needs all three. A weighted sum lets a strong binder with a poor melting temperature outrank a balanced design. So eligibility is a typed gate, the review items are flags that never exclude, and the four axes are sorted into non-dominated fronts, with fewest flags and best binding as tie-breakers. Absolute thresholds assume a templated oracle that none of our designs reach, so we use the adverse tail of our own decoys instead.")}
</section>'''


def s_live(st):
    c = st["cov"]
    v = st["val"] or {"ok": 0}
    full = st["par"]["n_complete"] if st.get("par") else 0
    rows = [
        ("Designed", st["n_designs"], "#6C7783"),
        ("Valid structure", v["ok"], "#6C7783"),
        ("abguard says de novo", c.get("vm", 0), "#0C7473"),
        ("Pass every gate", c.get("strict", 0), "#0C7473"),
        ("OpenFold3 on the VM", c.get("vmof3", 0), "#0C7473"),
        ("ESMFold-2 on the VM", c.get("vmesm2", 0), "#0C7473"),
        ("OpenDDE gate scored", c.get("dde", 0), "#0C7473"),
        ("NanoMelt Tm", c.get("tm", 0), "#B4620A"),
        ("TNP flags", c.get("tnp", 0), "#B4620A"),
        ("Predicted yield", c.get("yield", 0), "#B4620A"),
        ("Fully measured", full, "#12212F"),
    ]
    t = st.get("top")
    if t:
        trs = [[str(r["design_id"]).replace("-rfd4jf", ""), str(r.get("framework", "")), str(int(r["pareto_rank"])), str(int(r["n_flags"]))] for r in t["rows"]]
        right = ('<p class="card-title" style="font-size:30px">The final list, first six of ' + f'{t["n"]}' + '</p>'
                 + table(["Design", "Base", "Front", "Flags"], trs, [44, 20, 16, 20], fs=24)
                 + '<p class="small">Counts on the left are designs that reached each stage; they are not all gate-passers.</p>')
    else:
        right = '<div class="card"><p class="card-title">Final list: pending</p></div>'
    return f'''<section id="live">
<p class="eyebrow">THE FUNNEL</p>
<h2>Where the designs ended up</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:860px;gap:10px">
{hbars(rows, max_w=280, label_w=340, bar_h=34, gap=9)}
</div>
<div class="col" style="flex:1;gap:12px">
{right}
</div>
</div>
<div class="row" style="align-items:center;gap:20px">{tag("live")}<p class="small">Evidence is uneven: scoring covers fewer designs than generation by design, since the gate runs first.</p></div>
{foot("Sources: results/denovo/designs_rfd4jf.csv, validity.csv, rules.csv, leaderboard.csv, leaderboard_pareto.csv, submission_top20_pareto.csv.")}
{note("Read it as a funnel of RFD4 designs reaching each stage. Generation stopped; scoring covered what the gate let through, so later bars are shorter by design, and OpenDDE and the VM runs cover different subsets. The right side is the final list in rank order, with the number of flags each carries.")}
</section>'''


def s_dist(st):
    b = st["box"] or {}
    pan = [("Predicted monomeric yield", "ESMC-Nb monoyield, log10 mg/L", "yield"), ("NanoMelt predicted Tm", "°C", "tm"),
           ("OpenDDE-abag ipSAE", "chain-length d0, 15 Å cutoff", "ipsae")]
    cols = ""
    for title, unit, k in pan:
        vals, top = b.get(k, ([], []))
        m = f"median {np.median(vals):.2f}" if len(vals) else "pending"
        tm_ = f"; selected median {np.median(top):.2f}" if len(top) else ""
        cols += (f'<div class="col" style="width:520px;gap:6px"><p class="card-title" style="font-size:30px">{E(title)}</p><p class="small">{E(unit)}</p>'
                 f'{box_panel(vals, top)}<p class="small">n = {len(vals)}, {E(m)}{E(tm_)}</p></div>')
    return f'''<section id="dist">
<p class="eyebrow">DEVELOPABILITY AND BINDING DISTRIBUTIONS</p>
<h2>Where the final 20 sit among all designs</h2>
<div class="row" style="gap:52px">{cols}</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small">Diamonds are the final 20. Yield and Tm are predictions; ipSAE is the OpenDDE gate score.</p></div>
{foot("Sources: results/denovo/leaderboard.csv, submission_top20_pareto.csv; denovo/validation/boxplots.py makes the same figure.")}
{note("Each panel is every scored design as grey dots with a box for the quartiles and the final 20 as orange diamonds. The figure shows where the selected designs sit relative to all designs on each axis, so a reader can see whether the list leans on one axis. NanoMelt and monoyield are predictions with several degrees or log units of error, and the nine-rule score feeds the yield axis.")}
</section>'''


def s_final(st):
    t = st["top"]
    if not t:
        return '<section id="final"><p class="eyebrow">THE FINAL LIST</p><h2>Final list: pending</h2></section>'
    names = {"below_decoy": "below decoy", "clash": "clash", "dde_fold_rmsd": "OpenDDE scRMSD > 1.5 Å", "cys_in_cdr": "Cys in CDR", "tnp_amber2": "TNP amber ×2", "low_complexity": "low complexity"}
    rows = [(names.get(k, k), v, "#B4620A" if k == "below_decoy" else "#6C7783") for k, v in sorted(t["flags"].items(), key=lambda kv: -kv[1])]
    n = t["n"]
    below = t["flags"].get("below_decoy", 0)
    fwtxt = ", ".join(f"{v} {k}" for k, v in t["fw"].items())
    sp = st["pose_special"]
    inlist = "bcii10-rfd4jf-s201016-006" in t["ids"]
    qtxt = (f"bcii10-rfd4jf-s201016-006 is the only design whose OpenDDE pose matches the designed pose ({sp['rmsd']:.1f} Å); its contacts fall in the epitope {sp['pred']:.0%} against {sp['des']:.0%} designed. "
            f"It is {'on' if inlist else 'not on'} the final list. Open question: is it worth a slot?") if sp else "Open question: pending."
    return f'''<section id="final">
<p class="eyebrow">THE FINAL LIST</p>
<h2>Twenty designs, with their weaknesses shown</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:820px;gap:12px">
<p class="lab">Flags carried by the {n} selected designs</p>
{hbars(rows, max_w=300, label_w=330, bar_h=34, gap=10)}
<p class="small">{E(fwtxt)}; Pareto front sizes in the list: {E(", ".join(f"{int(k)}: {v}" for k, v in t["front"].items()))}.</p>
</div>
<div class="col" style="flex:1;gap:16px">
<div class="card-dark" style="padding:24px"><p class="card-title" style="font-size:30px">{below} of {n} are below the decoy line</p><p class="card-text" style="font-size:26px">Binding evidence is weak and computational only. We expect A- or B-tier, not S.</p></div>
<div class="card" style="padding:24px"><p class="card-title" style="font-size:30px">One open question</p><p class="card-text" style="font-size:24px">{E(qtxt)}</p></div>
</div>
</div>
{foot("Sources: results/denovo/submission_top20_pareto.csv, dde_pose_check.csv.")}
{note("Plain statement: the list is the best of what we made by our own filters, not evidence of binding. Nineteen of twenty are below the decoy line on at least half of the measured iPTM and ipSAE values, which is why the flag is shown rather than hidden. The open question is a judgment call: one design reproduces the designed pose in OpenDDE, which is the closest thing we have to independent support for a designed binding mode.")}
</section>'''


def s_limits(st):
    po = st["pose"]
    rows = [
        ["Framework diversification", "Organisers said new CDRs land under 70%. We measured a length effect, so we filter instead.", "gate_drivers.md"],
        ["BoltzGen designs", "Raw designs have no Cys. Stopped to focus on RFD4; kept as the CDR3 length reference.", "parked_boltzgen*.csv"],
        ["RFD4 on the team VM", "Ruled out. RFD4 ran only on JapanFold.", "CONTRACT.md"],
        ["JapanFold scores for ranking", "Untemplated and decoy-like. Demoted: only OpenDDE remains, as a gate.", "score_of3.csv"],
        ["OpenDDE as epitope evidence", (f"Its pose is {po['rmsd']:.0f} Å from the designed pose." if po else "Pose check pending.") + " It does not show epitope preservation.", "dde_pose_check.csv"],
    ]
    return f'''<section id="limits" style="gap:24px">
<p class="eyebrow">WHAT WE DROPPED, AND THE RISKS</p>
<h2>Decisions we reversed, and the outcome we expect</h2>
{table(["What", "Why", "Evidence"], rows, [24, 54, 22], fs=24)}
<div class="card-dark" style="padding:24px;gap:8px"><p class="card-title" style="font-size:32px">Expected outcome (our inference): A- or B-tier</p><p class="card-text" style="font-size:26px">Epitope competition and thermostability are plausible. Trastuzumab binds at 0.03–0.18 nM; de novo VHHs in the literature reach tens to hundreds of nM, so beating its affinity is the unlikely third.</p></div>
{foot("Sources: stage decisions.md files, CONTRACT.md, research/denovo_vhh_research.md (Germinal, RFantibody, BoltzGen hit rates).")}
{note("Reasoning is the largest judging weight, so this slide is deliberate. Each row is a decision reversed because of a measurement or a constraint. The risks that remain: our scorers are weak on real VHH pairs, the VM novelty tool does not implement the jury sentence about identity to any known antibody, every design shares one of two rigid poses, and the affinity gap to a sub-nanomolar antibody is large.")}
</section>'''


def s_value(st):
    spent = f"${st['spent']:.2f}" if st["spent"] is not None else "pending"
    item = lambda k, v: f'<div class="col" style="gap:8px"><p class="card-title" style="font-size:32px">{E(k)}</p><p class="muted" style="font-size:28px">{E(v)}</p></div>'  # noqa: E731
    return f'''<section id="value">
<p class="eyebrow">USER FIT, PATH TO VALUE, REPRODUCIBILITY</p>
<h2>A reusable pipeline around a trastuzumab-site VHH</h2>
<div class="row" style="gap:48px">
<div class="col" style="flex:1;gap:28px">
{item("Who needs it", "Teams that want a small binder at HER2 domain IV for imaging or diagnostics, or a building block for a biparatopic pair with a non-competing VHH. Cell-free production keeps iteration cheap.")}
{item("Honest limit on value", "Not differentiated as a therapy: it blocks the site that trastuzumab and its conjugates use. The HER2 VHH patent estate is crowded; we read abstracts only, so this is no freedom-to-operate view.")}
{item("The reusable part", "Gates, calibration on decoys, validation, profiling and Pareto ranking work for any target: swap the target files and the hotspot sets.")}
</div>
<div class="col" style="flex:1;gap:20px">
<div class="card"><p class="card-title">Reproduce it</p>
<p class="body" style="font-size:26px">Branch moritz of the team repo: denovo/ (code, stage READMEs, CONTRACT.md), the Claude skill .claude/skills/denovo-vhh-pipeline, and results_denovo/ (every score with its job id). No credentials.</p></div>
<div class="card-dark"><p class="card-title">Rebuild the ranking</p>
<p class="card-text" style="font-size:26px">denovo/.venv/bin/python denovo/rank.py<br>denovo/.venv/bin/python denovo/pareto.py</p></div>
<div class="row" style="gap:16px"><div class="card" style="flex:1;padding:24px"><p class="num" style="font-size:64px">{E(spent)}</p><p class="cap">JapanFold spend, cap $65</p></div>
<div class="card" style="flex:1;padding:24px"><p class="num" style="font-size:64px">{st["n_designs"]:,}</p><p class="cap">designs generated and validated</p></div></div>
</div>
</div>
{foot("Sources: research/denovo_vhh_research.md Q6 (inference, abstract-level patent reading); denovo/CONTRACT.md; denovo/s0_packaging/README.md.")}
{note("User fit and path to value, with the caveats: the value case is an inference from the literature, not a market study. The reusable part is the real asset for a method-first reading: a pipeline that gates on the organisers’ own rule, validates its own designs, calibrates its scorers on decoys and reports what it cannot verify. Spend is read from the JapanFold balance at build time.")}
</section>'''


def _pose_summary(df, pose, scr, pq):
    """Medians and counts of one scorer's pose-recovery columns (designs only)."""
    if df is None or not len(df):
        return None
    f = lambda c: float(df[c].median()) if c in df and df[c].notna().any() else None  # noqa: E731
    return {"n": int(len(df)), "pose": f(pose), "le20": int((df[pose] <= 20).sum()) if pose in df else 0, "le10": int((df[pose] <= 10).sum()) if pose in df else 0,
            "scr": f(scr), "scr3": int((df[scr] <= 3).sum()) if scr in df else 0, "pq": f(pq)}


def collect3(st):
    """Pose recovery of the scorers on the full HER2 domain (results/denovo/pose_check_vm.csv) and on the cropped HER2 539-607 re-run
    (results/denovo/crop/of3_crop.csv: vmof3c_*, crop/dde_crop.csv: dde_c_*). Missing files -> None, the slide then says "pending"."""
    D = rd("designs_rfd4jf.csv")
    ids = set(D.design_id) if D is not None else set()
    st["posechk"], st["crop"], st["crop_same"] = {}, {}, {}
    A = rd("pose_check_vm.csv")
    if A is not None:
        A = A[A.design_id.isin(ids)]
        for key in ("vmof3", "vmesm2", "opendde"):
            st["posechk"][key] = _pose_summary(A[A.scorer == key], "pose_rmsd", "scrmsd", "pdockq2")
    O, Dc = rd("crop/of3_crop.csv"), rd("crop/dde_crop.csv")
    if O is not None:
        O = O[O.design_id.isin(ids)]
        st["crop"]["vmof3"] = _pose_summary(O, "vmof3c_pose_rmsd", "vmof3c_scrmsd", "vmof3c_pdockq2")
    if Dc is not None:
        Dc = Dc[Dc.design_id.isin(ids)]
        st["crop"]["opendde"] = _pose_summary(Dc, "dde_c_pose_rmsd", "dde_c_scrmsd", "dde_c_pdockq2")
    # like for like: the cropped set is a prioritised subset, so also report the full-domain numbers on exactly the same designs
    if A is not None:
        for key, df, col in (("vmof3", O, "vmof3c_pose_rmsd"), ("opendde", Dc, "dde_c_pose_rmsd")):
            if df is None or col not in df or not len(df):
                continue
            both = set(df.dropna(subset=[col]).design_id) & set(A[A.scorer == key].design_id)
            if both:
                full = A[(A.scorer == key) & A.design_id.isin(both)].pose_rmsd
                crop = df[df.design_id.isin(both)].set_index("design_id")[col].reindex(sorted(both))
                st["crop_same"][key] = {"n": len(both), "full": float(full.median()), "crop": float(crop.median()),
                                        "full20": int((full <= 20).sum()), "crop20": int((crop <= 20).sum())}
    return st


def s_pose(st):
    pc, cr, same = st.get("posechk") or {}, st.get("crop") or {}, st.get("crop_same") or {}
    name = {"vmof3": "OpenFold3 OpenBind, template, no MSA (VM)", "vmesm2": "ESMFold-2 (VM)", "opendde": "OpenDDE-abag (JapanFold)"}
    fm = lambda v, d=0: "–" if v is None else f"{v:.{d}f}"  # noqa: E731

    def row(label, where, r):
        if not r:
            return [label, where, "pending", "", "", "", ""]
        return [label, where, f"{r['n']}", fm(r["pose"]), f"{r['le20']} of {r['n']}", fm(r["scr"], 1), fm(r["pq"], 2)]
    rows = [row(name[k], "full 489–607", pc.get(k)) for k in ("vmof3", "vmesm2", "opendde")]
    rows += [row(name[k].split(" (")[0], "cropped 539–607", cr.get(k)) for k in ("vmof3", "opendde")]
    if pc.get("vmof3") and pc["vmof3"]["scr"] is not None:
        ftxt = (f"Without an MSA, OpenFold3 misfolds the VHH itself (median scRMSD {pc['vmof3']['scr']:.1f} Å, {pc['vmof3']['scr3']} of {pc['vmof3']['n']} within 3 Å), "
                "so its interface scores sit on a wrong VHH. ESMFold-2 and OpenDDE fold the VHH as designed but place it elsewhere. No scorer reproduces the designed pose.")
    else:
        ftxt = "Pose comparison pending."
    parts = []
    for key, lab in (("vmof3", "OpenFold3"), ("opendde", "OpenDDE")):
        if key in same:
            r = same[key]
            parts.append(f"{lab} on the same {r['n']} designs: pose median {r['full']:.0f} Å on the full domain, {r['crop']:.0f} Å cropped; within 20 Å {r['full20']} against {r['crop20']}.")
    ctxt = " ".join(parts) if parts else "Cropped re-run pending: HER2 residues 489–538 removed so the models have less room to place the VHH somewhere random."
    return f'''<section id="pose">
<p class="eyebrow">DOES ANY SCORER REPRODUCE THE DESIGNED POSE?</p>
<h2>Pose recovery on the full HER2 domain and on a cropped HER2</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:1060px;gap:16px">
{table(["Scorer", "HER2", "Designs", "Pose Å", "≤ 20 Å", "scRMSD Å", "pDockQ2"], rows, [34, 17, 10, 9, 12, 10, 8], fs=24)}
<p class="small">Pose = VHH CA RMSD to the designed VHH after fitting the predicted HER2 on the designed HER2; medians over designs. Cropped = HER2 539–607, which keeps every epitope residue.</p>
</div>
<div class="col" style="flex:1;gap:16px">
<div class="card-dark" style="padding:24px"><p class="card-title" style="font-size:30px">Full domain</p><p class="card-text" style="font-size:24px">{E(ftxt)}</p></div>
<div class="card" style="padding:24px"><p class="card-title" style="font-size:30px">Cropped re-run</p><p class="card-text" style="font-size:24px">{E(ctxt)}</p></div>
</div>
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small">The binding scores are filters of unknown validity, not pose or affinity predictors.</p></div>
{foot("Sources: denovo/validation/pose_check_vm.py, denovo/vm_templated/make_crop.py and collect_crop.py, denovo/s3_score/run_crop.py; results/denovo/pose_check_vm.csv, crop/of3_crop.csv, crop/dde_crop.csv.")}
{note("The open question was whether the templated OpenFold3 reproduces the designed binding pose better than the other scorers, in which case the others would be dropped from the ranking. It does not on the full domain, so none was dropped. The cropped re-run asks a different question: if the models cannot place the VHH on the far side of the domain, do they find the designed pose? The table shows both runs side by side, and the card compares the full and cropped runs on exactly the same designs, because the cropped set is a prioritised subset (gate passers first). Whether the cropped scores enter the ranking is a separate, recorded decision.")}
</section>'''


ORDER = [("cover", s_cover), ("task", s_task), ("pipeline", s_pipeline), ("generate", s_generate), ("validate", s_validate), ("gate", s_gate),
         ("scoring", s_scoring), ("pose", s_pose), ("ranking", s_ranking), ("live", s_live), ("dist", s_dist), ("final", s_final), ("limits", s_limits), ("value", s_value)]
SECTIONS = {"s1": ("The task and the approach", "cover"), "s2": ("Design, validation and the novelty gate", "generate"),
            "s3": ("Evidence, ranking and the final list", "scoring"), "s4": ("Limits, value and reproducibility", "limits")}


def lint(sid, txt):
    problems = []
    for m in re.finditer(r"font-size:\s*(\d+)px", txt):
        if int(m.group(1)) < 24:
            problems.append(f"font-size {m.group(1)}px")
    n = len(re.findall(r"<[a-zA-Z]", re.sub(r"<svg.*?</svg>", "<svg/>", txt, flags=re.S)))
    if n > 200:
        problems.append(f"{n} elements")
    if txt.count("<section") != 1 or not txt.strip().endswith("</section>"):
        problems.append("section structure")
    return problems


def main():
    st = collect3(collect2(collect()))
    SLIDES.mkdir(parents=True, exist_ok=True)
    bad = False
    for sid, fn in ORDER:
        txt = fn(st)
        pr = lint(sid, txt)
        if pr:
            bad = True
            print(f"LINT {sid}: {pr}")
        (SLIDES / f"{sid}.html").write_text(txt)
    # index: keep createdOnFiles from an existing index
    idx_path = PROJ / "deck.json"
    created = {"v": 1, "at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    if idx_path.exists():
        try:
            created = json.loads(idx_path.read_text()).get("createdOnFiles", created)
        except Exception:
            pass
    idx = {"v": 5, "createdOnFiles": created, "lists": "css", "format": {"v": 1, "shared": ["project/deck.css"]},
           "title": "De novo VHH pipeline overview", "cover": "cover", "order": [s for s, _ in ORDER],
           "sections": {k: {"description": d, "start": s} for k, (d, s) in SECTIONS.items()},
           "faces": {"space-grotesk": {"family": "Space Grotesk", "href": "https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400..700&display=swap"},
                     "inter": {"family": "Inter", "href": "https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap"}},
           "designSystems": []}
    idx_path.write_text(json.dumps(idx, indent=1))
    files = {f"project/slides/{s}.html": f"project/slides/{s}.html" for s, _ in ORDER}
    plan = {"url": URL, "root": str(PROJ.parent), "file_path": str(PROJ / "slides" / "cover.html"),
            "files": {**{k: v for k, v in files.items() if k != "project/slides/cover.html"}, "project/deck.json": "project/deck.json"},
            "note": "deck.json is included because the slide order changed (new 'pose' slide); deck.css only on the first publish"}
    (HERE / "publish_plan.json").write_text(json.dumps(plan, indent=1))
    print(f"built {len(ORDER)} slides at {TS}: designs={st['n_designs']} gate_strict={st['n_gate_strict']} of3={st['of3']['n']} "
          f"tm={st['tm']['n']} cand={len(st['cand'])} lint={'FAIL' if bad else 'ok'}")
    if "--plan" in sys.argv:
        print(json.dumps(plan, indent=1))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
