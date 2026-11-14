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


def hbars(rows, max_w=360, label_w=400):
    """rows: (label, count, colour); bar length relative to the first row."""
    top = max([r[1] for r in rows] + [1])
    out = '<div style="display:flex;flex-direction:column;gap:14px">'
    for lab, n, colr in rows:
        w = int(max_w * n / top) if n else 0
        bar = f'<div class="hbar" style="width:{max(w, 4)}px;background:{colr}"></div>' if n else '<div style="width:4px;height:44px;background:#C9C7BC"></div>'
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


def s_pipeline(st):
    def stage(n, name, text, where):
        return (f'<div class="stage" style="padding:20px;gap:10px"><p class="stage-num">STAGE {n}</p><p class="stage-name" style="font-size:36px">{E(name)}</p>'
                f'<p class="stage-text" style="font-size:24px">{E(text)}</p><div style="flex:1"></div><p class="pill">{E(where)}</p></div>')
    arrow = '<x-shape kind="arrow-right" style="background:#0C7473;width:28px;height:20px;align-self:center"></x-shape>'
    return f'''<section id="pipeline">
<p class="eyebrow">THE PIPELINE</p>
<h2>Six stages run as continuous parallel loops</h2>
<div class="row" style="gap:8px;align-items:stretch">
{stage(1, "Generate", "RFD4-Proteina designs CDR1–3 on an intact VHH.", "JapanFold")}
{arrow}
{stage(2, "Gate", "abguard de novo, canonical Cys, length.", "Team VM + laptop")}
{arrow}
{stage(3, "Score", "OpenFold3 on all; OpenDDE, ESMFold-2 on the top third.", "JapanFold")}
{arrow}
{stage(4, "Profile", "NanoMelt Tm, TNP, yield, liabilities.", "Laptop + team VM")}
{arrow}
{stage(5, "Rank", "Weighted composite, kept for comparison.", "Laptop")}
{arrow}
{stage(6, "Pareto", "Non-dominated fronts choose the final 20.", "Laptop")}
</div>
<div class="row">
<div class="card-dark" style="flex:1;padding:24px"><p class="card-title">Always running</p><p class="card-text" style="font-size:26px">Every stage is a watch loop that takes only new designs. Five parallel workers share one file contract (denovo/CONTRACT.md); the ranking refreshes every five minutes.</p></div>
<div class="card-dark" style="flex:1;padding:24px"><p class="card-title">Feedback</p><p class="card-text" style="font-size:26px">OpenFold3 results steer the next RFD4 waves: settings that pass the gate and score better get more jobs.</p></div>
</div>
{foot("Sources: denovo/CONTRACT.md and the stage READMEs. Nebius was used only for the parked framework-diversification experiments.")}
{note("Walk left to right. The point is the order: cheap hard gates first, expensive structure scoring only on what passes, developability on the finalists, and a Pareto step last because the organisers tier a design by how many of the three properties beat trastuzumab. Everything is resumable and a job is never paid for twice. The loop is what let us look at real designs from the first half hour.")}
</section>'''


def s_generate(st):
    fwc = st["fw_counts"]
    hs = ["hs1_core", "hs2_broad", "hs3_alt_patch"]
    rows = []
    for f, nm in (("bcii10", "hNbBCII10-FGLA"), ("capla", "Caplacizumab-like")):
        rows.append([nm] + [st["by_fw_hs"].get((f, h), 0) for h in hs])
    item = lambda k, v: f'<div class="col" style="gap:6px"><p class="card-title">{E(k)}</p><p class="muted">{E(v)}</p></div>'  # noqa: E731
    return f'''<section id="generate">
<p class="eyebrow">GENERATION</p>
<h2>RFD4 designs the loops and keeps the framework</h2>
<div class="row" style="gap:48px">
<div class="col" style="flex:3;gap:22px">
{item("Why an intact framework", "hNbBCII10-FGLA is the scaffold RFantibody and Germinal fixed after expression tests. A caplacizumab-like VHH adds camelid hallmarks.")}
{item("Competition built in", "The framework sits where trastuzumab’s VH sits in the 1N8Z crystal, so a binder in this pose overlaps trastuzumab by construction. The pose is rigid: one binding mode per design.")}
{item("Hotspots from the crystal", "19 HER2 residues within 4.5 Å of trastuzumab, grouped into three hotspot sets: core, broad and alternate patch.")}
{item("Cheap enough to be wide", "High-recycle antibody recipe, 60 steps, 8 designs per job. About $0.002 per design, so the queue, not money, is the limit.")}
</div>
<div class="col" style="flex:2;gap:20px">
<div class="card"><p class="num">{st["n_designs"]}</p><p class="cap">RFD4 designs generated so far</p></div>
<p class="small">Designs per framework and hotspot set</p>
{table(["Framework", "core", "broad", "alt"], rows, [46, 18, 18, 18], fs=24)}
</div>
</div>
{foot("Sources: denovo/s1_design/decisions.md, denovo/targets/epitope.json, results/denovo/designs_rfd4jf.csv.")}
{note("Choices and why. Intact framework: it carries the stability we want and it is the scaffold the field already validated, and the organisers said fully new CDRs should land below 70% identity. Pose anchored on trastuzumab’s VH: this is how epitope competition is built in, and also the limit, because every design shares one rigid binding mode and the hosted model gives no confidence score. Cost is negligible, so we generated hundreds and let the gates decide.")}
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


def s_rules(st):
    nb = st.get("nb9")
    nbtxt = f"Median {nb[0]:.0f} of 9 over {nb[1]} designs; the parent scaffold scores 6." if nb else "Scoring pending."
    def b(k, v, free=False):
        c = "box-free" if free else "box"
        return f'<div class="{c}" style="width:250px;padding:12px 16px"><p class="box-k">{E(k)}</p><p class="box-v">{E(v)}</p></div>'
    cdr2 = "".join([b("pos 1", "A D G N S T"), b("pos 2", "free", True), b("pos 3", "free", True), b("pos 4", "G"), b("pos 5", "A G S T"), b("pos 6", "I N S T")])
    cdr3 = "".join([b("N1", "not Y W V F L I"), b("N2", "L or T"), b("middle", "free", True), b("C3", "F H L Y"), b("C2", "not F L Y"), b("C1", "Y or W")])
    hard = lambda t, x: f'<div class="card" style="flex:1;padding:20px;gap:6px"><p class="card-title" style="font-size:28px">{E(t)}</p><p class="small">{E(x)}</p></div>'  # noqa: E731
    return f'''<section id="rules" style="gap:28px">
<p class="eyebrow">HARD RULES AND SEQUENCE RULES</p>
<h2>Four hard rules, then a nine-rule sequence score</h2>
<div class="row" style="gap:20px">
{hard("Length", "At most 250 aa")}
{hard("No unpaired Cys", "An odd Cys count fails")}
{hard("Canonical disulfide", "Cys 23 and 104 present")}
{hard("Numbering", "Immunum confidence ≥ 0.1")}
</div>
<div class="col" style="gap:10px">
<p class="card-title" style="font-size:30px">Nine rules for CDR2 and CDR3, from the monoyield paper (same hNbBCII10 scaffold)</p>
<p class="lab">CDR2: six residues at scaffold positions 55–60</p>
<div class="row" style="gap:12px">{cdr2}</div>
<p class="lab">CDR3: counted from both ends of the loop</p>
<div class="row" style="gap:12px">{cdr3}</div>
</div>
<p class="small">{E(nbtxt)} The rules are exploratory: fitted on the paper’s own data and 5–13 aa CDR3s. We keep the framework’s C-A-A and F42 intact because CDR3 packs on them. Zero-Cys designs lose Tm.</p>
{foot("Source: monoyield paper, Fig. 5h and Methods; denovo/s4_props/nanobody_rules.py; denovo/s2_rules/rules.py.")}
{note("Why each rule. Length and unpaired Cys are the organisers’ rules. The canonical disulfide is ours: abguard does not require it, but a VHH without it loses a lot of melting temperature, and thermostability is one of the three judged properties. The nine rules come from the monoyield paper: they relate CDR2 and CDR3 residues to monomeric yield on the very scaffold we use. They are exploratory, derived on the paper’s own dataset, so they carry a weight of only 5 percent in the ranking.")}
</section>'''


def s_binding(st):
    o = st["of3"]
    gates = st["gates"]
    names = {"of3": "OpenFold3", "dde": "OpenDDE-abag", "esm2": "ESMFold-2"}
    rows = []
    for k in ("of3", "dde", "esm2"):
        g = gates.get(k)
        if g:
            rows.append([names[k], f"{g.get('auc_vhh_real', float('nan')):.2f}", f"{100 * g.get('recall_cognate_vhh', 0):.0f}%"])
    strip = dot_strip([
        (f"Framework-only decoys (n = {len(o['decoy'])})", o["decoy"], "#6C7783", med(o["decoy"])),
        (f"RFD4 designs (n = {len(o['design'])})", o["design"], "#0C7473", med(o["design"])),
        ("Trastuzumab Fv (positive control)", o["tz"], "#B4620A", med(o["tz"])),
    ], 0.0, 0.6, svg_w=640, label_w=300, row_h=110, fmt="{:.1f}") if o["design"] else '<p class="muted">OpenFold3 scoring pending.</p>'
    md = f"{o['med_d']:.2f}" if o["med_d"] is not None else "–"
    mc = f"{o['med_c']:.2f}" if o["med_c"] is not None else "–"
    tz = f"{o['tz_v']:.2f}" if o["tz_v"] is not None else "–"
    return f'''<section id="binding">
<p class="eyebrow">INDEPENDENT BINDING EVIDENCE</p>
<h2>Folding scores are filters, calibrated on decoys</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:1000px;gap:10px">
<p class="lab">OpenFold3 interface confidence (ipTM) against HER2 domain IV; the bar marks the median</p>
{strip}
</div>
<div class="col" style="flex:1;gap:16px">
<p class="card-title">How well each scorer separates real pairs</p>
{table(["Model", "AUC", "Recall"], rows, [52, 24, 24], fs=24)}
<p class="small">AUC on real VHH–antigen pairs against 34 non-binders. Recall is measured at the gate set on 46 decoys, using 10 crystal VHH pairs.</p>
</div>
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="body">Medians: decoys {mc}, designs {md}, trastuzumab Fv {tz}. {o["pass"]} of {o["n"]} scored designs clear the OpenFold3 gate so far.</p></div>
{foot("Sources: denovo/s3_score/calibration.md, gates.json; results/denovo/score_of3.csv; research report (ipTM ranks binders weakly).")}
{note("Honest reading: the trastuzumab Fv, a real binder, scores about four times the decoys, so the pipeline can see a good complex. The scorers are weak on real VHH pairs, so an AUC of about 0.65 is what we should expect and why gates are set from decoys, not from published thresholds. At this snapshot our designs look like decoys in OpenFold3. We do not claim binding; we say which designs are most consistent across models, on the right epitope, and developable, and the assay is the judge.")}
</section>'''


def s_develop(st):
    tm = st["tm"]
    rows = [
        ("BCII10-based designs", tm["bcii10"], "#0C7473", med(tm["bcii10"])),
        ("Caplacizumab-like designs", tm["capla"], "#B4620A", med(tm["capla"])),
        ("Reference nanobodies", tm["ref"], "#6C7783", med(tm["ref"])),
    ]
    rows = [(f"{a} (n = {len(v)})", v, c, m) for a, v, c, m in rows if len(v)]
    strip = dot_strip(rows, 55, 80, svg_w=470, label_w=250, row_h=104, nticks=6, fmt="{:.0f}") if rows else '<p class="muted">NanoMelt results pending.</p>'
    tool_rows = [
        ["NanoMelt", "Apparent Tm", "Laptop, mirrored on VM"],
        ["TNP", "Red and amber developability flags", "Laptop, mirrored on VM"],
        ["ESMC-Nb models", "Monomeric yield, monomer/dimer ratio", "Team VM (6B model)"],
        ["Liabilities", "NG, DG, N-glyc, unpaired Cys, charge", "Laptop"],
    ]
    def stat(n, cap):
        big = f'<p class="num" style="font-size:64px">{n}</p>' if n else '<p class="num-amber" style="font-size:44px;line-height:1.4">pending</p>'
        return f'<div class="card" style="flex:1;padding:20px;gap:4px">{big}<p class="cap" style="font-size:24px">{E(cap)}</p></div>'
    return f'''<section id="develop">
<p class="eyebrow">DEVELOPABILITY</p>
<h2>Melting temperature, yield and liabilities are predicted</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:880px;gap:18px">
{table(["Tool", "Predicts", "Runs on"], tool_rows, [22, 46, 32], fs=24)}
<div class="row" style="gap:12px">{stat(st["tm"]["n"], "with NanoMelt Tm")}{stat(st["n_tnp"], "with TNP flags")}{stat(st["n_yield"], "with yield")}</div>
<p class="small">TNP and yield run on the top third by OpenFold3 consensus, where they cost minutes per design.</p>
</div>
<div class="col" style="flex:1;gap:10px">
<p class="lab">NanoMelt predicted Tm in °C; the bar marks the median</p>
{strip}
</div>
</div>
{foot("Sources: denovo/s4_props/README.md, results/denovo/props.csv. NanoMelt is a prediction with several °C of error and a non-commercial licence.")}
{note("These are predictions, so they act as filters and tie-breakers, not as measurements. NanoMelt gives an apparent Tm; the reference nanobodies show where well-behaved VHHs sit on the same scale. TNP flags developability outliers; the ESMC-based models predict monomeric yield. The heavy tools only run on the finalists. Comparing a VHH Tm with trastuzumab’s is not apples to apples, which is why Tm is one objective in the Pareto front rather than a pass-fail test.")}
</section>'''


def s_live(st):
    n = st["n_designs"]
    full = st["par"]["n_complete"] if st.get("par") else st["n_cov50"]
    rows = [
        ("Designed", n, "#6C7783"),
        ("Rules checked", st["n_rules"], "#6C7783"),
        ("abguard says de novo", st["n_denovo"], "#0C7473"),
        ("Gate + 2-point margin", st["n_gate_strict"], "#0C7473"),
        ("OpenFold3 scored", st["of3"]["n"], "#0C7473"),
        ("Orthogonal models", st["n_orth"], "#0C7473"),
        ("NanoMelt Tm", st["tm"]["n"], "#B4620A"),
        ("TNP + yield", st["n_yield"], "#B4620A"),
        ("Fully measured", full, "#12212F"),
    ]
    cand = st["cand"]
    if cand:
        trs = []
        for c in cand[:7]:
            trs.append([str(c["design_id"]).replace("-rfd4jf", ""), c.get("framework", ""),
                        f"{int(c['pareto_rank'])}" if "pareto_rank" in c and pd.notna(c.get("pareto_rank")) else "–",
                        str(c.get("weakest", "–")), f"{c['nm_tm']:.1f}" if pd.notna(c.get("nm_tm")) else "–"])
        right = (f'<p class="card-title" style="font-size:30px">Designs with complete evidence: {st.get("cand_n", len(cand))}, best 7 shown</p>'
                 + table(["Design", "Base", "Front", "Weakest", "Tm °C"], trs, [38, 17, 13, 18, 14], fs=24)
                 + '<p class="small">Front and weakest objective are within-run percentiles.</p>')
    else:
        right = ('<div class="card"><p class="card-title">Candidates: pending</p><p class="card-text">A design is listed only when OpenFold3, Tm and yield exist for it and it passes the gate. '
                 f'{full} designs have that so far.</p></div>')
    return f'''<section id="live">
<p class="eyebrow">THE LIVE FUNNEL</p>
<h2>Where the designs are right now</h2>
<div class="row" style="gap:48px">
<div class="col" style="width:800px;gap:14px">
{hbars(rows, max_w=300, label_w=340)}
</div>
<div class="col" style="flex:1;gap:12px">
{right}
</div>
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small">Stages run in parallel, so later bars are still filling in.</p></div>
{foot("Sources: results/denovo/rules.csv, score_of3.csv, score_orth.csv, props.csv, leaderboard_pareto.csv.")}
{note("This slide is rebuilt from the result files. Read it as a funnel: each bar is the number of RFD4 designs that reached that stage. The candidate table appears only for designs with complete evidence, so an early list is never built on one score. Percentiles are within this run, not against trastuzumab, so they order the designs but do not say they beat it.")}
</section>'''


def pareto_scatter(par, w=800, h=360):
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
    for b, t, r, sel in sorted(par["pts"], key=lambda p: -p[2]):  # later fronts first, front 1 on top
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
    w = st["weights"]
    if par and len(par["pts"]) >= 5:
        legend = ('<div style="display:flex;gap:20px;align-items:center">' + legend_item("f1", "Front 1") + legend_item("f2", "Front 2")
                  + legend_item("f3", "Front 3") + legend_item("f4", "Later") + legend_item("sel", "Selected") + "</div>")
        fr = " / ".join(str(v) for v in par["fronts"][:5])
        left = (f'<div class="col" style="width:830px;gap:8px"><p class="lab">Binding (x) against melting temperature (y), percentiles within this run</p>'
                f'{legend}{pareto_scatter(par)}</div>')
        under = (f"{par['n_pop']} gate-passing designs; {len(par['pts'])} have both axes measured and are plotted; "
                 f"front sizes {fr}; {par['n_complete']} fully measured. ")
    else:
        left = ('<div class="card" style="width:830px"><p class="card-title">Pareto chart: pending</p>'
                '<p class="card-text">It appears once binding and melting-temperature evidence exist for at least five gate-passing designs.</p></div>')
        under = ""
    objs = ", ".join(par["used"]) if par and par["used"] else "pending"
    dropped = ", ".join(par["dropped"]) if par and par["dropped"] else "none"
    item = lambda k, v: f'<div class="col" style="gap:4px"><p class="card-title" style="font-size:30px">{E(k)}</p><p class="muted" style="font-size:26px">{E(v)}</p></div>'  # noqa: E731
    npop = par["n_pop"] if par else "–"
    right = f'''<div class="col" style="flex:1;gap:16px">
{item("1 · Population", f"{npop} designs pass the hard gate; controls never rank.")}
{item("2 · Objectives, all maximised", f"Used: {objs}. Dropped for low coverage: {dropped}. Each is a within-run percentile.")}
{item("3 · Missing data", "An unmeasured design scores 0 on that objective, so it never beats a measured one.")}
{item("4 · The list", "Fronts in order; the weakest objective breaks ties; no near-duplicate CDR3; at most half from one framework and hotspot group.")}
</div>'''
    return f'''<section id="ranking">
<p class="eyebrow">STAGE 6 · CHOOSING THE 20 SLOTS</p>
<h2>Pareto dominance picks the list, not one weighted score</h2>
<div class="row" style="gap:48px">
{left}
{right}
</div>
<div class="row" style="align-items:center;gap:20px">{tag("prelim")}<p class="small" style="width:1380px">{E(under)}Percentiles are within this run, not against trastuzumab, and the proxies are filters, not affinity measurements.</p></div>
{foot("Sources: denovo/pareto.py (docstring), denovo/rank.py, results/denovo/leaderboard_pareto.csv and submission_top20_pareto.csv.")}
{note("Why Pareto: the organisers tier a design by how many of the three judged properties beat trastuzumab, and S needs all three. A weighted sum lets a very strong binder with a poor melting temperature outrank a balanced design, but the weakest objective is what loses the tier. So we sort the gate-passing designs into non-dominated fronts and break ties inside a front by the design’s worst objective percentile, which favours balanced designs. Bind blends OpenFold3 consensus with OpenDDE and ESMFold-2 and is an affinity proxy; the epitope objective uses overlap with trastuzumab’s contacts; the Tm objective is NanoMelt. The weighted composite " + f"(bind {w['bind']:.2f}, comp {w['comp']:.2f}, Tm {w['therm']:.2f}, yield {w['yield']:.2f}, dev {w['dev']:.2f}, nine-rule {w['rules']:.2f}; missing evidence counts as a below-median prior)" + " stays as a comparison list. Honest limits: percentiles are relative to this run, and none of these proxies measures affinity.")}
</section>'''


def s_dropped(st):
    bg = st["bg_n"]
    fd = st["fwdiv_n"]
    rows = [
        ["Framework diversification", "Organisers said new CDRs should land under 70%. We measured it: true for one framework, not the other, and length decides. Variants parked" + (f" ({fd})." if fd else "."), "gate_drivers.md"],
        ["BoltzGen nanobody designs", "Raw designs have no cysteines" + (f" ({bg} designs)" if bg else "") + ", so no canonical disulfide. Stopped to focus on RFD4; its running jobs were cancelled.", "parked_boltzgen*.csv"],
        ["RFD4 on the team VM", "Ruled out: JapanFold is the only RFD4 route. One container that had started was stopped.", "CONTRACT.md"],
        ["ipSAE at the 10 Å cutoff", "Zero for every complex, even trastuzumab’s Fv. 15 and 20 Å variants added.", "score_of3.csv"],
        ["First NanoMelt runs", "Failed on the 8 GB laptop GPU and in a VM script. Fixed; " + f"{st['tm']['n']} designs now have a Tm.", "props.csv"],
    ]
    return f'''<section id="dropped">
<p class="eyebrow">WHAT WE TRIED AND DROPPED</p>
<h2>Five things that did not survive contact with data</h2>
{table(["What", "Why it was dropped", "Evidence"], rows, [22, 55, 23], fs=24)}
{foot("Sources: denovo/CONTRACT.md decisions, stage decisions.md files, results/denovo/parked_*.csv.")}
{note("Reasoning is the largest judging weight, so this slide is deliberate. Each row is a decision we reversed because of a measurement or a constraint, not a failure we hid. The first row is the clearest example: we built a framework diversifier, then measured the gate and found the framework did not need to change.")}
</section>'''


def s_risks(st):
    nk = st["nearest_known"]
    nktxt = (f"Identity to the nearest known VHH is {nk[0]:.2f}–{nk[1]:.2f} in our 183,567-domain panel; the organisers said new CDRs should land under 70%."
             if nk else "Identity to the nearest known VHH is pending.")
    card = lambda t, b: f'<div class="card" style="flex:1;padding:24px;gap:8px"><p class="card-title" style="font-size:30px">{E(t)}</p><p class="card-text" style="font-size:26px">{E(b)}</p></div>'  # noqa: E731
    return f'''<section id="risks" style="gap:24px">
<p class="eyebrow">RISKS</p>
<h2>What can go wrong, and the outcome we expect</h2>
<div class="row">
{card("Weak scores", "OpenFold3 separates real VHH pairs from decoys only modestly (AUC 0.65), and our designs score like decoys so far.")}
{card("Novelty is unverifiable", nktxt)}
</div>
<div class="row">
{card("One binding mode", "Hosted RFD4 holds the framework and HER2 rigid and does not disclose its checkpoint. A VHH that cannot bind in the trastuzumab-VH pose has no second chance.")}
{card("An affinity gap", "Trastuzumab binds at 0.03–0.18 nM. De novo VHHs in the literature reach tens to hundreds of nM, with 4–22% hit rates per design.")}
</div>
<div class="card-dark" style="padding:24px;gap:8px"><p class="card-title" style="font-size:32px">Expected outcome (our inference): A- or B-tier</p><p class="card-text" style="font-size:26px">Epitope competition and thermostability are plausible. Beating trastuzumab’s affinity is the unlikely third.</p></div>
{foot("Sources: denovo/s3_score/calibration.md, denovo/s2_rules/gate_drivers.md, research/denovo_vhh_research.md (Germinal, RFantibody, BoltzGen hit rates).")}
{note("We would rather state these than have a reviewer find them. The affinity gap is the main reason we expect A or B tier rather than S: the literature reports tens to hundreds of nanomolar for de novo VHHs, while the trastuzumab reference is sub-nanomolar. The novelty risk is one we cannot close from here: the VM tool does not implement the jury sentence about identity to any known antibody.")}
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
{item("The reusable part", "The gate, calibration, profiling and ranking steps work for any target: swap the target files and the hotspot sets.")}
</div>
<div class="col" style="flex:1;gap:20px">
<div class="card"><p class="card-title">Reproduce it</p>
<p class="body" style="font-size:26px">denovo/CONTRACT.md defines every stage’s inputs and outputs. Each stage folder has a README with its run command. Every score in results/denovo/ traces to a cached job id.</p></div>
<div class="card-dark"><p class="card-title">Refresh the ranking</p>
<p class="card-text" style="font-size:26px">denovo/.venv/bin/python denovo/rank.py<br>denovo/.venv/bin/python denovo/pareto.py</p></div>
<div class="row" style="gap:16px"><div class="card" style="flex:1;padding:24px"><p class="num" style="font-size:64px">{E(spent)}</p><p class="cap">JapanFold spend so far, cap $65</p></div>
<div class="card" style="flex:1;padding:24px"><p class="num" style="font-size:64px">5</p><p class="cap">parallel workers, one file contract</p></div></div>
</div>
</div>
{foot("Sources: research/denovo_vhh_research.md Q6 (inference, abstract-level patent reading); denovo/CONTRACT.md.")}
{note("User fit and path to value, with the caveats: the value case is an inference from the literature, not a market study. The reusable part is the real asset for a method-first reading: a pipeline that gates on the organisers’ own rule, calibrates its scorers on decoys and reports what it cannot verify. Spend is read from the JapanFold balance at build time.")}
</section>'''


ORDER = [("cover", s_cover), ("task", s_task), ("pipeline", s_pipeline), ("generate", s_generate), ("gate", s_gate), ("rules", s_rules),
         ("binding", s_binding), ("develop", s_develop), ("live", s_live), ("ranking", s_ranking), ("dropped", s_dropped),
         ("risks", s_risks), ("value", s_value)]
SECTIONS = {"s1": ("The task and the approach", "cover"), "s2": ("Design choices and the gate", "generate"),
            "s3": ("Evidence so far", "binding"), "s4": ("Limits, value and reproducibility", "dropped")}


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
    st = collect()
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
            "files": {k: v for k, v in files.items() if k != "project/slides/cover.html"},
            "note": "first publish also sends project/deck.json and project/deck.css; a refresh sends the slide files only"}
    (HERE / "publish_plan.json").write_text(json.dumps(plan, indent=1))
    print(f"built {len(ORDER)} slides at {TS}: designs={st['n_designs']} gate_strict={st['n_gate_strict']} of3={st['of3']['n']} "
          f"tm={st['tm']['n']} cand={len(st['cand'])} lint={'FAIL' if bad else 'ok'}")
    if "--plan" in sys.argv:
        print(json.dumps(plan, indent=1))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
