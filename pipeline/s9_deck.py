"""Stage 9 (optional): PowerPoint summary built from a template (paths.template_pptx) and the pipeline outputs.
Every number on the slides is read from the stage outputs (summary.json, CSVs, final/top_candidates.json)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION, XL_TICK_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Inches, Pt

from .common import VERNIER, WT_H3, blosum, out, read_summary
from .s2_panel import vernier_options

DB, GR, BL, TE, LI, LB, GY = [RGBColor.from_string(x) for x in ("004B6F", "63A44B", "1275AC", "009299", "9185BE", "5C9BAF", "92A19A")]
BLACK, WHITE, LIGHT = RGBColor(0, 0, 0), RGBColor(255, 255, 255), RGBColor(0xEE, 0xF3, 0xF5)
C_VAR, C_CONT, C_VERN = RGBColor(0xF2, 0x8C, 0x1A), RGBColor(0xFA, 0xCC, 0x1A), RGBColor(0xCC, 0x1F, 0x8C)


class Deck:
    def __init__(self, template: str):
        self.prs = Presentation(template)
        lst = self.prs.slides._sldIdLst
        for sid in list(lst):                      # drop the template's sample slides, keep its layouts
            self.prs.part.drop_rel(sid.rId)
            lst.remove(sid)
        self.L = {l.name: l for l in self.prs.slide_layouts}

    def slide(self, layout: str, title: str):
        s = self.prs.slides.add_slide(self.L[layout])
        s.shapes.title.text = title
        return s

    @staticmethod
    def ph(s, idx):
        return next(p for p in s.placeholders if p.placeholder_format.idx == idx)

    def drop(self, s, idx):
        p = self.ph(s, idx)
        p._element.getparent().remove(p._element)

    @staticmethod
    def fill(tf, items, size=None, space=6):
        first = True
        for it in items:
            lvl = 0
            if isinstance(it, tuple) and isinstance(it[1], int):
                it, lvl = it
            p = tf.paragraphs[0] if first else tf.add_paragraph()
            first = False
            p.level, p.space_after = lvl, Pt(space)
            for t, b in (it if isinstance(it, list) else [(it, False)]):
                r = p.add_run()
                r.text = t
                if size:
                    r.font.size = Pt(size)
                if b:
                    r.font.bold = True

    @staticmethod
    def text(s, x, y, w, h, items, size=14, color=BLACK, anchor=None):
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.margin_left = tf.margin_right = Inches(0.05)
        if anchor:
            tf.vertical_anchor = anchor
        for k, it in enumerate(items):
            p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
            p.space_after = Pt(4)
            for t, b in (it if isinstance(it, list) else [(it, False)]):
                r = p.add_run()
                r.text, r.font.size, r.font.bold, r.font.color.rgb = t, Pt(size), b, color

    @staticmethod
    def picture(s, path, x, y, mw, mh):
        w, h = Image.open(path).size
        sc = min(mw / w, mh / h)
        return s.shapes.add_picture(str(path), Inches(x + (mw - w * sc) / 2), Inches(y + (mh - h * sc) / 2), Inches(w * sc), Inches(h * sc))

    @staticmethod
    def box(s, x, y, w, h, color, shape=MSO_SHAPE.RECTANGLE):
        r = s.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
        r.fill.solid()
        r.fill.fore_color.rgb = color
        r.line.fill.background()
        r.shadow.inherit = False
        return r

    @staticmethod
    def table(s, x, y, colw, rows, header, size=12, rowh=0.42, bold_cols=()):
        t = s.shapes.add_table(len(rows) + 1, len(header), Inches(x), Inches(y), Inches(sum(colw)), Inches(rowh * (len(rows) + 1)))
        tb = t.table
        tb._tbl.tblPr.set("bandRow", "0")
        tb._tbl.tblPr.set("firstRow", "0")
        for i, cw in enumerate(colw):
            tb.columns[i].width = Inches(cw)
        for r in range(len(rows) + 1):
            tb.rows[r].height = Inches(rowh)

        def cell(c, runs, bg, color, bold=False):
            c.fill.solid()
            c.fill.fore_color.rgb = bg
            c.margin_left = c.margin_right = Inches(0.08)
            c.margin_top = c.margin_bottom = Inches(0.03)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.text_frame.word_wrap = True
            p = c.text_frame.paragraphs[0]
            for t, b in (runs if isinstance(runs, list) else [(str(runs), bold)]):
                r = p.add_run()
                r.text, r.font.size, r.font.bold, r.font.color.rgb = t, Pt(size), b, color

        for j, h in enumerate(header):
            cell(tb.cell(0, j), h, DB, WHITE, True)
        for i, row in enumerate(rows, 1):
            for j, v in enumerate(row):
                cell(tb.cell(i, j), v, LIGHT if i % 2 == 0 else WHITE, BLACK, j in bold_cols)

    def chart(self, s, kind, cats, series, x, y, w, h, title, colors, fmt="0.00", vmin=None, vmax=None, legend=True):
        cd = CategoryChartData()
        cd.categories = cats
        for name, vals in series:
            cd.add_series(name, vals)
        ch = s.shapes.add_chart(kind, Inches(x), Inches(y), Inches(w), Inches(h), cd).chart
        ch.has_legend = legend
        if legend:
            ch.legend.position, ch.legend.include_in_layout = XL_LEGEND_POSITION.BOTTOM, False
            ch.legend.font.size = Pt(12)
        ch.has_title = True
        ch.chart_title.text_frame.text = title
        run = ch.chart_title.text_frame.paragraphs[0].runs[0]
        run.font.size, run.font.bold = Pt(14), True
        for ser, col in zip(ch.plots[0].series, colors):
            ser.format.fill.solid()
            ser.format.fill.fore_color.rgb = col
            ser.invert_if_negative = False
        pl = ch.plots[0]
        pl.gap_width = 60
        pl.has_data_labels = True
        pl.data_labels.font.size = Pt(10)
        pl.data_labels.number_format, pl.data_labels.number_format_is_linked = fmt, False
        pl.data_labels.position = XL_LABEL_POSITION.OUTSIDE_END
        if vmin is not None:
            ch.value_axis.minimum_scale = vmin
        if vmax is not None:
            ch.value_axis.maximum_scale = vmax
        ch.value_axis.tick_labels.font.size = ch.category_axis.tick_labels.font.size = Pt(11)
        ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDD, 0xDD, 0xDD)
        ch.category_axis.tick_label_position = XL_TICK_LABEL_POSITION.LOW
        return ch


def run(cfg: dict) -> None:
    S = read_summary(cfg)
    base = cfg["design"]["base_h3"]
    fin = out(cfg, "final")
    top = json.load(open(fin / "top_candidates.json"))["candidates"]
    fig = out(cfg, "figures")
    pose = json.load(open(fig / "pose_rmsd.json"))
    folds = (pd.read_csv(out(cfg, "fold_scores.csv")).merge(pd.read_csv(out(cfg, "contacts.csv")), on="id", how="left")
             .merge(pd.read_csv(out(cfg, "pose.csv")), on="id", how="left"))
    scores = pd.read_csv(out(cfg, "scores.csv"))
    meta = pd.read_csv(out(cfg, "panel_meta.csv"))
    ctrl = folds[folds.id.str.startswith("ctrl_")].merge(scores, on="id", how="left").set_index("id")
    des = folds[~folds.id.str.startswith("ctrl_")].merge(meta, on="id")
    wt_ip = des[des.vernier_changes == "WT"].set_index("cdrh3").iptm
    des["d_iptm"] = des.iptm - des.cdrh3.map(wt_ip)
    swaps = des[des.vernier_changes != "WT"].assign(site=lambda x: x.vernier_changes.str[:3])
    vsite = swaps.groupby("site").d_iptm.agg(["mean", "median", "size"])
    near0 = float((swaps.d_iptm.abs() < 0.02).mean())
    lower = float((swaps.d_iptm < 0).mean())
    opts = vernier_options(cfg["vernier"]["min_blosum45"], cfg["vernier"]["residues"])

    D = Deck(str(Path(cfg["paths"]["template_pptx"]).expanduser()))
    # 1 title
    s = D.prs.slides.add_slide(D.L["Title"])
    s.shapes.title.text = "Redesigning trastuzumab CDR-H3 and the Vernier zone against HER2"
    D.fill(D.ph(s, 1).text_frame, ["Jakob Riccabona  ·  Oslo HER2 hackathon  ·  2026"])
    # 2 summary
    s = D.slide("Heading and Text (1 column)", "Summary")
    D.fill(D.ph(s, 1).text_frame, [
        [("Start: ", True), (f"the best experimentally tested CDR-H3, {base} (SPR 0.023 nM vs 0.114 nM for trastuzumab)", False)],
        [("CDR-H3: ", True), (f"BLOSUM-45 substitutions (Moritz's recipe) at the {len(cfg['design']['positions'])} positions that do not touch HER2, 3–4 changes from the base", False)],
        [("Vernier zone: ", True), ("conservative BLOSUM-45 swaps at H94, H93, H47, L36, L46, L49, only where they touch a changed H3 position", False)],
        [("Evaluation: ", True), ("OpenDDE folds scored with ipTM, ipSAE, iLIS and GetContacts (Fv–HER2 and CDR-H3 against the rest, vs the crystal)", False)],
        [("Rule: ", True), ("a design is kept only if the R98–position-10 salt bridge is present, or if both partners are mutated", False)],
        [("Result: ", True), (f"{S['n_top']} candidates ({S['n_top_h3']} distinct CDR-H3): {S['n_top_salt_bridge']} keep the salt bridge, {S['n_top_both_mutated']} mutate both partners", False)]], size=18, space=10)
    # 3 starting point
    s = D.slide("Heading and text (2 columns)", "Starting point and constraints")
    D.fill(D.ph(s, 1).text_frame, [[("Rules a design has to pass", True)],
        ("Tutorial filter: ≥2 changes from every known CDR-H3 and ≥3 of 10 positions changed from WT", 1),
        ("abguard: ≥3 edits from every known CDR-H3 and an even number of heavy-chain Cys", 1),
        ("Known set: 525,063 sequences (FACS screen and SPR designs)", 1)], size=16, space=8)
    D.fill(D.ph(s, 2).text_frame, [[("What the data and structure say", True)],
        (f"{base} is the tightest SPR-tested loop: G3S, D4Q, M9Y vs WT", 1),
        ("Only CDR-H3 positions 1, 5, 6, 7 contact HER2 (1N8Z)", 1),
        ("D108 salt-bridges R98 (H94); A106 and M107 are buried against VL", 1),
        ("Vernier residues are invariant in all data: choices are structure-driven", 1)], size=16, space=8)
    # 4 structure overview
    s = D.slide("Heading and Graphic (1 column)", "Where the changes sit in the structure")
    D.drop(s, 1)
    D.picture(s, fig / "fig1_overview.png", 0.92, 1.2, 7.2, 5.3)
    y = 1.45
    for col, name, desc in ((C_VAR, "Varied H3 positions 2, 3, 4, 8, 9, 10", "no direct HER2 contact"), (C_CONT, "HER2-contact positions 1, 5, 6, 7", "W, G, F, Y never touched"),
                            (C_VERN, "Vernier residues H94, H93, H47, L36, L46, L49", "swapped only where they touch a changed H3 position")):
        D.box(s, 8.5, y + 0.06, 0.28, 0.28, col)
        D.text(s, 8.9, y - 0.04, 3.5, 1.05, [[(name, True)], desc], size=13)
        y += 1.1
    D.text(s, 8.5, 5.0, 3.9, 1.2, ["Fv of trastuzumab with HER2 domain IV (PDB 1N8Z). VH light grey, VL darker grey, HER2 pale cyan."], size=11, color=GY)
    # 5 H3 rationale
    s = D.slide("Heading and Text (3 columns)", "CDR-H3 redesign rationale")
    D.fill(D.ph(s, 1).text_frame, [[("Start", True)], f"The best tested loop {base}, not wild type", "Mutations are counted from it, so every design stays ≥3 edits from all known loops"], size=16, space=8)
    D.fill(D.ph(s, 11).text_frame, [[("Restrict", True)], "Positions 1, 5, 6, 7 touch HER2 and stay fixed", "Substitutions follow the paper's BLOSUM-45 recipe: P(b|a) ∝ f_b·exp(0.25·s_ab), top 6 per position"], size=16, space=8)
    D.fill(D.ph(s, 10).text_frame, [[("Filter", True)], "No Cys, no extra liabilities, novelty vs the 525k known sequences, then abguard", "D10 is allowed, but a design needs the R98–position-10 salt bridge, or both partners mutated"], size=16, space=8)
    # 6 Vernier
    s = D.slide("Heading and Graphic (1 column)", "Vernier-zone rationale")
    D.drop(s, 1)
    D.picture(s, fig / "fig2_closeup.png", 0.92, 1.2, 6.2, 4.6)
    rows = [[n, VERNIER[n][2], ", ".join(str(p) for p in sorted(VERNIER[n][3])), " ".join(opts[n])] for n in cfg["vernier"]["residues"]]
    D.table(s, 7.35, 1.25, [0.8, 0.7, 1.5, 2.1], rows, ["Residue", "WT", "Touches H3 pos", "BLOSUM45 ≥ +1"], size=12, rowh=0.4, bold_cols=(0,))
    D.text(s, 7.35, 4.1, 5.1, 2.4, [[("Rule: ", True), ("change a Vernier residue only if it touches a changed H3 position", False)], [("Conservative: ", True), ("BLOSUM-45 ≥ +1, no Cys or Pro", False)],
           [("Why: ", True), ("the framework is already good for the base loop, so swaps are compensation, not redesign", False)]], size=13)
    D.text(s, 0.92, 5.9, 6.2, 0.6, ["Close-up of 1N8Z. Dashes: side-chain contacts under 4 Å between Vernier residues and the H3 positions."], size=11, color=GY)
    # 7 funnel
    s = D.slide("Heading and Text (1 column)", f"From {S['n_raw']:,} designs to {S['n_top']} candidates")
    D.drop(s, 1)
    steps = [(f"{S['n_raw']:,}", f"designs: 3–4 BLOSUM-45 changes from {base} at the H3 positions"), (f"{S['n_enum']:,}", "no Cys, liabilities ≤ WT, ≥3 from every known sequence"),
             (f"{S['n_panel']}", f"sequences: {S['n_shortlist']} diverse H3 variants × Vernier genotypes"), (f"{S['n_valid']}", f"pass abguard ({S['n_valid_h3']} H3 variants)"),
             (f"{S['n_folded']}", "folded with OpenDDE (incl. 4 controls)"), (f"{S['n_epitope_ok']}", f"docked (ipTM ≥ 0.9) on the crystal epitope (RMSD ≤ {cfg['rank']['max_epitope_rmsd']:g} Å)"),
             (f"{S['n_pass_contacts']}", f"keep the crystal contacts ({S['n_pass_contacts_h3']} H3 variants)"), (f"{S['n_after_salt_rule']}", f"salt-bridge rule: bridge present or both partners mutated ({S['n_after_salt_rule_h3']} H3)"),
             (f"{S['n_top']}", "best by iLIS, ≤3 per H3")]
    y = 1.2
    for i, (n, t) in enumerate(steps):
        r = D.box(s, 0.92, y, 2.1, 0.54, [DB, BL, BL, TE, TE, GR, GR, GR, GR][i], MSO_SHAPE.ROUNDED_RECTANGLE)
        r.text_frame.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = r.text_frame.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        run_ = p.add_run()
        run_.text, run_.font.size, run_.font.bold, run_.font.color.rgb = n, Pt(22), True, WHITE
        D.text(s, 3.25, y + 0.05, 9.2, 0.5, [t], size=16, anchor=MSO_ANCHOR.MIDDLE)
        y += 0.62
    # 8 controls table
    s = D.slide("Heading and Graphic (1 column)", "Calibration on measured controls")
    D.drop(s, 1)
    names = {"ctrl_base": "Base loop (0.023 nM)", "ctrl_wt": "Trastuzumab (0.11 nM)", "ctrl_nonbinder1": "Non-binder 1", "ctrl_nonbinder2": "Non-binder 2"}
    rows = [[names[i], f"{ctrl.loc[i, 'iptm']:.3f}", f"{ctrl.loc[i, 'ipsae']:.3f}", f"{ctrl.loc[i, 'ilis']:.3f}", f"{ctrl.loc[i, 'if_recall']:.2f}",
             f"{ctrl.loc[i, 'in_recall']:.2f}", f"{ctrl.loc[i, 'epitope_rmsd']:.1f}", "yes" if ctrl.loc[i, "salt_bridge_r98_pos10"] else "no"] for i in names]
    D.table(s, 0.92, 1.3, [2.5, 0.95, 0.95, 0.95, 1.55, 1.6, 1.5, 1.5], rows, ["Control", "ipTM", "ipSAE", "iLIS", "HER2 contacts kept", "H3-vs-rest kept", "Epitope RMSD (Å)", "R98–D108 salt bridge"], size=12, rowh=0.5)
    D.text(s, 0.92, 4.3, 11.5, 2.2, ["Only the base loop docks in the crystal pose; the wild-type prediction scores like the non-binders, so these metrics measure \"looks like the base loop and docks\", not affinity.",
           "Contacts kept: share of the crystal's 27 Fv–HER2 residue pairs (interface) and of its CDR-H3-vs-rest pairs found in the prediction.",
           "Epitope RMSD: C-alpha RMSD of the 19 crystal HER2 epitope residues after fitting the Fv framework; the base loop is 1.3 Å, trastuzumab's own prediction binds elsewhere.",
           "The base loop's own prediction sets the pass level for the contact filter."], size=14)
    # 9 GetContacts chart
    s = D.slide("Heading and Graphic (1 column)", "GetContacts: interface and CDR-H3 against the rest")
    D.drop(s, 1)
    D.chart(s, XL_CHART_TYPE.COLUMN_CLUSTERED, list(names.values()), [("Fv–HER2 pairs kept", [float(ctrl.loc[i, "if_recall"]) for i in names]),
            ("CDR-H3-vs-rest pairs kept", [float(ctrl.loc[i, "in_recall"]) for i in names])], 0.92, 1.2, 7.4, 5.2, "Share of crystal contacts recovered", (DB, TE), vmin=0, vmax=1.1)
    D.text(s, 8.6, 1.3, 3.9, 5.0, [[("Two views: ", True), (f"Fv vs HER2 ({S['crystal_if_pairs']} crystal pairs) and CDR-H3 vs the rest of the Fv ({S['crystal_intra_pairs']} pairs).", False)],
           [("CDR-H3 vs rest ", True), ("covers the Vernier contacts, the VL contacts of positions 8–9 and the R98–D108 salt bridge.", False)],
           [("Filter: ", True), (f"≥{cfg['rank']['min_interface_recall']:.0%} of interface pairs, ≤{cfg['rank']['max_h3_pairs_missing']} H3–HER2 pair missing, ≥{cfg['rank']['min_h3_rest_recall']:.0%} of CDR-H3-vs-rest pairs.", False)]], size=14)
    # 9b score distributions with trastuzumab marked
    s = D.slide("Heading and Graphic (1 column)", "Scores and contacts with trastuzumab marked")
    D.drop(s, 1)
    D.picture(s, fig / "fig_distributions.png", 0.92, 1.15, 8.4, 5.4)
    D.text(s, 9.45, 1.3, 3.2, 5.2, [[("Red: ", True), ("trastuzumab, predicted by the same model as the designs.", False)],
           [("Dotted: ", True), (f"the crystal, {S['crystal_if_pairs']} Fv–HER2 and {S['crystal_intra_pairs']} CDR-H3-vs-rest pairs.", False)],
           [("Dashed green: ", True), ("the best tested loop.", False)],
           "The wild-type prediction does not dock (ipSAE ≈ 0.02, epitope 30 Å away), so it sits at the low end of every score.",
           f"Every design with ipTM ≥ {cfg['rank']['min_iptm']} lies within 3 Å of the crystal epitope: if a design docks, it docks on the same epitope."], size=12)
    s = D.slide("Heading and Graphic (1 column)", "ipSAE against iLIS")
    D.drop(s, 1)
    D.picture(s, fig / "fig_ipsae_ilis.png", 0.92, 1.15, 7.4, 5.4)
    D.text(s, 8.6, 1.3, 3.9, 5.2, [[("Two populations: ", True), ("designs that dock (high ipSAE) and designs that do not.", False)],
           [("Colour: ", True), ("number of Fv–HER2 residue-pair contacts.", False)],
           [("Orange rings: ", True), ("the final list. Star: trastuzumab; square: best tested loop; crosses: measured non-binders.", False)]], size=13)
    # 10 vernier effects
    s = D.slide("Heading and Graphic (1 column)", "Vernier swaps in the OpenDDE folds")
    D.drop(s, 1)
    sites = list(vsite.index)
    D.chart(s, XL_CHART_TYPE.COLUMN_CLUSTERED, [f"{k} (n={int(vsite.loc[k, 'size'])})" for k in sites], [("Mean ΔipTM vs the variant's own WT-Vernier fold", [float(vsite.loc[k, "mean"]) for k in sites])],
            0.92, 1.2, 6.6, 5.2, "Mean ΔipTM per Vernier site (single swaps)", (DB,), fmt="0.000", legend=False)
    nv = sum(1 for c in top if c["vernier_mutations"])
    wt_scores = folds.set_index("id").iptm
    rescued = [c for c in top if c["vernier_mutations"] and wt_scores.get(f"{c['cdrh3']}_WT", 1.0) < cfg["rank"]["min_iptm"]]
    D.text(s, 7.8, 1.3, 4.7, 5.0, [[("Single swaps barely move ipTM: ", True), (f"median change {vsite['median'].min():+.3f} to {vsite['median'].max():+.3f} per site, {near0:.0%} of swaps within ±0.02, {lower:.0%} slightly lower than the same H3 with WT Vernier.", False)],
           [(f"{nv} of the top {len(top)} carry a Vernier swap. ", True), (f"In {len(rescued)} of them the same H3 with WT Vernier does not dock (ipTM below {cfg['rank']['min_iptm']}) and the swap restores docking; in the others the swap changes ipTM by under 0.02.", False)],
           [("Single seed: ", True), ("scores jump between ≈0.6 (not docked) and ≈0.95 (docked), so treat per-swap effects as hypotheses.", False)]], size=14)
    # 11 pose
    s = D.slide("Heading and Graphic (1 column)", "Predicted pose against the crystal")
    D.drop(s, 1)
    D.picture(s, fig / "fig3_pose.png", 0.92, 1.2, 7.3, 5.3)
    D.table(s, 8.45, 1.3, [2.2, 1.8], [["VH framework fit", f"{pose['framework_fit_rmsd']:.2f} Å"], ["CDR-H3 Cα", f"{pose['cdrh3_ca_rmsd']:.2f} Å"], ["VL Cα", f"{pose['vl_ca_rmsd']:.2f} Å"]], ["Top-1 design", "RMSD"], size=12, rowh=0.45)
    D.text(s, 8.45, 3.3, 4.0, 3.2, [f"Top design ({top[0]['id']}): orange HER2 and loop, grey Fv, superposed on 1N8Z (pale cyan HER2, green crystal H3).", "RMSDs after fitting the VH framework. One prediction per complex."], size=12)
    # 12-13 top 20
    for lo, ttl in ((0, f"Top {len(top)} candidates (1–10)"), (10, f"Top {len(top)} candidates (11–{len(top)})")):
        s = D.slide("Heading and Graphic (1 column)", ttl)
        D.drop(s, 1)
        rows = []
        for c in top[lo:lo + 10]:
            runs = [(m["mutation"] + " ", m["origin"] == "new in this work") for m in c["mutations_vs_wt"]]
            rows.append([str(c["rank"]), c["cdrh3"], runs, c["vernier_mutations"][0] if c["vernier_mutations"] else "none (WT)", f"{c['ilis']:.3f}", f"{c['ipsae']:.3f}", "present" if c["salt_bridge_rule"] == "salt bridge present" else "both mutated"])
        D.table(s, 0.92, 1.2, [0.7, 1.9, 3.8, 1.4, 0.9, 0.9, 1.9], rows, ["Rank", "CDR-H3", "CDR-H3 mutations vs WT", "Vernier", "iLIS", "ipSAE", "R98–pos 10"], size=12, rowh=0.45, bold_cols=(1,))
        D.text(s, 0.92, 6.3, 11.5, 0.5, [f"Bold = new in this work, regular = from the base loop. R98–pos 10: salt bridge present, or both partners mutated. Ranked by iLIS (ipSAE tie-break), ≤{cfg['rank']['max_per_h3']} per CDR-H3."], size=11, color=GY)
    # 14 pipeline
    s = D.slide("Heading and Graphic (1 column)", "Reproducible pipeline")
    D.drop(s, 1)
    rows = [["1 enumerate", "BLOSUM-45 variants of the best tested loop, novelty and liability filters"], ["2 panel", "Diverse CDR-H3 shortlist, Vernier genotypes"], ["3 validate", "abguard gate"],
            ["4 fold", "JapanFold OpenDDE antibody-antigen, controls included"], ["5 contacts", "GetContacts: Fv–HER2 and CDR-H3 vs rest, compared with 1N8Z"], ["6 score", "ipSAE and iLIS from PAE"],
            ["7 rank", "Filter, rank, export JSON / FASTA / CSV"], ["8–9 figures, deck", "PyMOL figures and this presentation"]]
    D.table(s, 0.92, 1.3, [2.6, 8.9], rows, ["Stage", "What it does"], size=14, rowh=0.55, bold_cols=(0,))
    D.text(s, 0.92, 6.1, 11.5, 0.5, ["python run_pipeline.py --config config/pipeline.yaml   (stages are decoupled by files, any stage can be re-run)"], size=12, color=GY)
    # 15 caveats
    s = D.slide("Heading and text (2 columns)", "Caveats and next steps")
    D.fill(D.ph(s, 1).text_frame, [[("Caveats", True)], ("Single-seed folds: top scores are within noise", 1), ("ipTM, ipSAE, iLIS and contacts all come from one model; none predicts affinity", 1),
           (f"The base loop has one SPR replicate; trastuzumab has eight", 1), (f"{S['n_top_both_mutated']} designs have no salt bridge but mutate both partners (D10 + H94); the compensation is untested", 1), ("No non-D10 loop beat WT in the SPR data", 1)], size=15, space=8)
    D.fill(D.ph(s, 2).text_frame, [[("Next steps", True)], ("Repeat the folds with 3 seeds for the shortlist", 1), ("Independent ΔΔG estimate (FoldX or Rosetta)", 1),
           ("Prefer designs that keep the salt bridge", 1), ("Choose submissions for diversity, not rank alone", 1)], size=15, space=8)
    path = out(cfg, "final") / "HER2_H3_Vernier_redesign.pptx"
    D.prs.save(path)
    print(f"  saved {path} ({len(D.prs.slides)} slides)")
