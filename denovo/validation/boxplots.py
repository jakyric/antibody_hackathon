"""Boxplots of predicted monomeric yield (ESMC-Nb monoyield), NanoMelt Tm and OpenDDE-abag ipSAE over all designs, with the final
top-20 list marked. Run: denovo/.venv/bin/python denovo/validation/boxplots.py  -> results/denovo/boxplots_yield_tm_ipsae.png"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L = pd.read_csv("results/denovo/leaderboard.csv")
L = L[L.is_control == 0]
top = pd.read_csv("results/denovo/submission_top20_pareto.csv")
topids = list(top.design_id)
INK, MUTED, GRID, BOX, ACCENT = "#1f2933", "#6b7280", "#e5e7eb", "#9aa5b1", "#d55e00"
panels = [("my_yield", "Predicted monomeric yield\n(ESMC-Nb monoyield, log10 mg/L)", "monoyield"),
          ("nm_tm", "NanoMelt predicted Tm (°C)", "NanoMelt"),
          ("orth_dde_ipsae_chn15", "OpenDDE-abag ipSAE\n(chain-length d0, 15 Å cutoff)", "OpenDDE ipSAE")]
rng = np.random.default_rng(1)
fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.6))
for ax, (col, ylab, name) in zip(axes, panels):
    v = L[["design_id", "framework", col]].dropna()
    ax.set_facecolor("white")
    bp = ax.boxplot(v[col], positions=[0], widths=0.42, showfliers=False, patch_artist=True,
                    medianprops=dict(color=INK, linewidth=1.6), whiskerprops=dict(color=MUTED, linewidth=1.2),
                    capprops=dict(color=MUTED, linewidth=1.2), boxprops=dict(facecolor="#e8ecf0", edgecolor=MUTED, linewidth=1.2))
    ax.scatter(rng.normal(0, 0.09, len(v)), v[col], s=9, color=BOX, alpha=0.35, linewidths=0, zorder=2)
    m = v[v.design_id.isin(topids)]
    ax.scatter(rng.normal(0, 0.05, len(m)), m[col], s=70, marker="D", color=ACCENT, edgecolor="white", linewidths=1.6, zorder=4)
    span = v[col].max() - v[col].min()
    labs, last = [], None
    for _, r in m.sort_values(col, ascending=False).head(3).iterrows():
        y = r[col] if last is None else min(r[col], last - 0.045 * span)   # keep labels 4.5% of the axis range apart
        labs.append((r.design_id.replace("-rfd4jf", ""), r[col], y)); last = y
    for name_, yv, yl in labs:
        ax.annotate(name_, (0.05, yv), xytext=(0.30, yl), fontsize=8.5, color=INK, va="center",
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_xlim(-0.55, 0.95); ax.set_xticks([0]); ax.set_xticklabels([f"all designs\nn = {len(v)}   (top-20 with a value: {len(m)})"], fontsize=9, color=MUTED)
    ax.set_ylabel(ylab, fontsize=10, color=INK)
    ax.set_title(name, fontsize=12, color=INK, loc="left", fontweight="bold")
    ax.yaxis.grid(True, color=GRID, linewidth=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
h1 = plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#e8ecf0", markeredgecolor=MUTED, markersize=11, label="all designs (box = quartiles, dots = designs)")
h2 = plt.Line2D([0], [0], marker="D", color="w", markerfacecolor=ACCENT, markeredgecolor="white", markersize=9, label="final top-20 list (Pareto selection)")
fig.legend(handles=[h1, h2], loc="lower center", ncol=2, frameon=False, fontsize=10, labelcolor=INK)
fig.suptitle("Developability and binding scores of all RFD4 designs, with the selected top 20 marked", fontsize=13, color=INK, x=0.01, ha="left", fontweight="bold")
fig.text(0.01, 0.915, "Predictions, not measurements. Coverage differs per metric: not every design has every score yet.", fontsize=9.5, color=MUTED, ha="left")
fig.tight_layout(rect=(0, 0.06, 1, 0.90))
fig.savefig("results/denovo/boxplots_yield_tm_ipsae.png", dpi=150)
for col, _, name in panels:
    v = L[["design_id", col]].dropna(); t = v[v.design_id.isin(topids)]
    print(f"{name}: n={len(v)} median={v[col].median():.2f} q1={v[col].quantile(.25):.2f} q3={v[col].quantile(.75):.2f} | top-20 with value: {len(t)}, their median {t[col].median():.2f}" if len(t) else f"{name}: n={len(v)} median={v[col].median():.2f} | no top-20 design has a value")
