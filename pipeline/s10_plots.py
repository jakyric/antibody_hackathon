"""Stage 10 (optional): score distributions and the ipSAE-iLIS plane, with trastuzumab marked.

fig_distributions.png : iLIS, ipSAE, epitope / HER2 RMSD to the crystal and the number of residue-pair contacts (Fv-HER2, CDR-H3 vs rest) for every
                        folded design; designs that pass the contact + salt-bridge filters in colour, final list as ticks.
fig_ipsae_ilis.png    : ipSAE vs iLIS, colour = Fv-HER2 contacts, final designs outlined.
Trastuzumab is shown as the *prediction* of its sequence (ctrl_wt, same model and pipeline as the designs) and, for the
contact counts, also as the crystal value (1N8Z). The base loop and the two measured non-binders are shown for reference.
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .common import out, read_summary  # noqa: E402

BLUE, TEAL, GREEN, GREY, RED, ORANGE = "#004B6F", "#009299", "#63A44B", "#B8C2C7", "#D62728", "#F28C1A"


def load(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, set, set]:
    d = (pd.read_csv(out(cfg, "fold_scores.csv")).merge(pd.read_csv(out(cfg, "contacts.csv")), on="id")
         .merge(pd.read_csv(out(cfg, "scores.csv")), on="id").merge(pd.read_csv(out(cfg, "pose.csv")), on="id"))
    top = json.load(open(out(cfg, "final") / "top_candidates.json"))["candidates"]
    final = {c["id"] for c in top}
    ok = set(pd.read_csv(out(cfg, "final") / "all_filtered_designs.csv").query("salt_bridge_rule != 'excluded'").id) \
        if "salt_bridge_rule" in pd.read_csv(out(cfg, "final") / "all_filtered_designs.csv").columns else set()
    return d[~d.id.str.startswith("ctrl_")], d[d.id.str.startswith("ctrl_")].set_index("id"), final, ok


def run(cfg: dict) -> None:
    d, ctrl, final, kept = load(cfg)
    S = read_summary(cfg)
    figs = out(cfg, "figures")
    wt, base = ctrl.loc["ctrl_wt"], ctrl.loc["ctrl_base"]
    sel = d[d.id.isin(kept)]
    panels = [("ilis", "iLIS (Fv–HER2)", np.linspace(0, 0.85, 35), None), ("ipsae", "ipSAE (Fv–HER2)", np.linspace(0, 1, 35), None),
              ("epitope_rmsd", "HER2 epitope Cα RMSD to the crystal (Å)", np.arange(0, 38, 1), None),
              ("if_pairs", "Fv–HER2 contacts (residue pairs)", np.arange(0, 80, 2), S.get("crystal_if_pairs")),
              ("in_pairs", "CDR-H3 vs rest contacts (residue pairs)", np.arange(0, 60, 2), S.get("crystal_intra_pairs")),
              ("her2_rmsd", "HER2 Cα RMSD to the crystal (Å)", np.arange(0, 45, 1.5), None)]
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    for ax, (col, lab, bins, crystal) in zip(axes.ravel(), panels):
        ax.hist(d[col], bins=bins, color=GREY, label=f"all folded designs (n={len(d)})")
        if len(sel):
            ax.hist(sel[col], bins=bins, color=TEAL, label=f"pass contact + salt-bridge filters (n={len(sel)})")
        ymax = ax.get_ylim()[1]
        for x in d[d.id.isin(final)][col]:
            ax.plot([x, x], [0, ymax * 0.06], color=BLUE, lw=1.4)
        ax.axvline(wt[col], color=RED, lw=2.2, label=f"trastuzumab, predicted ({wt[col]:.2f})" if col in ("ilis", "ipsae", "epitope_rmsd", "her2_rmsd") else f"trastuzumab, predicted ({wt[col]:.0f})")
        ax.axvline(base[col], color=GREEN, lw=2, ls="--", label=f"best tested loop, predicted ({base[col]:.2f})" if col in ("ilis", "ipsae", "epitope_rmsd", "her2_rmsd") else f"best tested loop, predicted ({base[col]:.0f})")
        if col == "epitope_rmsd":
            ax.axvline(cfg["rank"]["max_epitope_rmsd"], color="black", lw=1, ls="-.", label=f"filter ({cfg['rank']['max_epitope_rmsd']:g} Å)")
        if crystal:
            ax.axvline(crystal, color="black", lw=2, ls=":", label=f"trastuzumab crystal 1N8Z ({crystal})")
        ax.set_xlabel(lab)
        ax.set_ylabel("designs")
        ax.legend(fontsize=8, frameon=True, framealpha=0.9, edgecolor='none')
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle("Score distributions of the folded designs; ticks along the x-axis = final list", fontsize=13)
    fig.tight_layout()
    fig.savefig(figs / "fig_distributions.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    sc = ax.scatter(d.ipsae, d.ilis, c=d.if_pairs, cmap="viridis", s=26, alpha=0.85, edgecolor="none")
    fin = d[d.id.isin(final)]
    ax.scatter(fin.ipsae, fin.ilis, facecolor="none", edgecolor=ORANGE, s=95, lw=1.8, label=f"final list (n={len(fin)})")
    for name, row, mk, col, lab in (("ctrl_wt", wt, "*", RED, "trastuzumab (predicted)"), ("ctrl_base", base, "s", GREEN, "best tested loop"),
                                    ("ctrl_nonbinder1", ctrl.loc["ctrl_nonbinder1"], "X", "black", "measured non-binders"),
                                    ("ctrl_nonbinder2", ctrl.loc["ctrl_nonbinder2"], "X", "black", None)):
        ax.scatter(row.ipsae, row.ilis, marker=mk, s=230 if mk == "*" else 90, color=col, edgecolor="white", lw=0.8, zorder=5, label=lab)
    ax.annotate(f"trastuzumab\n{int(wt.if_pairs)} pairs (crystal: {S.get('crystal_if_pairs', '?')})", (wt.ipsae, wt.ilis), xytext=(0.12, 0.30),
                arrowprops=dict(arrowstyle="->", color=RED), color=RED, fontsize=9)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Fv–HER2 contacts (residue pairs)")
    ax.set_xlabel("ipSAE (Fv–HER2)")
    ax.set_ylabel("iLIS (Fv–HER2)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(figs / "fig_ipsae_ilis.png", dpi=200)
    plt.close(fig)
    print(f"  wrote fig_distributions.png and fig_ipsae_ilis.png to {figs}")
