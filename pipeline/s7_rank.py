"""Stage 7: final ranking and export (JSON, FASTA, CSV).

Salt-bridge rule (hard filter): a design is kept if the R98 (H94) - position 10 salt bridge is present in its fold, or, if
it is absent, only when both partners are mutated; one mutated partner without the bridge excludes the design.

Filters before ranking: docked (ipTM), same epitope (C-alpha RMSD of the crystal HER2 epitope after fitting the Fv
framework), crystal contacts kept, salt-bridge rule.

Ranking: iLIS (it is the only metric that spreads designs ipTM sees as tied), ipSAE as tie-break, at most `max_per_h3`
designs per CDR-H3 so the list is not one loop with many Vernier swaps. Mutations are listed against wild type
(trastuzumab) and against the best tested loop; Vernier changes are listed separately.
"""
from __future__ import annotations

import csv
import json

import pandas as pd

from .common import WT_H3, mutations, out, save_summary
from .s6_score import prefilter


def apply_salt_bridge_rule(bridge: pd.Series, pos10_mutated: pd.Series, h94_mutated: pd.Series) -> pd.Series:
    """'salt bridge present' | 'both partners mutated' | 'excluded' (bridge absent and fewer than both partners mutated)."""
    out_ = pd.Series("excluded", index=bridge.index)
    out_[bridge] = "salt bridge present"
    out_[~bridge & pos10_mutated & h94_mutated] = "both partners mutated"
    return out_


def write_filter_report(cfg: dict, selected: set, path) -> None:
    """One row per folded design: which criteria it passes and why it was not selected."""
    r = cfg["rank"]
    a = (pd.read_csv(out(cfg, "fold_scores.csv")).merge(pd.read_csv(out(cfg, "contacts.csv")), on="id")
         .merge(pd.read_csv(out(cfg, "pose.csv")), on="id").merge(pd.read_csv(out(cfg, "scores.csv")), on="id", how="left")
         .merge(pd.read_csv(out(cfg, "panel_meta.csv"))[["id", "cdrh3", "vernier_changes"]], on="id"))
    a["salt_bridge_rule"] = apply_salt_bridge_rule(a.salt_bridge_r98_pos10, a.cdrh3.str[9] != WT_H3[9], a.vernier_changes.str.startswith("H94"))
    checks = {"not docked (ipTM)": a.iptm >= r["min_iptm"], "off the crystal epitope": a.epitope_rmsd <= r["max_epitope_rmsd"],
              "Fv-HER2 contacts lost": a.if_recall >= r["min_interface_recall"], "CDR-H3-HER2 pairs lost": a.if_h3_missing <= r["max_h3_pairs_missing"],
              "CDR-H3-vs-rest contacts lost": (a.in_recall >= r["min_h3_rest_recall"]) & (a.in_polar_recall >= r["min_h3_rest_polar_recall"]),
              "salt-bridge rule": a.salt_bridge_rule != "excluded"}
    a["selected"] = a.id.isin(selected)
    a["dropped_because"] = [", ".join(k for k, v in checks.items() if not v.loc[i]) or ("per-CDR-H3 cap or rank" if not sel else "") for i, sel in zip(a.index, a.selected)]
    cols = ["id", "cdrh3", "vernier_changes", "selected", "dropped_because", "iptm", "ipsae", "ilis", "epitope_rmsd", "if_recall", "if_h3_missing", "in_recall", "in_polar_recall", "salt_bridge_rule"]
    a[cols].sort_values(["selected", "ilis"], ascending=False).to_csv(path, index=False)


def run(cfg: dict) -> None:
    r, base = cfg["rank"], cfg["design"]["base_h3"]
    meta = pd.read_csv(out(cfg, "panel_meta.csv"))
    d = prefilter(cfg).merge(pd.read_csv(out(cfg, "scores.csv")), on="id")
    ctrl = d[d.id.str.startswith("ctrl_")]
    d = d[~d.id.str.startswith("ctrl_")].merge(meta, on="id")
    n_contacts_ok = len(d)
    # Salt-bridge rule (R98/H94 - position 10): keep a design if the salt bridge is present in its fold; if it is absent,
    # keep it only when BOTH partners are mutated (H94 swap and a changed position 10); exclude it when just one partner is.
    d["pos10_mutated"] = d.cdrh3.str[9] != WT_H3[9]
    d["h94_mutated"] = d.vernier_changes.str.startswith("H94")
    d["salt_bridge_rule"] = apply_salt_bridge_rule(d.salt_bridge_r98_pos10, d.pos10_mutated, d.h94_mutated)
    if r.get("salt_bridge_rule", True):
        d = d[d.salt_bridge_rule != "excluded"]
    print(f"  salt-bridge rule: {len(d)} of {n_contacts_ok} contact-filtered designs kept "
          f"({int((d.salt_bridge_rule == 'salt bridge present').sum())} with the bridge, {int((d.salt_bridge_rule == 'both partners mutated').sum())} with both partners mutated)")
    d = d.sort_values(["ilis", "ipsae"], ascending=False)
    chosen, count = [], {}
    for row in d.itertuples():
        if count.get(row.cdrh3, 0) >= r["max_per_h3"]:
            continue
        chosen.append(row)
        count[row.cdrh3] = count.get(row.cdrh3, 0) + 1
        if len(chosen) == r["top_n"]:
            break
    panel = pd.read_csv(out(cfg, "panel_valid.csv")).set_index("sequence_id")
    final = out(cfg, "final")
    final.mkdir(exist_ok=True)
    recs, fasta = [], []
    for i, row in enumerate(chosen, 1):
        vh, vl = panel.loc[row.id, "heavy_sequence"], panel.loc[row.id, "light_sequence"]
        h3m = [{"mutation": m, "origin": "best tested loop" if row.cdrh3[int(m[1:-1]) - 1] == base[int(m[1:-1]) - 1] else "new in this work"}
               for m in mutations(row.cdrh3, WT_H3)]
        recs.append(dict(
            rank=i, id=row.id, cdrh3=row.cdrh3, mutations_vs_wt=h3m, mutations_vs_base=mutations(row.cdrh3, base),
            vernier_mutations=[] if row.vernier_changes == "WT" else [row.vernier_changes],
            opendde_iptm=round(row.iptm, 3), ipsae=round(row.ipsae, 3), ilis=round(row.ilis, 3),
            interface_recall=round(row.if_recall, 2), h3_her2_pairs_missing=int(row.if_h3_missing),
            h3_vs_rest_recall=round(row.in_recall, 2), h3_vs_rest_polar_recall=round(row.in_polar_recall, 2),
            epitope_rmsd=round(row.epitope_rmsd, 2), her2_rmsd=round(row.her2_rmsd, 2), h3_rmsd=round(row.h3_rmsd, 2),
            salt_bridge_r98_pos10=bool(row.salt_bridge_r98_pos10), salt_bridge_rule=row.salt_bridge_rule, blosum45_logp=round(row.blosum_logp, 2),
            ridge_sum=round(row.ridge_sum, 3), vh_sequence=vh, vl_sequence=vl))
        fasta.append(f">{i:02d}_{row.id}_VH cdrh3={row.cdrh3} ilis={row.ilis:.3f}\n{vh}\n>{i:02d}_{row.id}_VL\n{vl}\n")
    json.dump({"base_h3": base, "wild_type_h3": WT_H3, "candidates": recs}, open(final / "top_candidates.json", "w"), indent=1)
    (final / "top_candidates.fasta").write_text("".join(fasta))
    with open(final / "top_mutations.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "id", "cdrh3", "mutations_vs_wt", "mutations_vs_base", "vernier", "iptm", "ipsae", "ilis",
                    "interface_recall", "h3_vs_rest_recall", "epitope_rmsd", "her2_rmsd", "salt_bridge_r98_pos10", "salt_bridge_rule"])
        for x in recs:
            w.writerow([x["rank"], x["id"], x["cdrh3"], " ".join(m["mutation"] for m in x["mutations_vs_wt"]),
                        " ".join(x["mutations_vs_base"]), " ".join(x["vernier_mutations"]) or "none", x["opendde_iptm"], x["ipsae"],
                        x["ilis"], x["interface_recall"], x["h3_vs_rest_recall"], x["epitope_rmsd"], x["her2_rmsd"], x["salt_bridge_r98_pos10"], x["salt_bridge_rule"]])
    write_filter_report(cfg, {x["id"] for x in recs}, final / "filter_report.csv")
    ctrl.to_csv(final / "controls.csv", index=False)
    d.to_csv(final / "all_filtered_designs.csv", index=False)
    print(f"  {len(recs)} designs written to {final} ({len({x['cdrh3'] for x in recs})} distinct CDR-H3; "
          f"{sum(x['salt_bridge_r98_pos10'] for x in recs)} keep the R98-pos10 salt bridge, "
          f"{sum(x['salt_bridge_rule'] == 'both partners mutated' for x in recs)} have both partners mutated)")
    save_summary(cfg, n_after_salt_rule=int(len(d)), n_after_salt_rule_h3=int(d.cdrh3.nunique()), n_top=len(recs), n_top_h3=len({x["cdrh3"] for x in recs}),
                 n_top_salt_bridge=int(sum(x["salt_bridge_r98_pos10"] for x in recs)),
                 n_top_both_mutated=int(sum(x["salt_bridge_rule"] == "both partners mutated" for x in recs)))
