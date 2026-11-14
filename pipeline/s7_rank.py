"""Stage 7: final ranking and export (JSON, FASTA, CSV).

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


def run(cfg: dict) -> None:
    r, base = cfg["rank"], cfg["design"]["base_h3"]
    meta = pd.read_csv(out(cfg, "panel_meta.csv"))
    d = prefilter(cfg).merge(pd.read_csv(out(cfg, "scores.csv")), on="id")
    ctrl = d[d.id.str.startswith("ctrl_")]
    d = d[~d.id.str.startswith("ctrl_")].merge(meta, on="id")
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
            salt_bridge_r98_pos10=bool(row.salt_bridge_r98_pos10), blosum45_logp=round(row.blosum_logp, 2),
            ridge_sum=round(row.ridge_sum, 3), vh_sequence=vh, vl_sequence=vl))
        fasta.append(f">{i:02d}_{row.id}_VH cdrh3={row.cdrh3} ilis={row.ilis:.3f}\n{vh}\n>{i:02d}_{row.id}_VL\n{vl}\n")
    json.dump({"base_h3": base, "wild_type_h3": WT_H3, "candidates": recs}, open(final / "top_candidates.json", "w"), indent=1)
    (final / "top_candidates.fasta").write_text("".join(fasta))
    with open(final / "top_mutations.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "id", "cdrh3", "mutations_vs_wt", "mutations_vs_base", "vernier", "iptm", "ipsae", "ilis",
                    "interface_recall", "h3_vs_rest_recall", "salt_bridge_r98_pos10"])
        for x in recs:
            w.writerow([x["rank"], x["id"], x["cdrh3"], " ".join(m["mutation"] for m in x["mutations_vs_wt"]),
                        " ".join(x["mutations_vs_base"]), " ".join(x["vernier_mutations"]) or "none", x["opendde_iptm"], x["ipsae"],
                        x["ilis"], x["interface_recall"], x["h3_vs_rest_recall"], x["salt_bridge_r98_pos10"]])
    ctrl.to_csv(final / "controls.csv", index=False)
    d.to_csv(final / "all_filtered_designs.csv", index=False)
    print(f"  top {len(recs)} written to {final} ({len({x['cdrh3'] for x in recs})} distinct CDR-H3, "
          f"{sum(x['salt_bridge_r98_pos10'] for x in recs)} keep the R98-pos10 salt bridge)")
    save_summary(cfg, n_top=len(recs), n_top_h3=len({x["cdrh3"] for x in recs}),
                 n_top_salt_bridge=int(sum(x["salt_bridge_r98_pos10"] for x in recs)))
