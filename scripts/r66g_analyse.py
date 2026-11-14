"""Does VL R66G move the epitope? Matched pairs, measured against seed noise.

HOW TO READ THIS
  Each entry was folded with and without R66G, by the same model, at the same
  seeds, from the same template. Two numbers come out of that:

    pair overlap   the epitope of the stabilised structure against the
                   unstabilised one, same seed. 1.0 = identical contact set.

    seed overlap   the epitope of a sequence against ITSELF at a different
                   seed. This is the noise floor. The model does not give the
                   same answer twice, and a pair overlap no worse than the
                   seed overlap is a difference we cannot distinguish from
                   the model talking to itself.

  The comparison that matters is pair overlap versus seed overlap. A pair
  overlap of 0.85 sounds like a shift until you see that the same sequence
  only agrees with itself at 0.83.

  Trastuzumab is in the panel both ways. It is the backbone we know binds, so
  the epitope shift R66G causes there is the best available estimate of what
  the mutation does on its own, uncontaminated by a redesigned loop.

USAGE
  python scripts/r66g_analyse.py --runs runs --out results/r66g
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from epitope import (  # noqa: E402
    AB_CHAINS,
    AG_CHAINS,
    AG_OFFSET,
    CUTOFF,
    contact_residues,
    jaccard,
    read_cif_atoms,
    read_pdb_atoms,
)


def find_r66g_runs(runs_dir: Path) -> dict[int, Path]:
    """seed -> predictions dir, for the r66g runs only."""
    found: dict[int, Path] = {}
    for run in sorted(runs_dir.glob("*-r66g-boltz2-*")):
        manifest, pred = run / "run.json", run / "predictions"
        if not manifest.exists() or not pred.is_dir():
            continue
        meta = json.loads(manifest.read_text())
        seed = meta.get("seed")
        if seed is not None:
            found[int(seed)] = pred
    return found


def confidence(pred_dir: Path, qid: str, seed: int) -> dict[str, float]:
    """Boltz confidence for one prediction.

    `iptm_ab_ag` is the one that matters: the interface pTM between the
    antibody chains (A, B) and the antigen chains (C, D), averaged over the
    chain-pair matrix. Overall `iptm` mixes in the VH-VL interface, which is
    not what we are asking about.
    """
    import numpy as np

    out: dict[str, float] = {}
    hits = sorted((pred_dir / qid).glob(f"seed_{seed}/*_scores.json"))
    if not hits:
        return out
    try:
        d = json.loads(hits[0].read_text())
    except Exception:
        return out
    if "interface_predicted_tm_score" in d:
        out["iptm"] = float(d["interface_predicted_tm_score"])
    if "average_predicted_local_distance_difference_test" in d:
        out["plddt"] = float(d["average_predicted_local_distance_difference_test"])
    m = d.get("chain_chain_interface_predicted_tm_score")
    if m:
        m = np.asarray(m, dtype=float)
        # chain order in the queries is A (VL), B (VH), C, D (HER2)
        pairs = [(a, g) for a in (0, 1) for g in (2, 3)]
        out["iptm_ab_ag"] = float(np.mean([(m[a, g] + m[g, a]) / 2 for a, g in pairs]))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=Path, default=Path("runs"))
    ap.add_argument("--panel", type=Path, default=Path("runs/r66g/panel.csv"))
    ap.add_argument("--crystal", type=Path, default=Path("data/1n8z.pdb"))
    ap.add_argument("--out", type=Path, default=Path("results/r66g"))
    ap.add_argument("--cutoff", type=float, default=CUTOFF)
    args = ap.parse_args(argv)

    panel = {r["id"]: r for r in csv.DictReader(args.panel.open())}
    bases = [i for i in panel if not i.endswith("-stab")]

    crystal: set[int] = set()
    if args.crystal.exists():
        crystal = contact_residues(
            read_pdb_atoms(args.crystal, ("A", "B", "C")),
            ("A", "B"), ("C",), cutoff=args.cutoff,
        )

    runs = find_r66g_runs(args.runs)
    if not runs:
        raise SystemExit(f"no r66g runs under {args.runs}")
    print(f"seeds found: {sorted(runs)}")

    # epitopes[(qid, seed)] = set of HER2 residues
    epitopes: dict[tuple[str, int], set[int]] = {}
    confs: dict[tuple[str, int], dict[str, float]] = {}
    for seed, pred in runs.items():
        for cif in sorted(pred.glob("*/seed_*/*_model.cif")):
            qid = cif.parent.parent.name
            epitopes[(qid, seed)] = contact_residues(
                read_cif_atoms(cif), AB_CHAINS, AG_CHAINS, AG_OFFSET, args.cutoff
            )
            confs[(qid, seed)] = confidence(pred, qid, seed)

    seeds = sorted(runs)

    # ---- the noise floor: same sequence, different seed ----
    seed_overlaps: list[float] = []
    noise_rows = []
    for qid in panel:
        vals = [
            jaccard(epitopes[(qid, a)], epitopes[(qid, b)])
            for a, b in combinations(seeds, 2)
            if (qid, a) in epitopes and (qid, b) in epitopes
        ]
        if vals:
            seed_overlaps += vals
            noise_rows.append({
                "id": qid,
                "n_seed_pairs": len(vals),
                "seed_overlap_median": round(statistics.median(vals), 3),
                "seed_overlap_min": round(min(vals), 3),
            })
    floor = statistics.median(seed_overlaps) if seed_overlaps else float("nan")
    floor_lo = min(seed_overlaps) if seed_overlaps else float("nan")

    # ---- the question: same seed, +/- R66G ----
    pair_rows = []
    for base in bases:
        stab = f"{base}-stab"
        for seed in seeds:
            if (base, seed) not in epitopes or (stab, seed) not in epitopes:
                continue
            ew, es = epitopes[(base, seed)], epitopes[(stab, seed)]
            cw, cs = confs.get((base, seed), {}), confs.get((stab, seed), {})
            pair_rows.append({
                "id": base,
                "seed": seed,
                "n_epitope_wt": len(ew),
                "n_epitope_stab": len(es),
                "pair_overlap": round(jaccard(ew, es), 3),
                "lost_vs_wt": len(ew - es),
                "gained_vs_wt": len(es - ew),
                "overlap_wt_vs_crystal": round(jaccard(ew, crystal), 3) if crystal else "",
                "overlap_stab_vs_crystal": round(jaccard(es, crystal), 3) if crystal else "",
                "iptm_ab_ag_wt": round(cw["iptm_ab_ag"], 4) if "iptm_ab_ag" in cw else "",
                "iptm_ab_ag_stab": round(cs["iptm_ab_ag"], 4) if "iptm_ab_ag" in cs else "",
                "plddt_wt": round(cw["plddt"], 2) if "plddt" in cw else "",
                "plddt_stab": round(cs["plddt"], 2) if "plddt" in cs else "",
            })

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "pairs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pair_rows[0]))
        w.writeheader(); w.writerows(pair_rows)
    with (args.out / "seed_noise.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(noise_rows[0]))
        w.writeheader(); w.writerows(noise_rows)

    # ---- report ----
    print(f"\ncrystal epitope (1N8Z): {len(crystal)} HER2 residues")
    print(f"\nNOISE FLOOR -- same sequence, different seed, {len(seed_overlaps)} comparisons")
    print(f"  median epitope overlap of a sequence with itself: {floor:.3f}")
    print(f"  worst: {floor_lo:.3f}")
    print("\nR66G -- same seed, with vs without the mutation")
    print(f'{"id":18s} {"seed":>5s} {"overlap":>8s} {"lost":>5s} {"gain":>5s} '
          f'{"iptmAg wt":>10s} {"iptmAg stab":>12s}')
    for r in pair_rows:
        print(f'{r["id"]:18s} {r["seed"]:>5} {r["pair_overlap"]:>8.3f} '
              f'{r["lost_vs_wt"]:>5} {r["gained_vs_wt"]:>5} '
              f'{str(r["iptm_ab_ag_wt"]):>10s} {str(r["iptm_ab_ag_stab"]):>12s}')

    # ---- the binding proxy: does R66G cost interface confidence? ----
    # Same treatment as the epitope. The seed-to-seed spread of iptm_ab_ag for
    # ONE sequence is the noise floor; the +/- R66G difference has to clear it.
    iptm = {k: v["iptm_ab_ag"] for k, v in confs.items() if "iptm_ab_ag" in v}
    seed_spreads = []
    for qid in panel:
        vals = [iptm[(qid, s)] for s in seeds if (qid, s) in iptm]
        if len(vals) >= 2:
            seed_spreads.append(max(vals) - min(vals))
    iptm_floor = statistics.median(seed_spreads) if seed_spreads else float("nan")

    print(f"\nINTERFACE CONFIDENCE (iptm antibody-antigen)")
    print(f"  seed-to-seed spread for one sequence: median {iptm_floor:.4f}, "
          f"max {max(seed_spreads):.4f}" if seed_spreads else "  no data")
    deltas_all = []
    iptm_lines = []
    for base in bases:
        ds = [
            iptm[(f"{base}-stab", s)] - iptm[(base, s)]
            for s in seeds
            if (base, s) in iptm and (f"{base}-stab", s) in iptm
        ]
        if not ds:
            continue
        deltas_all += ds
        md = statistics.median(ds)
        sign = "lower" if md < 0 else "higher"
        flag = "clears noise" if abs(md) > iptm_floor else "within noise"
        iptm_lines.append((base, md, ds, flag))
        print(f"  {base:18s} delta(stab - wt) = {md:+.4f}  "
              f"({', '.join(f'{d:+.4f}' for d in ds)})  stab is {sign}, {flag}")
    overall = statistics.median(deltas_all) if deltas_all else float("nan")
    n_neg = sum(1 for d in deltas_all if d < 0)
    print(f"  ALL PAIRS: median delta {overall:+.4f}; "
          f"{n_neg} of {len(deltas_all)} comparisons are lower with R66G")

    print("\nper design, median across seeds, against the noise floor")
    verdicts = []
    for base in bases:
        mine = [r["pair_overlap"] for r in pair_rows if r["id"] == base]
        if not mine:
            continue
        med = statistics.median(mine)
        call = "within noise" if med >= floor else "BELOW the noise floor"
        verdicts.append((base, med, call))
        print(f"  {base:18s} pair overlap {med:.3f}   vs floor {floor:.3f}   {call}")

    lines = [
        "# VL R66G and the epitope: matched pairs, Boltz-2",
        "",
        f"Seeds: {', '.join(map(str, seeds))}. "
        f"{len(bases)} matched pairs, {len(pair_rows)} paired predictions.",
        "Each pair is the same sequence with and without R66G, folded by the same",
        "model at the same seed from the same template.",
        "",
        "## The noise floor comes first",
        "",
        f"Folded twice at different seeds, the SAME sequence reproduces its own",
        f"epitope at a median Jaccard overlap of **{floor:.3f}** "
        f"(worst case {floor_lo:.3f}, {len(seed_overlaps)} comparisons).",
        "",
        "That is the resolution of this experiment. Any R66G effect smaller than",
        "it is not measurable here, whatever the pair numbers look like.",
        "",
        "## The pairs",
        "",
        "| Design | Median pair overlap | Verdict |",
        "| --- | ---: | --- |",
    ] + [f"| {b} | {m:.3f} | {c} |" for b, m, c in verdicts] + [
        "",
        f"Noise floor: {floor:.3f}.",
        "",
        "## The binding proxy: interface confidence",
        "",
        "Boltz-2's antibody-antigen interface pTM is not an affinity, and the",
        "calibration work in `results/folding/` showed it does not order variants",
        "whose KD is known. It is reported here only because the published trade",
        "says R66G costs affinity, and a model that agreed would be weak",
        "corroboration. Read the sign, not the size.",
        "",
        f"Seed-to-seed spread of this number for a single sequence: "
        f"median {iptm_floor:.4f}.",
        "",
        "| Design | Median delta (stab - wt) | Per seed | |",
        "| --- | ---: | --- | --- |",
    ] + [
        f"| {b} | {m:+.4f} | {', '.join(f'{d:+.4f}' for d in ds)} | {f} |"
        for b, m, ds, f in iptm_lines
    ] + [
        "",
        f"Across all pairs the median change is **{overall:+.4f}**, and "
        f"{n_neg} of {len(deltas_all)} comparisons are lower with R66G.",
        "",
        "Per-seed detail in `pairs.csv`, per-sequence seed reproducibility in",
        "`seed_noise.csv`.",
        "",
    ]
    (args.out / "R66G.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"\nwrote {args.out}/pairs.csv, seed_noise.csv, R66G.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
