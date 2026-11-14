#!/usr/bin/env python3
"""Read the scores from the fold runs and answer one question: do they rank affinity?

We know KD for the calibration part of the panel. A fold model that cannot
rank the measured variants cannot rank our designs either. That is measured
here, per model, instead of assumed.

    python scripts/fold_analyse.py --runs runs --out results/folding
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.metrics import auroc, bootstrap_ci, spearman  # noqa: E402

# chain order in queries.json: A=VL, B=VH, C and D = HER2 domain IV
AB_CHAINS = (0, 1)
AG_CHAINS = (2, 3)


def interface_iptm(score: dict) -> dict[str, float]:
    """ipTM split up: antibody against antigen, and VH against antigen."""
    m = score.get("chain_chain_interface_predicted_tm_score")
    out: dict[str, float] = {}
    if not m:
        return out
    m = np.asarray(m, dtype=float)
    pairs = [(a, g) for a in AB_CHAINS for g in AG_CHAINS]
    out["iptm_ab_ag"] = float(np.mean([(m[a, g] + m[g, a]) / 2 for a, g in pairs]))
    vh_pairs = [(1, g) for g in AG_CHAINS]
    out["iptm_vh_ag"] = float(np.mean([(m[a, g] + m[g, a]) / 2 for a, g in vh_pairs]))
    return out


def _prediction_files(run_dir: Path) -> list[tuple[str, Path, int | None]]:
    """(query_id, score file, seed). Uses result.json when it exists.

    A run that was interrupted has no result.json, but the score files are
    there. They should be readable, otherwise half a run becomes no run.
    """
    manifest = run_dir / "result.json"
    if manifest.exists():
        res = json.loads(manifest.read_text())
        return [
            (p["query_id"], run_dir / p["scores_path"], p.get("seed"))
            for p in res["result"]["predictions"]
        ]
    out = []
    for sp in sorted(run_dir.glob("*/seed_*/*_scores.json")):
        out.append((sp.parent.parent.name, sp, int(sp.parent.name.split("_")[-1])))
    return out


def read_predictions(run_dir: Path) -> dict[str, dict]:
    """query_id -> scores, for one model run."""
    rows: dict[str, dict] = {}
    for query_id, sp, seed in _prediction_files(run_dir):
        if not sp.exists():
            continue
        s = json.loads(sp.read_text())
        row = {
            "plddt": s.get("average_predicted_local_distance_difference_test"),
            "ptm": s.get("predicted_tm_score"),
            "iptm": s.get("interface_predicted_tm_score"),
            "pde": s.get("global_predicted_distance_error"),
            "ranking_score": s.get("model_ranking_score"),
            "seed": seed,
        }
        row.update(interface_iptm(s))
        rows.setdefault(query_id, row)
    return rows


def panel_of(manifest: dict) -> str | None:
    """Which panel the run used, read from what the run itself recorded.

    WHY THIS MUST BE READ AND NOT ASSUMED
      Two panels number the variants the same way — `calib-05` exists in
      both — but they are different sequences with different KD. If we read
      the scores from one run against the KDs from another panel, we get a
      number that looks real and is not. That happened. Now the panel
      follows the run.
    """
    for key in ("inputs", "input"):
        inputs = manifest.get(key)
        if isinstance(inputs, dict):
            names = list(inputs)
        elif isinstance(inputs, list):
            names = [str(x) for x in inputs]
        else:
            continue
        for name in names:
            if name.endswith("panel.csv"):
                return name
    return None


def find_runs(runs_dir: Path) -> dict[str, tuple[Path, str | None]]:
    """model name -> (the folder with predictions, the panel the run used)."""
    found: dict[str, tuple[Path, str | None]] = {}
    for run in sorted(runs_dir.glob("*-fold-*")):
        manifest_path = run / "run.json"
        pred = run / "predictions"
        if not manifest_path.exists() or not _prediction_files(pred):
            continue
        manifest = json.loads(manifest_path.read_text())
        key = manifest.get("name", run.name).replace("fold-", "", 1)
        found[key] = (pred, panel_of(manifest))
    return found


def load_panel(path: Path) -> dict[str, dict]:
    with path.open(newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


METRICS = ["iptm_ab_ag", "iptm_vh_ag", "iptm", "ptm", "plddt", "ranking_score"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, default=Path("runs"))
    ap.add_argument("--panel", type=Path, default=Path("runs/panel/panel.csv"))
    ap.add_argument("--out", type=Path, default=Path("results/folding"))
    args = ap.parse_args()

    runs = find_runs(args.runs)
    if not runs:
        raise SystemExit(f"found no fold runs in {args.runs}")

    # Each run is read against its own panel. If we fall back to --panel for
    # a run that used a different one, we say so: the ids collide between
    # panels, so a silent mix-up gives a number that means nothing.
    default_panel = load_panel(args.panel)
    panels: dict[str, dict] = {}
    for model, (_, panel_rel) in runs.items():
        if panel_rel and (Path(panel_rel).exists() or (args.runs.parent / panel_rel).exists()):
            path = Path(panel_rel) if Path(panel_rel).exists() else args.runs.parent / panel_rel
            panels[model] = load_panel(path)
        else:
            if panel_rel and panel_rel != str(args.panel):
                print(f"  WARNING: {model} used {panel_rel}, which does not exist. Skipping.")
                continue
            panels[model] = default_panel
    runs = {m: v for m, v in runs.items() if m in panels}

    args.out.mkdir(parents=True, exist_ok=True)

    # ── one long table with everything ───────────────────────────────────────
    long_rows = []
    per_model: dict[str, dict[str, dict]] = {}
    for model, (pred_dir, _) in sorted(runs.items()):
        rows = read_predictions(pred_dir)
        per_model[model] = rows
        panel = panels[model]
        for qid, r in rows.items():
            p = panel.get(qid, {})
            long_rows.append(
                {"model": model, "id": qid, "group": p.get("group", "?"),
                 "cdrh3": p.get("cdrh3", ""), "kd_nM": p.get("kd_nM", ""), **r}
            )
    cols = ["model", "id", "group", "cdrh3", "kd_nM", "seed"] + METRICS + ["pde"]
    with (args.out / "scores.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(long_rows)

    # ── question 1: does the model rank measured affinity? ───────────────────
    calib_lines = []
    for model, rows in sorted(per_model.items()):
        panel = panels[model]
        ids = [q for q in rows if panel.get(q, {}).get("group") in ("calibration", "baseline")]
        kd = np.array([float(panel[q]["kd_nM"]) for q in ids])
        ok = np.isfinite(kd)
        for metric in METRICS:
            vals = np.array([rows[q].get(metric) if rows[q].get(metric) is not None
                             else np.nan for q in ids], dtype=float)
            m = ok & np.isfinite(vals)
            if m.sum() < 8:
                continue
            # a high score should mean low KD, so a negative Spearman against log KD
            rho = spearman(vals[m], -np.log10(kd[m]))
            lo, hi = bootstrap_ci(
                lambda a, b: spearman(a, -np.log10(b)), vals[m], kd[m], seed=0
            )
            # can it separate the strong from the weak? (better than median KD)
            strong = (kd[m] < np.median(kd[m])).astype(float)
            au = auroc(vals[m], strong)
            calib_lines.append(
                {"model": model, "metric": metric, "n": int(m.sum()),
                 "spearman_vs_logKD": round(rho, 3),
                 "ci_low": round(lo, 3), "ci_high": round(hi, 3),
                 "auroc_strong_half": round(au, 3)}
            )
    with (args.out / "calibration.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(calib_lines[0]) if calib_lines else
                           ["model", "metric", "n", "spearman_vs_logKD", "ci_low", "ci_high",
                            "auroc_strong_half"])
        w.writeheader()
        w.writerows(calib_lines)

    # ── the same variants, ranked with what we had from before ───────────────
    # Without this we compare fold scores on one sample with edit distance
    # on another, and then the difference means nothing.
    from abx.data import TZ_CDRH3_10  # noqa: PLC0415
    from abx.gate import hamming  # noqa: PLC0415

    panel = default_panel
    calib_ids = [q for q in panel if panel[q]["group"] in ("calibration", "baseline")]
    kd = np.array([float(panel[q]["kd_nM"]) for q in calib_ids])
    ok = np.isfinite(kd)
    baselines = []
    dist = np.array([hamming(panel[q]["cdrh3"], TZ_CDRH3_10) for q in calib_ids], dtype=float)
    cand = {"edit_distance": -dist}  # close to trastuzumab = high score
    model_npz = Path("results/binder_model.npz")
    if model_npz.exists():
        try:
            from abx.model import load as load_model  # noqa: PLC0415
            from abx.model import onehot  # noqa: PLC0415
            net = load_model(model_npz)
            cand["binder_net"] = net.logits(
                onehot(np.array([panel[q]["cdrh3"] for q in calib_ids]))
            )
        except Exception as exc:  # noqa: BLE001
            print(f"  (the binder model could not be used: {exc})")
    for name, vals in cand.items():
        vals = np.asarray(vals, dtype=float)
        m = ok & np.isfinite(vals)
        rho = spearman(vals[m], -np.log10(kd[m]))
        lo, hi = bootstrap_ci(lambda a, b: spearman(a, -np.log10(b)), vals[m], kd[m], seed=0)
        baselines.append({"model": f"baseline:{name}", "metric": name, "n": int(m.sum()),
                          "spearman_vs_logKD": round(rho, 3), "ci_low": round(lo, 3),
                          "ci_high": round(hi, 3),
                          "auroc_strong_half": round(
                              auroc(vals[m], kd[m] < np.median(kd[m])), 3)})
    calib_lines += baselines
    with (args.out / "calibration.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "metric", "n", "spearman_vs_logKD",
                                           "ci_low", "ci_high", "auroc_strong_half"])
        w.writeheader()
        w.writerows(calib_lines)

    # ── question 2: do the models agree on the designs? ──────────────────────
    agree = []
    models = sorted(per_model)
    design_ids = [q for q in panel if panel[q]["group"] == "design"]
    for i, a in enumerate(models):
        for b in models[i + 1:]:
            shared = [q for q in design_ids if q in per_model[a] and q in per_model[b]]
            if len(shared) < 8:
                continue
            for metric in ("iptm_ab_ag", "plddt"):
                va = np.array([per_model[a][q].get(metric, np.nan) for q in shared], dtype=float)
                vb = np.array([per_model[b][q].get(metric, np.nan) for q in shared], dtype=float)
                m = np.isfinite(va) & np.isfinite(vb)
                if m.sum() < 8:
                    continue
                agree.append({"model_a": a, "model_b": b, "metric": metric,
                              "n": int(m.sum()), "spearman": round(spearman(va[m], vb[m]), 3)})
    with (args.out / "agreement.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model_a", "model_b", "metric", "n", "spearman"])
        w.writeheader()
        w.writerows(agree)

    # ── question 3: how do the designs stand against trastuzumab? ────────────
    base = []
    for model, rows in sorted(per_model.items()):
        tz = rows.get("trastuzumab")
        if not tz:
            continue
        for metric in ("iptm_ab_ag", "iptm", "plddt"):
            t = tz.get(metric)
            vals = np.array([rows[q].get(metric, np.nan) for q in design_ids
                             if q in rows], dtype=float)
            vals = vals[np.isfinite(vals)]
            if t is None or not len(vals):
                continue
            base.append({"model": model, "metric": metric, "trastuzumab": round(float(t), 4),
                         "designs_median": round(float(np.median(vals)), 4),
                         "designs_max": round(float(np.max(vals)), 4),
                         "n_above_trastuzumab": int((vals > t).sum()), "n_designs": len(vals)})
    with (args.out / "vs_trastuzumab.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "metric", "trastuzumab", "designs_median",
                                           "designs_max", "n_above_trastuzumab", "n_designs"])
        w.writeheader()
        w.writerows(base)

    print(f"models: {', '.join(models)}")
    print(f"wrote {args.out}/scores.csv, calibration.csv, agreement.csv, vs_trastuzumab.csv")
    for r in calib_lines:
        if r["metric"] == "iptm_ab_ag":
            print(f"  {r['model']:18s} ipTM(ab-ag) vs log KD: rho={r['spearman_vs_logKD']:+.3f} "
                  f"[{r['ci_low']:+.2f},{r['ci_high']:+.2f}]  n={r['n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
