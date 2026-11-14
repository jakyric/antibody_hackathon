#!/usr/bin/env python3
"""Rank the de novo designs on what we have actually measured, and say how much that ranking is worth.

THREE AXES, as the lab measures them:
  binding    interface ipTM (antibody against HER2) from each folding model
  epitope    overlap between the HER2 residues the design takes and trastuzumab's own
  stability  ddG from a stability model, as change FROM trastuzumab

TWO MODELS, and they are not mixed: each number is turned into a percentile
within its own model run, and the percentiles are averaged. A model with a
template sees the HER2 structure; one without a template does not. That they
agree is evidence; that one of them is confident is not.

CALIBRATION. The same model, template and seed were run on 24 variants with
MEASURED KD. The Spearman between the score and log KD there says how much the
ranking is worth. If it is near zero, that is said outright instead of being
hidden behind a sort.

    python scripts/denovo_rank.py --panel runs/denovo-panel/panel.csv \
        --fold-run boltz2=runs/<mine>/predictions --calib-run boltz2=runs/<theirs>/predictions \
        --calib-panel runs/panel/panel.csv --epitope boltz2=results/denovo/epitope_boltz2.csv \
        --out results/denovo
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.metrics import spearman  # noqa: E402

_spec = importlib.util.spec_from_file_location("fold_analyse", ROOT / "scripts" / "fold_analyse.py")
fold_analyse = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fold_analyse)

SCORE_KEY = "iptm_ab_ag"      # interface ipTM, antibody against antigen


def pairs(arg: list[str]) -> dict[str, Path]:
    out = {}
    for item in arg:
        name, _, path = item.partition("=")
        if not path:
            raise SystemExit(f"write model=path, not {item!r}")
        out[name] = Path(path)
    return out


def percentiles(values: dict[str, float]) -> dict[str, float]:
    """Percentile within one model run. Higher score = higher percentile."""
    items = [(k, v) for k, v in values.items() if v is not None and not math.isnan(v)]
    order = sorted(items, key=lambda kv: kv[1])
    n = len(order)
    return {k: round(100 * (i + (1 if n == 1 else 0)) / max(n - 1, 1), 1)
            for i, (k, _) in enumerate(order)} if n else {}


def read_panel(path: Path) -> dict[str, dict]:
    with path.open(newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


def calibration(calib_scores: dict[str, dict], calib_panel: dict[str, dict]) -> dict:
    """Spearman between the score and log10(KD) on the variants with measured KD.

    A negative Spearman is what we hope for: low KD (strong binder) should give a high score.
    """
    xs, ys = [], []
    for qid, row in calib_scores.items():
        rec = calib_panel.get(qid)
        if not rec:
            continue
        try:
            kd = float(rec.get("kd_nM", "nan"))
        except ValueError:
            continue
        v = row.get(SCORE_KEY)
        if v is None or math.isnan(kd) or kd <= 0:
            continue
        xs.append(v)
        ys.append(math.log10(kd))
    if len(xs) < 5:
        return {"n": len(xs), "spearman_score_vs_logKD": float("nan")}
    return {"n": len(xs), "spearman_score_vs_logKD": round(float(spearman(np.array(xs), np.array(ys))), 3)}


def read_epitope(path: Path) -> dict[str, dict]:
    with path.open(newline="") as fh:
        return {r["query_id"]: r for r in csv.DictReader(fh)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, required=True)
    ap.add_argument("--fold-run", action="append", default=[], metavar="MODEL=DIR")
    ap.add_argument("--calib-run", action="append", default=[], metavar="MODEL=DIR")
    ap.add_argument("--calib-panel", type=Path)
    ap.add_argument("--epitope", action="append", default=[], metavar="MODEL=CSV")
    ap.add_argument("--thermo", type=Path,
                    help="CSV with design_id and ddg_vs_trastuzumab (lower = more stable)")
    ap.add_argument("--epitope-model", default="",
                    help="which model the epitope requirement is read from; empty = the one without a template")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    panel = read_panel(args.panel)
    folds = pairs(args.fold_run)
    if not folds:
        raise SystemExit("give at least one --fold-run model=path")
    epi = {m: read_epitope(p) for m, p in pairs(args.epitope).items()}

    # ── scores per model, and percentile within the model ────────────────
    scores: dict[str, dict[str, dict]] = {m: fold_analyse.read_predictions(p) for m, p in folds.items()}
    for model, rows in scores.items():
        if not rows:
            raise SystemExit(f"no predictions for {model} in {folds[model]}")
    pct = {m: percentiles({q: r.get(SCORE_KEY) for q, r in rows.items()})
           for m, rows in scores.items()}

    # ── calibration against measured KD, per model ───────────────────────
    calib_panel = read_panel(args.calib_panel) if args.calib_panel else {}
    calib = {}
    for model, path in pairs(args.calib_run).items():
        calib[model] = calibration(fold_analyse.read_predictions(path), calib_panel)
    # the anchors in our own batch tell whether the two batches can be compared
    for model, rows in scores.items():
        own = {q: r for q, r in rows.items() if panel.get(q, {}).get("group") in ("anchor", "baseline")}
        own_cal = calibration(own, panel)
        calib.setdefault(model, {})["own_batch_anchors"] = len(own)
        calib[model]["own_batch_spearman"] = own_cal["spearman_score_vs_logKD"]

    thermo = {}
    if args.thermo and args.thermo.exists():
        with args.thermo.open(newline="") as fh:
            thermo = {r["design_id"]: r for r in csv.DictReader(fh)}

    epi_model = args.epitope_model or (sorted(epi)[0] if epi else "")

    rows = []
    for qid, rec in panel.items():
        row = {"design_id": qid, "group": rec.get("group", ""), "cdrh3": rec.get("cdrh3", ""),
               "kd_nM_measured": rec.get("kd_nM", "")}
        pcts = []
        for model in sorted(scores):
            s = scores[model].get(qid, {})
            row[f"{model}_{SCORE_KEY}"] = round(s[SCORE_KEY], 4) if s.get(SCORE_KEY) is not None else ""
            row[f"{model}_plddt"] = round(s["plddt"], 2) if s.get("plddt") is not None else ""
            p = pct[model].get(qid)
            row[f"{model}_percentile"] = p if p is not None else ""
            if p is not None:
                pcts.append(p)
            if model in epi:
                e = epi[model].get(qid, {})
                row[f"{model}_epitope_jaccard"] = e.get("jaccard_to_predicted_trastuzumab", "")
        row["models_scored"] = len(pcts)
        row["binding_percentile_mean"] = round(float(np.mean(pcts)), 1) if pcts else ""
        row["binding_percentile_spread"] = round(float(max(pcts) - min(pcts)), 1) if len(pcts) > 1 else ""
        ej = row.get(f"{epi_model}_epitope_jaccard", "") if epi_model else ""
        row["epitope_jaccard"] = ej
        t = thermo.get(qid, {})
        row["ddg_vs_trastuzumab"] = t.get("ddg_vs_trastuzumab", "")
        rows.append(row)

    designs = [r for r in rows if r["group"] == "denovo"]
    ref = next((r for r in rows if r["group"] == "baseline"), None)

    def sort_key(r):
        return -(r["binding_percentile_mean"] if r["binding_percentile_mean"] != "" else -1)

    designs.sort(key=sort_key)
    for i, r in enumerate(designs, 1):
        r["rank"] = i

    args.out.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) + ["rank"]
    with (args.out / "scores.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(designs + [r for r in rows if r["group"] != "denovo"])
    with (args.out / "calibration.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model", "n_with_measured_kd", "spearman_score_vs_logKD",
                    "own_batch_anchors", "own_batch_spearman", "note"])
        for model, c in sorted(calib.items()):
            w.writerow([model, c.get("n", 0), c.get("spearman_score_vs_logKD", ""),
                        c.get("own_batch_anchors", ""), c.get("own_batch_spearman", ""),
                        "negative Spearman is what we want: low KD should score high"])

    print(f"{len(designs)} designs ranked -> {args.out / 'scores.csv'}")
    if ref:
        print(f"  trastuzumab (baseline): {SCORE_KEY} "
              + ", ".join(f"{m}={ref.get(f'{m}_{SCORE_KEY}', '')}" for m in sorted(scores)))
    for model, c in sorted(calib.items()):
        print(f"  calibration {model}: n={c.get('n', 0)} "
              f"spearman(score, log KD)={c.get('spearman_score_vs_logKD')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
