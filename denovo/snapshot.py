#!/usr/bin/env python3
"""Copy small, public-safe result snapshots of a pipeline run into results_denovo/ (no structures, no job ids, no logs).

    python3 denovo/snapshot.py [--src results/denovo] [--dst results_denovo] [--max-mb 1.5]

Files above --max-mb are cut to their first rows (the leaderboards are sorted best-first, so the head is what matters).
Writes results_denovo/README.md with the timestamp, row counts and what each file is.
"""
import argparse
import csv
import datetime as dt
import pathlib
import shutil

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="results/denovo")
ap.add_argument("--dst", default="results_denovo")
ap.add_argument("--max-mb", type=float, default=1.5)
a = ap.parse_args()
src, dst = pathlib.Path(a.src), pathlib.Path(a.dst)
WHAT = {
    "submission_top20_pareto.csv": "FINAL LIST: top 20 designs by Pareto dominance (diversity-aware), with per-design reasons",
    "leaderboard_pareto.csv": "all gate-passing designs with Pareto rank, crowding, weakest objective and every metric",
    "leaderboard.csv": "all designs, weighted composite and gates (the comparison to the Pareto list)",
    "status.md": "funnel per stage, controls vs designs, top-10 with all metrics, Pareto summary, spend",
    "designs_rfd4jf.csv": "every generated design (RFD4-Proteina on JapanFold): sequence, CDRs, settings",
    "rules.csv": "hard rules and novelty gate per design (abguard class, identities, Cys checks)",
    "score_of3.csv": "OpenFold3 interface metrics per design (JapanFold recipe: MSA on, no template)",
    "score_orth.csv": "orthogonal predictors (OpenDDE, ESMFold-2, ...) per design",
    "score_ptx.csv": "Protenix-v2 (Nebius) metrics for the top designs",
    "props.csv": "NanoMelt Tm, TNP flags, monoyield, sequence liabilities, nine-rule score per design",
}
dst.mkdir(parents=True, exist_ok=True)
lines = []
for name, what in WHAT.items():
    f = src / name
    if not f.exists():
        continue
    out = dst / name
    n_rows = None
    if f.suffix == ".csv":
        rows = list(csv.reader(open(f, newline="")))
        n_rows = len(rows) - 1
        # structures and job ids are not part of the snapshot: drop the columns that point at them
        drop = [i for i, h in enumerate(rows[0]) if h == "job_id" or h.endswith("_cif") or h == "design_cif"]
        if drop:
            rows = [[c for i, c in enumerate(r) if i not in drop] for r in rows]
        if f.stat().st_size > a.max_mb * 1e6:
            keep = max(100, int(len(rows) * a.max_mb * 1e6 / f.stat().st_size))
            rows = rows[:keep]
            what += f" (truncated to the first {keep - 1} of {n_rows} rows)"
        with open(out, "w", newline="") as fh:
            csv.writer(fh).writerows(rows)
    else:
        shutil.copy(f, out)
    lines.append(f"| `{name}` | {n_rows if n_rows is not None else ''} | {what} |")
(dst / "README.md").write_text(
    "# Snapshot of the de novo pipeline results\n\n"
    f"Written {dt.datetime.now().strftime('%Y-%m-%d %H:%M')} (local time) by `denovo/snapshot.py` from the pipeline's `results/denovo` folder. "
    "Scores are filters, not binding evidence: nothing here has been measured in the lab.\n\n"
    "| file | rows | what |\n|---|---|---|\n" + "\n".join(lines) + "\n")
print(f"snapshot: {len(lines)} files -> {dst}")
