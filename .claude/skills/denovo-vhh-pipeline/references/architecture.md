# Architecture and data contract

Full text: `denovo/CONTRACT.md`. Summary.

## Stages, folders, outputs (all under `results/denovo/`, ignored by git; snapshots go to `results_denovo/`)

| Stage | Folder / entry point | Reads | Writes |
|---|---|---|---|
| S0 targets | `denovo/targets/make_targets.py` | `optimization/data/her2/her2_fv_1n8z.cif` (PDB 1N8Z) | `targets/*` (committed) |
| S1 design | `denovo/s1_design/run_design.py --auto` | targets, framework inputs, `score_*.csv` (feedback) | `designs_rfd4jf.csv`, `designs/*.cif`, `geom_rfd4jf.csv` |
| S2 rules | `denovo/s2_rules/rules.py --watch` | `designs_*.csv` | `rules.csv` |
| S3 scoring | `denovo/s3_score/run_s3.py watch` (`start_watch.sh`) | `designs_*.csv`, `rules.csv` | `score_of3.csv`, `score_orth.csv` |
| S3d | `denovo/s3d_protenix/` | Pareto list | `score_ptx.csv` |
| S4 props | `denovo/s4_props/props.py --watch` | `designs_*.csv`, `score_of3.csv`, `rules.csv` | `props.csv` |
| S5 | `denovo/rank.py` | everything above | `leaderboard.csv` |
| S6 | `denovo/pareto.py` | S5 output | `leaderboard_pareto.csv`, `submission_top20_pareto.csv` |
| status | `denovo/status.py`, `denovo/watch_rank.sh` | everything | `status.md` |

## Contract rules that make the loops safe to run in parallel

* One writer per file; every writer replaces files atomically (tmp + rename), so readers never see half a file.
* The design table is `designs_<generator>.csv` (one file per generator, no shared writers); `design_id` is unique and filesystem-safe; `denovo/dn_common.py:load_designs()` globs
  `designs_*.csv`. Cancelled generators are moved out of the glob (`parked_*.csv`).
* Every stage output is keyed by `design_id`, one row per design, with `is_control` and `stage_ts`; metric columns carry a stage prefix (`of3_`, `dde_`, `esm2_`, `ptx_`, `nm_`, `my_`, `tnp_`, `liab_`, `nb9_`).
* Every paid job has a stable name; its id is cached in `results/denovo/jobs/<name>.json`; re-running follows the cached job (`dn_common.JapanFold.submit`).
* Loops poll every 1.5-5 min and stop at their deadline or when `denovo/STOP` exists.
* Controls use ids starting with `ctl-` and generator `control`; they are scored for calibration and never ranked.

## Timeline and budgets used at the hackathon (all configurable)

Generation until 15:10, scoring until 15:40, props/rules until 15:45, final list 15:50, deadline 16:02. Caps: design $15, scoring $25 of a $100 JapanFold balance (actual design spend was a few dollars).

## Where each decision lives

`denovo/<stage>/decisions.md` (evidence first, choice second, what was rejected last), `denovo/s3_score/calibration.md` (scorer calibration), `denovo/s2_rules/gate_drivers.md` (novelty gate).
