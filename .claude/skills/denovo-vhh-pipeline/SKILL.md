---
name: denovo-vhh-pipeline
description: Set up, run, monitor, resume, stop and finish the de novo VHH design pipeline against HER2 domain IV (RFD4-Proteina on JapanFold, abguard novelty gate, OpenFold3/OpenDDE/ESMFold-2/Protenix scoring, NanoMelt/TNP/monoyield profiling, Pareto top-20). Use when the user asks to set up, install, run, start, resume, check on, tune, stop or summarise the de novo VHH / nanobody pipeline in this repository (the denovo/ folder), or to read its leaderboard, Pareto list, status or gate analysis.
---

# De novo VHH pipeline (denovo/)

Continuous, resumable pipeline: generate VHHs with RFD4-Proteina on an intact framework, gate them with the hackathon's hard rules and the abguard
track class, score the gate passers with OpenFold3 and independent predictors, profile Tm / yield / developability, rank by weighted composite and by
Pareto dominance. Scores are filters, never binding evidence: say so whenever you report results.

Read `denovo/README.md` for the architecture and `references/` (below) for evidence. Everything is run from the repository root.

## 1. Ask the user first (one short message, not a questionnaire)

* **Deadline**: when must the final list exist? (`PIPELINE_UNTIL=HH:MM`, default now + 2 h.) Generation stops 50 min before it.
* **Money**: JapanFold balance and caps (defaults: design $12, scoring $20; RFD4 costs about $0.003 per design, scoring dominates). Check `GET /v1/credits` via
  `bash denovo/run_pipeline.sh status`, which prints the balance.
* **Shared VM**: is the user allowing the shared team VM at all (monoyield on its GPU, VM abguard)? Default is NO. Never assume; the GPU is shared with other teams. Without it the gate runs on a local abguard replica (same track decision), and monoyield (yield objective) is replaced by the nine-rule score.
* **Nebius**: do they have `NEBIUS_CREDS` and the `nebius-scientific-ai` skill (needed only for Protenix-v2)?
* **Resume or fresh**: is `results/denovo/` already populated? Then resume; never restart from scratch (it would not re-pay: every job is cached by name).

## 2. Prerequisites

* `python3 >= 3.10`; `uv` is used if present. Optional: `micromamba` for NanoMelt/TNP (CPU).
* `JAPANFOLD_API_KEY` in the shell that launches Claude Code or the scripts (get one at https://japanfold.aiand.com/account/token). The user sets it; you never ask for it in
  chat, never print it, never write it to a file. Test with `[ -n "$JAPANFOLD_API_KEY" ]` only.
* Optional: `NEBIUS_CREDS` (path to the hackathon credentials file) and the `nebius-scientific-ai` skill; `PIPELINE_VM_HOST` / `PIPELINE_VM_KEY` for the VM.

## 3. Set up

```bash
bash denovo/setup.sh --check     # read-only report: what exists, what is missing (safe to run anytime)
bash denovo/setup.sh             # creates denovo/.venv from denovo/requirements.txt, validates targets and inputs, no-cost smoke test
bash denovo/setup.sh --tools     # optional: micromamba envs `nanomelt` and `tnp` (large download; Tm/TNP are skipped without them)
bash denovo/setup.sh --panel     # optional: known-VHH panel for the nearest-known identity column (public downloads, about 1 GB)
```

`setup.sh` exits non-zero only on hard failures (no python, broken targets, syntax errors). Warnings for optional parts are normal. Targets are committed
(`denovo/targets/`, derived from public PDB 1N8Z); they can be rebuilt with `denovo/targets/make_targets.py` if `optimization/data/her2/her2_fv_1n8z.cif` is present.

## 4. Run

```bash
bash denovo/run_pipeline.sh dry-run      # ALWAYS first: shows deadlines, caps, the exact command of each loop; spends nothing
bash denovo/run_pipeline.sh start        # rules gate, RFD4 design, scoring, properties, leaderboard + Pareto loops (background, resumable)
bash denovo/run_pipeline.sh status       # alive loops, funnel table, spend, last log lines
bash denovo/run_pipeline.sh ptx          # optional, once enough designs are scored: Protenix-v2 on Nebius for the top designs
bash denovo/run_pipeline.sh stop         # graceful (writes denovo/STOP); jobs already submitted keep running server-side and are re-adopted on restart
bash denovo/run_pipeline.sh final        # drain, rank.py + pareto.py + status.py, results_denovo/ snapshots, public-repo safety scan
```

Settings (environment variables on the same command): `PIPELINE_UNTIL`, `PIPELINE_GEN_MARGIN`, `PIPELINE_BUDGET_DESIGN`, `PIPELINE_BUDGET_S3`, `PIPELINE_TIER2`
(`dde,esm2,of3s1,af2ig`), `PIPELINE_ALLOW_VM=1` (only with the user's explicit yes). Example:
`PIPELINE_UNTIL=17:30 PIPELINE_BUDGET_S3=15 bash denovo/run_pipeline.sh start`.

Expected timeline after `start`: first designs in 10-20 minutes (JapanFold queue), rules rows within minutes, first OpenFold3 scores 15-30 minutes after designs,
Pareto list meaningful once about 50 designs are fully measured. The JapanFold fleet can be saturated ("account_backlog" in the S3 log): that is queueing, not a fault.

## 5. Safety rules (non-negotiable)

1. **Secrets**: never print, log, commit or copy `JAPANFOLD_API_KEY`, the contents of `$NEBIUS_CREDS`, ssh keys or VM addresses. The repository is public. Run
   `python3 denovo/scan_public.py .` before any commit or push and fix every finding.
2. **Spend**: respect the caps; check the balance before and after large batches; never resubmit a job (job ids are cached in `results/denovo/jobs/`; a repeated
   command follows the existing job). After a timeout use status/resume, not a new submit.
3. **No VM GPU without the user's explicit permission.** Monoyield needs the VM's GPU and stays paused unless `PIPELINE_ALLOW_VM=1` was granted; even then it runs one
   bf16 batch at a time and skips a cycle if less than 28 GB GPU memory is free. Never kill or touch other people's processes or containers.
4. **No RFD4-Proteina (or other design models) on a shared VM** without permission; design runs on JapanFold.
5. **No git push** unless the user asks; commit only specific paths.
6. Report honestly: scores are filters; OpenFold3 on JapanFold is untemplated with MSA and near decoy level on these designs, the VM's templated OpenFold3 (no MSA) misfolds the VHH, and **no scorer reproduces the designed pose** (`results_denovo/pose_check_vm.csv`; see `references/failure-modes.md`). Never present iPTM/ipSAE as binding evidence.
7. **Never run `props.py --ids ...` on the default output** (it rewrites `props.csv` and drops all other rows). Add profiled designs with `run_vm.sh` + `denovo/s4_props/merge_side.py`.

## 6. Reading the results

* `results/denovo/status.md` (refreshed every 5 min): funnel per stage (designs, rules, OpenFold3, orthogonal, props, gate, composite), controls vs designs, top-10 with all metrics, Pareto summary, JapanFold balance.
* `results/denovo/leaderboard.csv`: weighted composite (priors in `denovo/rank.py`: bind 0.40, epitope 0.15, Tm 0.20, yield 0.15, developability 0.10, nine rules 0.05; missing evidence counts as 0.3, so partly measured designs cannot top the list).
* `results/denovo/leaderboard_pareto.csv` and **`submission_top20_pareto.csv`**: the final list. Columns: `pareto_rank` (1 = front), `maximin` (the weakest objective's percentile), `weakest`, `pareto_complete` (all
  used objectives measured), `pareto_reason`, `rank_reason`. Objectives with under 30% coverage are dropped automatically and printed by `pareto.py`.
* Gates: `gate_vm` = abguard class de_novo (identity to trastuzumab VH < 70%, hard rules, canonical Cys); `gate_strict` adds a 2-point margin (<= 0.68). Nearest-known-VHH identity is reported, not gated.
* Calibration: `denovo/s3_score/calibration.md` and `gates.json`; evidence per stage in each `decisions.md`.

## 7. Tuning (what to change, where)

| Want | Change |
|---|---|
| More designs per hour | `--max-inflight` / `--jobs-per-wave` in `run_pipeline.sh` (`CMDS[design]`); the queue, not money, is the limit |
| Longer CDR3 (more designs pass the gate) | `s1_design/run_design.py` grid (`c3n` block = 15-20 designed residues, C104-A105-A106 fixed); read `s2_rules/gate_drivers.md` first |
| Other frameworks / hotspot sets | `s1_design/prep_inputs.py` (frameworks), `targets/epitope.json` (hotspot sets); rerun prep, then `run_design.py --auto --dry-run` |
| Ranking priorities | weights in `denovo/rank.py` (`--w-bind` ...), objectives in `denovo/pareto.py --objectives bind,comp,therm,yield,dev` |
| Stricter gates | `PIPELINE_TIER2`, `s3_score/gates.json` (decoy 95th percentiles), `rank.py` margin (0.68) |
| Add a predictor | write `results/denovo/score_<name>.csv` keyed by `design_id` with a 0..1 `<name>_score`, merge it in `rank.py` like `score_ptx.csv` |

## 8. Overview slides

`denovo/slides/build_slides.py --plan` regenerates the 14 slide files from the result files and prints the publish plan (see `denovo/slides/README.md`); publish with the
Artifact tool to the same URL. Numbers are read from the result files, anything without evidence shows "pending".

## 9. If something goes wrong

Read `references/failure-modes.md` first (JapanFold backlog, OpenFold3 decoy-level scores, NanoMelt out of memory on small GPUs, abguard identity depends on length,
zero-Cys designs, pandas read-only sort, sshd resets on the shared VM, and how each was handled). `bash denovo/run_pipeline.sh status` and the logs in `denovo/.run/` and
`denovo/s3_score/watch.log` come before any guesswork.

## 10. Optional re-scoring: VM OpenFold3 / ESMFold-2, pose check, cropped HER2

Only with the user's explicit yes for the shared VM (`PIPELINE_VM_HOST`, `PIPELINE_VM_KEY`; GPU shared, workers wait while < 20 GB is free).

* Pose check: `python denovo/validation/pose_check_vm.py` compares VM OpenFold3, VM ESMFold-2 and OpenDDE with the designed pose (`results_denovo/pose_check_vm.csv`). Run it before trusting any binding score.
* Cropped HER2 (539-607) re-run, to limit random placement: `denovo/vm_templated/make_crop.py`, workers `run_vm_crop.sh` on the VM, `collect_crop.py`, and `denovo/s3_score/run_crop.py` for OpenDDE on JapanFold.
  Commands and rationale are in `denovo/README.md`; results land in `results_denovo/crop/` (read them there, do not assume they improve the ranking).
* New profiled designs: `merge_side.py` (see rule 7).

## References

* `references/architecture.md`: stages, files, data contract between stages, timeline and budgets
* `references/gate-analysis.md`: what the abguard track gate measures and why length drives it
* `references/decision-log.md`: every choice, the evidence, what was rejected
* `references/failure-modes.md`: problems met during the hackathon and their fixes (including the pose-check and `props.py --ids` traps)
