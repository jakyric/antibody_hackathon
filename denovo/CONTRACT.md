# De novo VHH pipeline: contract between stages

Goal: de novo-track VHHs against HER2 domain IV (trastuzumab epitope), made by RFD4-Proteina on a solid
framework, scored by OpenFold3 plus orthogonal predictors, then profiled for yield, Tm and developability.
Hackathon day is 14 Nov 2026 and everything is time-boxed: every stage must be **resumable, idempotent and
useful at any moment** (partial outputs are fine, a late perfect output is not).

## Stages and owners

| Stage | Dir | Input | Output (under `results/denovo/`) |
|---|---|---|---|
| T  targets (orchestrator) | `denovo/targets/` | 1N8Z file in `data/her2/` | `targets/READY` marker + files below |
| S1a design, JapanFold `rfd4` | `denovo/s1_design/` | targets | `designs_rfd4jf.csv`, `designs/*.cif` |
| ~~S1b design, VM Foundry `antibody-2026-10-30`~~ | CANCELLED by the user: **no RFD4-Proteina on the VM** | | RFD4 runs only on JapanFold (S1a) |
| S1c design, JapanFold `boltzgen` `nanobody-anything` | `denovo/s1c_boltzgen/` | targets | `designs_boltzgen.csv`, `designs_boltzgen_metrics.csv` (generator-native scores keyed by design_id), `designs/*.cif` |
| S2  rules + framework diversification | `denovo/s2_rules/` | any `designs_*.csv` | `rules.csv`, `designs_fwdiv.csv` |
| S3  OF3 primary scoring + orthogonal binding (OpenDDE-abag, ESMFold2, ...) | `denovo/s3_score/` | any `designs_*.csv` | `score_of3.csv`, `score_orth.csv` |
| S4  yield / Tm / developability | `denovo/s4_props/` | any `designs_*.csv` | `props.csv` |
| S5  ranking (orchestrator) | `denovo/rank.py` | all of the above | `leaderboard.csv` |
| S6  Pareto scoring, the last step (orchestrator) | `denovo/pareto.py` | S5 output | `leaderboard_pareto.csv`, `submission_top20_pareto.csv` (the final list) |

Every worker only writes inside its own stage dir and its own output files. Never edit `common.py`,
`CONTRACT.md`, or another stage's files: report the need instead.

## Shared files (made by T)

`denovo/targets/`: `her2_domIV.cif` (HER2 residues 489-607 as ONE chain `A`, original numbering kept in
`auth_seq_id`, mmCIF), `her2_domIV.fasta`, `her2_domIV.a3m` (MSA, copy of the existing one), `ref_1n8z.cif`
(crystal Fv+HER2, chains A=VL B=VH C,D=HER2), `epitope.json` (residues of HER2 within 4.5 A of trastuzumab,
with per-residue contact counts and 3 ranked hotspot sets), `README.md` (numbering notes). Workers must wait for
`denovo/targets/READY` before relying on these (poll, do not recreate).

## Design table schema (`designs_<generator>.csv`, one file per generator, no shared writers)

`design_id, parent_id, generator, framework, seed, hotspot_set, vhh_seq, cdr1, cdr2, cdr3, cdr_scheme,
design_cif, job_id, params_json, created_utc`

* `design_id`: `<fw>-<gen>-s<seed>-<NNN>`, unique across everything, filesystem-safe, <= 48 chars.
* `vhh_seq`: the VHH only, no tags, no linker, one-letter, upper case. `parent_id` empty unless derived.
* `design_cif`: path relative to repo root of the design (VHH + target if the generator gives it), or empty.
* Controls (known binders, framework-only, decoys) use `generator = control` and are NEVER mixed into the
  de novo ranking by design_id prefix `ctl-`.
* Read all with `common.load_designs()` (globs `designs_*.csv`).

Per-stage outputs are keyed by `design_id`, one row per design (or per design x model), include a
`stage_ts` UTC column, and carry `is_control` (0/1). Column names are prefixed by stage metric family
(`of3_iptm`, `of3_ipae`, `of3_ipsae`, `dde_iptm`, `esm2_iptm`, `nm_tm`, `tnp_flags`, `my_yield`, ...).

## Rules for all workers

1. **Secrets**: never print, log, commit or copy `JAPANFOLD_API_KEY` or the file at `$NEBIUS_CREDS`. Test
   presence with `[ -n "$VAR" ]`. Never open the hackathon credentials file.
2. **No duplicate paid jobs**: every submit goes through a job-id cache file in `results/denovo/jobs/`
   (`common.cached_job`). Re-running a script follows the existing job. After a timeout, recover the job id,
   do not resubmit. Send long jobs early, then do other work while they run.
3. **Budget** (JapanFold balance about $93, shared by the team; the user lifted the earlier tight caps): total for
   this pipeline <= $65, keep >= $20 in reserve. Per worker: S1a $15, S1c $15, S3 $25, S2 $3, S4 $0 (use local + VM).
   `GET /v1/credits` before and after big batches; stop and report if a batch would exceed the worker's cap.
   rfd4 costs about $0.014 per job (up to 8 designs per job).
4. **VM etiquette** (`ssh -i $PIPELINE_VM_KEY $PIPELINE_VM_HOST`, H100 shared with teammates, already busy):
   work only under `~/moritz/` (S1b: `~/moritz/denovo_s1b`, S4: `~/moritz/monoyield-predictor`, S2: `~/moritz/denovo_s2`);
   one GPU container/process at a time per worker; check `nvidia-smi` first; use bf16 where available; do not kill
   or touch other users' processes or directories; no `docker system prune`; no sudo.
5. **No git**: do not commit or push. The project root is not a git repo; the orchestrator packages at the end.
6. **Reproducibility is 20% of the score and reasoning 40%**: every stage ships a `README.md` in its dir
   (what, why, exact commands, parameters, versions/commits of tools, known failure modes) and writes a
   `decisions.md` log (what you chose, the evidence, what you rejected). Be honest about weak evidence.
7. **Smoke first, then scale**: prove the path on 1-2 inputs, say so in README, then run the real batch.
   Do not claim a score means binding; scores are filters, and calibration against controls/decoys is the point.
8. Use `denovo/.venv/bin/python` (pandas/numpy/scipy/sklearn) unless you need a micromamba env
   (`tnp`, `nanomelt`, `flashabb` exist; see `tools/README.md`). Install nothing system-wide.
9. Final report (your last message, <= 300 words): what exists and works, exact run commands, output paths and
   row counts, key numbers, what failed or is weak, what you recommend the orchestrator does next. Nothing else.

## Timeline (local CET, 14 Nov 2026): final ranked list due 16:02

| Time | What |
|---|---|
| now - 14:30 | smoke tests pass; first real rows from S1a/S1c; S2/S3/S4 `--watch` loops started |
| continuous | each design flows generate -> OF3 (tier 1) -> orthogonal (tier 2, top ~35%) -> S2 diversifies top parents -> S3 re-scores variants -> S4 props (tiered) -> `denovo/rank.py` every 5 min (`denovo/watch_rank.sh`, funnel in `results/denovo/status.md`) |
| 15:10 | generators stop launching new jobs |
| 15:30 | S2 stops creating framework variants |
| 15:40 | S3 stops launching; 15:45 S4/S2 rules stop; drain in-flight work |
| 15:50 | final `rank.py`, `submission_top20.csv`, rationale; 16:02 hard deadline |
`touch denovo/STOP` stops every watch loop.

## DECISION 14:25: BoltzGen stopped, RFD4 only

The user stopped BoltzGen (S1c) to focus on RFD4. Its loop process was killed, its in-flight JapanFold jobs cancelled
(8 of 20 jobs were still queued/running; 12 had finished), and its files parked as `results/denovo/parked_boltzgen*.csv`
(raw designs have zero Cys, 160 designs, abguard-passable only after Cys restoration). Not ranked, not scored.
S1a (JapanFold `rfd4`) is the only generator.

## DECISION 14:15, supersedes the framework-diversification text below

The user and the organisers: keep the framework INTACT; fully new CDRs should land under 70%. No framework
diversification runs (S2's `diversify.py` is parked; its 8 variants are in `results/denovo/parked_fwdiv.csv`).
Measured on the first 96 RFD4 designs (intact hNbBCII10-FGLA, all CDRs redesigned): identity to trastuzumab VH 0.65-0.74
(mean 0.70), so the VM `abguard` class is de_novo for only 53%; the pass rate is driven by VHH length / CDR3 length
(CDR3 <= 14: 0-21%, 15-17: 50-62%, >= 18: 80-100%), not by the framework. The gate is therefore `abguard de_novo`
(+ 2-point margin), applied as a FILTER; generators are steered to longer CDR3 instead of mutating the framework.
Identity to the nearest known VHH in our panel (~0.80) is reported, not gated. If a good design is ever blocked
only by identity, the fallback repair is CDR-adjacent framework positions first, never IMGT 42 (F40 in the paper's
numbering) nor the C-A-A at the FR3/CDR3 junction. The nine-rule CDR2/CDR3 score from the monoyield paper (Fig. 5h) is
a column (`nb9_score`) and a small leaderboard prior.

## Tracks, and the rule that drives stage S2 (historical, see the decision above)

De novo track = under 70% identity to ANY known antibody (VHH vs known VHHs and heavy chains) AND at least
2 CDRH3 mutations from trastuzumab; <= 250 aa; no unpaired Cys; immunum numbering confidence >= 0.1; must pass
`abguard`. A fixed natural framework with brand-new CDRs is about 75% identical to its parent, so frameworks
must be diversified by roughly 10+ substitutions (hallmarks, core and Cys22/Cys92 held). The team VM's
`abguard` compares only against trastuzumab; the jury's rule is unconfirmed. Treat "under 70% to the nearest known
VHH" as the target and report both numbers.
