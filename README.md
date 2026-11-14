# AI x Antibodies hackathon (Oslo, 14 Nov 2026): improving trastuzumab against HER2

Two independent tracks, one folder each. Neither track stores any key or credential; JapanFold and Nebius access come from
environment variables (`JAPANFOLD_API_KEY`, `NEBIUS_CREDS`).

| Track | What | Start here |
|---|---|---|
| **De novo VHH** | A continuously running, filter-first pipeline that designs new VHHs against the trastuzumab epitope on HER2 domain IV with RFD4-Proteina (JapanFold), gates them with the hackathon's hard rules and `abguard`, scores them with OpenFold3 plus independent predictors (OpenDDE, ESMFold-2, Protenix-v2), profiles yield, Tm and developability (monoyield, NanoMelt, TNP), and picks the final list by Pareto dominance. | [`denovo/README.md`](denovo/README.md) and the Claude Code skill [`.claude/skills/denovo-vhh-pipeline`](.claude/skills/denovo-vhh-pipeline/SKILL.md) (ask Claude to "set up the de novo VHH pipeline"); quickstart `bash denovo/setup.sh && bash denovo/run_pipeline.sh dry-run`. Latest results snapshot: [`results_denovo/`](results_denovo/) (written by `denovo/run_pipeline.sh final`). |
| **Optimisation** | Trastuzumab CDRH3 variants: binding classifiers on the 524k-variant FACS screen, BLOSUM-45 library design (three variants A, B and C, the last following the published recipe), and refolding of the shortlists with OpenFold3, OpenDDE and Boltz-2 against the 1N8Z crystal. OpenFold3 was the most reliable refolder; variants A and B gave 7 and 3 vetted designs, variant C only 1 of 13. | [`optimization/README.md`](optimization/README.md) (explainer with workflow diagrams), code in `optimization/baseline/`, run records in `optimization/nebius-runs/` and `optimization/results/`. |

Repository layout

```
denovo/             de novo VHH pipeline: stages S0-S6, setup.sh, run_pipeline.sh, CONTRACT.md
results_denovo/     small result snapshots of the de novo pipeline (leaderboard, Pareto list, status)
optimization/       the earlier optimisation-track work (baseline/, nebius-runs/, results/, requirements.txt, README.md)
.claude/skills/     Claude Code skill that sets up, runs and explains the de novo pipeline
```
