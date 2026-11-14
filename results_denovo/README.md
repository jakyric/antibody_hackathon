# Snapshot of the de novo pipeline results

Written 2026-11-14 14:54 (local time) by `denovo/snapshot.py` from the pipeline's `results/denovo` folder. Scores are filters, not binding evidence: nothing here has been measured in the lab.

| file | rows | what |
|---|---|---|
| `submission_top20_pareto.csv` | 20 | FINAL LIST: top 20 designs by Pareto dominance (diversity-aware), with per-design reasons |
| `leaderboard_pareto.csv` | 639 | all gate-passing designs with Pareto rank, crowding, weakest objective and every metric (truncated to the first 605 of 639 rows) |
| `leaderboard.csv` | 1062 | all designs, weighted composite and gates (the comparison to the Pareto list) (truncated to the first 640 of 1062 rows) |
| `status.md` |  | funnel per stage, controls vs designs, top-10 with all metrics, Pareto summary, spend |
| `designs_rfd4jf.csv` | 1200 | every generated design (RFD4-Proteina on JapanFold): sequence, CDRs, settings |
| `rules.csv` | 1094 | hard rules and novelty gate per design (abguard class, identities, Cys checks) |
| `score_of3.csv` | 834 | OpenFold3 interface metrics per design (JapanFold recipe: MSA on, no template) |
| `score_orth.csv` | 834 | orthogonal predictors (OpenDDE, ESMFold-2, ...) per design |
| `props.csv` | 1449 | NanoMelt Tm, TNP flags, monoyield, sequence liabilities, nine-rule score per design |
