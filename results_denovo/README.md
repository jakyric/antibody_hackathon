# Snapshot of the de novo pipeline results

Written 2026-11-14 15:00 (local time) by `denovo/snapshot.py` from the pipeline's `results/denovo` folder. Scores are filters, not binding evidence: nothing here has been measured in the lab.

| file | rows | what |
|---|---|---|
| `submission_top20_pareto.csv` | 20 | FINAL LIST: top 20 designs by Pareto dominance (diversity-aware), with per-design reasons |
| `leaderboard_pareto.csv` | 283 | all gate-passing designs with Pareto rank, crowding, weakest objective and every metric |
| `leaderboard.csv` | 463 | all designs, weighted composite and gates (the comparison to the Pareto list) |
| `status.md` |  | funnel per stage, controls vs designs, top-10 with all metrics, Pareto summary, spend |
| `designs_rfd4jf.csv` | 449 | every generated design (RFD4-Proteina on JapanFold): sequence, CDRs, settings |
| `rules.csv` | 463 | hard rules and novelty gate per design (abguard class, identities, Cys checks) |
| `score_of3.csv` | 381 | OpenFold3 interface metrics per design (JapanFold recipe: MSA on, no template) |
| `score_orth.csv` | 381 | orthogonal predictors (OpenDDE, ESMFold-2, ...) per design |
| `props.csv` | 818 | NanoMelt Tm, TNP flags, monoyield, sequence liabilities, nine-rule score per design |
