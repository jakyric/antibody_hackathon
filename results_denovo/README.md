# Snapshot of the de novo pipeline results

Written 2026-11-14 16:53 (local time) by `denovo/snapshot.py` from the pipeline's `results/denovo` folder. Scores are filters, not binding evidence: nothing here has been measured in the lab.

| file | rows | what |
|---|---|---|
| `submission_top20_pareto.csv` | 20 | FINAL LIST: top 20 designs by Pareto dominance (diversity-aware), with per-design reasons |
| `leaderboard_pareto.csv` | 43 | all gate-passing designs with Pareto rank, crowding, weakest objective and every metric |
| `leaderboard.csv` | 1592 | all designs, weighted composite and gates (the comparison to the Pareto list) (truncated to the first 391 of 1592 rows) |
| `status.md` |  | funnel per stage, controls vs designs, top-10 with all metrics, Pareto summary, spend |
| `designs_rfd4jf.csv` | 1544 | every generated design (RFD4-Proteina on JapanFold): sequence, CDRs, settings |
| `rules.csv` | 1500 | hard rules and novelty gate per design (abguard class, identities, Cys checks) |
| `score_of3.csv` | 1235 | OpenFold3 interface metrics per design (JapanFold recipe: MSA on, no template) |
| `score_orth.csv` | 1235 | orthogonal predictors (OpenDDE, ESMFold-2, ...) per design |
| `score_ptx.csv` | 16 | Protenix-v2 (Nebius) metrics for the top designs |
| `props.csv` | 356 | NanoMelt Tm, TNP flags, monoyield, sequence liabilities, nine-rule score per design |

## Validation and re-run tables (snapshot_extra.py)

| file | rows | what |
|---|---|---|
| `pose_check_vm.csv` | 1413 | pose RMSD, scRMSD, epitope-contact fraction, pDockQ2 per design for VM OpenFold3 (templated, no MSA), VM ESMFold-2 and OpenDDE vs the designed pose |
| `dde_pose_check.csv` | 768 | independent re-check of the OpenDDE pose RMSD / scRMSD (HER2 matched by residue number) |
| `score_dde_extra.csv` | 847 | pDockQ2 of the OpenDDE predictions |
| `score_vmof3.csv` | 559 | VM OpenFold3 (HER2 template, no MSA) and ESMFold-2 interface scores per design (iPTM, pair iPTM, ipSAE) |
| `props_nm.csv` | 332 | NanoMelt Tm per design (side file read by rank.py) |
| `props_tnpvm.csv` | 302 | TNP flags per design (side file read by rank.py) |
| `props_my.csv` | 240 | monoyield per design (side file read by rank.py) |
| `validity.csv` | 1480 | structure / sequence validity checks of the RFD4 outputs |
| `gapcheck.csv` | 2080 | distance of each designed VHH to the unresolved HER2 loop 581-590 (AlphaFold loop) |
| `crop/dde_crop.csv` | 170 | cropped-HER2 (539-607) re-run: scores / metrics |
| `crop/of3_crop.csv` | 170 | cropped-HER2 (539-607) re-run: scores / metrics |
| `crop/set.csv` | 170 | cropped-HER2 (539-607) re-run: design set |
| `dde_compare.csv` | 10 | designed vs OpenDDE pose of the top 10 (the PyMOL session itself is not stored) |

## top20/

The 20 designs of `submission_top20_pareto.csv`: `top20.fasta` (VHH sequences) and the RFD4 design complexes (`NN_<design_id>.cif.gz`: chain A = HER2 489-580, chain B = designed VHH,
chain C = HER2 591-607). All 1,544 design sequences are in `designs_rfd4jf.csv`. These are designed structures, not predictions; no scorer reproduced the designed binding poses (see `pose_check_vm.csv`).
