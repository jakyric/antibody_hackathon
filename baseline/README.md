# Trastuzumab CDRH3 binding baseline (local)

A small, fast "will this CDRH3 be a high HER2 binder?" model, trained on the 524,346-variant
FACS screen, plus a gate script for filtering designs. Runs on CPU in about 40 seconds.

## Run

```bash
baseline/.venv/bin/python baseline/train_eval.py              # train, evaluate, save outputs/
baseline/.venv/bin/python baseline/score_cdrh3.py --seqs WGSQGFYAYD ARNDCQEGHI
baseline/.venv/bin/python baseline/score_cdrh3.py --csv designs.csv --col cdrh3 --out scored.csv
```

Environment: `baseline/.venv` (pandas, scikit-learn, scipy, h5py). Data lives in `data/` (git-ignored):
the FACS and SPR tables from the team bucket and the Mehta et al. split files from Zenodo 10831512.

## Files

| File | What it does |
|---|---|
| `common.py` | Encoding, data loading, novelty index, and the authors' CNN re-implemented in numpy (loads their `.h5`, no TensorFlow, no pickles) |
| `train_eval.py` | Trains three baselines, compares with the authors' CNN and a trivial score, tests on SPR and the Absci designs |
| `score_cdrh3.py` | The gate: P(high binder) from two models plus the novelty-rule check |
| `blosum_diversity.py` | Picks mutable CDRH3 positions from structure, FACS and SPR evidence, ranks substitutions with the paper's BLOSUM-45 recipe, filters for novelty, liabilities and the gate, runs an immunum-based "abguard-lite" check and writes a diverse shortlist |
| `outputs/blosum_candidates.csv` | The shortlist (also `blosum_all_candidates.csv`, `position_evidence.csv`) |
| `outputs/metrics.json`, `outputs/hgb_full.joblib` | Results and the final model (trained on all FACS data) |

Refolding the shortlist (Fv + HER2 domain IV, scored against the 1N8Z crystal): `refold_prepare.py` builds the 16
complexes (12 candidates + 4 controls), `refold_nebius.sh` runs Boltz-2 on Nebius, `refold_japanfold.py` runs
JapanFold batches, `structure_tools.py` compares with the crystal, `refold_analyze.py` and `refold_consensus.py`
summarise. Outputs: `outputs/refold_metrics.csv`, `outputs/refold_consensus.csv`. Finding so far: only OpenFold3
separated the measured binders from the measured non-binders; Boltz-2 reproduces the crystal pose for everything and
OpenDDE fails trastuzumab itself, so "pose recovered" is weak evidence, not proof of binding.

`blosum_diversity.py` also needs `data/paper/spr_submitted.csv` from the paper's GitHub repo
(oxpig/Tz_her2_affinity_and_beyond), which carries the design-method label of every SPR sequence, plus
`blosum` and `immunum` in the venv. "abguard-lite" is an approximation: the real `abguard` was not available.

## Results (authors' clonotype split, test set of 78,651)

| Model | AUC | AUC, test sequences with no train neighbour within 1 mutation |
|---|---|---|
| LR on one-hot positions (additive) | 0.944 | 0.932 |
| LR with pairwise interactions | 0.985 | 0.979 |
| Gradient-boosted trees (ours) | 0.995 | 0.993 |
| Authors' CNN (their weights) | 0.997 | 0.995 |
| Trivial: fewer mutations vs trastuzumab | 0.753 | 0.788 |

Label: 1 = "high" binder class (32.8% of the library), 0 = medium or low.

## External checks (sequences not in the FACS table)

| Model | SPR: AUC, measured binders vs no reportable KD | SPR: Spearman with -KD (152 binders) | Absci designs: share predicted high (198, all binders) |
|---|---|---|---|
| LR one-hot | 0.65 | 0.13 | 59% |
| LR pairwise | 0.72 | 0.40 | 60% |
| Gradient-boosted trees | 0.75 | 0.35 | 51% |
| Authors' CNN | 0.74 | 0.38 | 49% |
| Trivial: fewer mutations vs trastuzumab | 0.81 | 0.43 | n/a |

## Caveats

- In-library accuracy is near ceiling, but on designs made by other methods the learned models are only
  moderately predictive and barely beat "stay close to trastuzumab".
- The SPR designs were themselves pre-screened with a CNN from the same lineage (P(binder) above 90%, per
  Chinery et al. Methods), apart from controls and the length-9 AbLang designs. The SPR comparison is therefore
  range-restricted and unfair to learned scores, which already saturate there; it favours the "fewer mutations"
  score. Read the SPR rows as a weak test, not as proof the models fail on new designs.
- The score is a class probability (high vs not), not an affinity. Scores saturate near 1.0 for anything
  close to trastuzumab, so use ranks or logits, not raw probabilities.
- Only 4 SPR variants beat trastuzumab after removing FACS overlaps, so top-end ranking cannot be judged.
- The library is sparse at 3 to 4 mutations (7,900 variants), which is the range the novelty rule forces.
- `challenge_ok` follows the Apheris docs' definition (identity to trastuzumab CDRH3 at most 70%, and at least
  2 positions away from every FACS and SPR sequence). `abguard` is the official checker.
- `no_reportable_kd` is treated as "non-binder" in the SPR check. The docs describe it only as "measured with
  no KD reported", so that label is uncertain.

## Sources

- Chinery, Mehta, Greiff, Deane et al., "Baselining the Buzz. Trastuzumab-HER2 Affinity, and Beyond!",
  bioRxiv 2024, doi:10.1101/2024.03.26.586756
- Code: github.com/oxpig/Tz_her2_affinity_and_beyond
- Data and weights: Zenodo record 10831512 (CC BY 4.0)
