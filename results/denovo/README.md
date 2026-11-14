# De novo track: VHH binders against HER2

All numbers here are computational predictions. The lab measures binding, epitope and
thermostability after the event, so every statement below says *predicted*, never *binds*.

**The submission is `denovo_submission.csv`: six VHH designs, all six accepted by the official checker
(`team-validate` exit code 0, novelty class `de_novo`).**

## The finding that decided the track

The starting point was `configs/her2_all6cdr.yaml`: design all six CDRs on trastuzumab's own Fv
framework. It runs, and the designs pass the entry rules — but they are not de novo.

The official checker measures `identity_to_trastuzumab` **over the whole sequence**, not over the CDRs. For 47
designs from that spec the identity was 0.877 to 0.890, so the official checker put every one of them in the
`novel` class, meaning *a trastuzumab variant*, which is the optimisation track. No amount of CDR
redesign gets under 70 % while the framework is trastuzumab's.

So the de novo track needs a different framework. `configs/her2_vhh_denovo.yaml` designs the three
CDRs of a **VHH (nanobody)** on the caplacizumab framework from the organisers' `Data/03-vhh` file,
against HER2 domain IV, with the same four hotspot atoms as the organisers' HER2 spec — that is,
aimed at trastuzumab's epitope. The official checker's verdict on those: `de_novo`, identity to trastuzumab
0.63 to 0.70.

Getting the framework file usable took one step worth recording: the organisers' VHH file is built
to be run as-is against its own target, with placeholder residues where the CDRs go and the
conditions baked in as columns. Reused against HER2 it fails with
`Unresolved link valence at H/23/13-13/UNK`. `scripts/vhh_framework_cif.py` writes a clean
framework file from it.

## What was run

Everything ran on a single H100 through `record_run.py`, so each design traces back to a weight, a
seed, a command and a commit (`results/ledger.csv`).

| Run | Weight | Seed | Designs |
| --- | --- | --- | --- |
| `denovo-vhh-antibody-2026-10-30-s11` | `antibody-2026-10-30` | 11 | 40 |
| `denovo-vhh-antibody-2026-10-30-s12` | `antibody-2026-10-30` | 12 | 40 |
| `denovo-vhh-antibody-2026-10-30-s13` | `antibody-2026-10-30` | 13 | 40 |
| `denovo-antibody-2026-10-30-s1` (Fv, secondary) | `antibody-2026-10-30` | 1 | 33 |
| `denovo-her2-novel-2026-11-09-s2` (Fv, secondary) | `her2-novel-2026-11-09` | 2 | 25 |

The two HER2-tuned weights only make sense for the trastuzumab framework, so the VHH spec uses the
general antibody weight over three seeds instead. The two Fv runs were stopped early, once the
official checker's result showed they belong in the other track.

## Gate before ranking

Known binders are filtered out before anything is folded. Two layers, as `candidate-gate`
prescribes:

1. **Local** (`src/abx/gate.py`, `scripts/denovo_collect.py`): unpaired cysteines, length, identity
   to trastuzumab, and the distance from the CDR3 to all 525,063 provided CDRH3s. We require 3
   mutations where the de novo rule asks for 2, for margin against the larger reference the lab runs.
   118 of 120 VHH designs passed; the two that failed were at or above 70 % identity to trastuzumab.
2. **The official checker**, the authority. All 118 passed, class `de_novo`, nearest known CDR3 at
   least 9 mutations away.

On the Fv designs the two layers disagreed in a way worth keeping: the local gate passed every
design, while the official checker rejected two because their CDRH3 sits 2 mutations from a known antibody in
PLAbDab and the PDB (`pdb:10ie_H`, `plabdab_AQN27791_AQN27792`) — sequences our local index does not
contain. The trastuzumab control row failed at 0 mutations, as it must.

## How the folding scores were calibrated

Before ranking on an interface score we measured what it is worth. Twenty-five variants with
**measured** KD were folded with the same models, template and seed, and we correlated interface
ipTM with log KD:

| Folding model | Template | Spearman (score vs log KD), n=25 |
| --- | --- | --- |
| Boltz-2 | yes | **-0.41** |
| OpenFold3 OpenBind | yes | -0.24 |
| Boltz-1 | no | +0.23 |

Negative is the direction we want. Boltz-2 is weakly informative, OpenFold3 weaker, and the model
without a template points the wrong way and scores trastuzumab itself at 0.12 interface ipTM, so it
cannot model this complex at all. We rank on the two template models and keep the templateless one
as a control, not as evidence. **A Spearman of -0.41 orders designs loosely, not reliably**: the
ranking below separates the clearly good from the clearly bad, and should not be read as a precise
order.

For the VHH designs the organisers' `refold-score` applies, and unlike the HER2 Fv example its
`pass` column is populated, with `min_ipae` and the distance the binder moved from its designed pose.

## The three properties

- **Binding**: `min_ipae` from both folding models (lower is better), as a percentile within each
  model, averaged. 2 of 120 designs pass the organisers' oracle in at least one model.
- **Epitope**: HER2 residues within 4.5 Å of the binder in each predicted complex, against the
  trastuzumab interface in the experimental structure 1N8Z (19 residues), as Jaccard overlap. Median
  across the 120 designs is 0.16, best 0.67. The measurement was verified first on the Fv panel,
  where trastuzumab's own predicted interface recovers 19 of 19 1N8Z residues.
- **Thermostability**: ESM-2 650M mean per-residue score relative to trastuzumab's scFv
  (`thermo.csv`). A foldability proxy, not a melting temperature, and the two formats differ in
  length, so the signal is the ordering among the designs.

## Selection rule, in code

`scripts/denovo_select.py`, in this order: the design must have passed both gates; it must refold
within 10 Å of the pose it was designed in (a good `min_ipae` 40 Å away is a different binding site,
not ours); the epitope overlap must be at least 0.2; then rank by the averaged binding percentile;
then at most one design per CDR3 cluster, so we do not submit six variants of one idea.

## The six designs

| id | CDR3 | min_ipae b2 / of3 | epitope b2 | pLDDT | RMSD to designed pose | ESM-2 vs trastuzumab | identity to trastuzumab | CDR3 distance to nearest known |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vhh-01 | AMGYGGSREDGTVY | 2.62 / 1.53 | 0.27 | 92 | 2.7 Å | -0.051 | 0.658 | 11 |
| vhh-02 | NASRPHVYAANY | 5.63 / 15.33 | 0.50 | 88 | 4.5 Å | -0.084 | 0.683 | 10 |
| vhh-03 | AARKTGTSSTDFDY | 8.09 / 7.38 | 0.33 | 94 | 8.2 Å | +0.010 | 0.675 | 10 |
| vhh-04 | AAALSSNFYNARLHTASKYDY | 8.91 / 8.88 | 0.56 | 87 | 7.7 Å | -0.226 | 0.648 | 16 |
| vhh-05 | ARYDRSYSGGSAHSISAYDY | 0.51 / 18.83 | 0.67 | 97 | 1.4 Å | -0.185 | 0.654 | 14 |
| vhh-06 | AARKSYDYSSETSGFDY | 11.29 / 15.56 | 0.43 | 90 | 9.5 Å | -0.012 | 0.659 | 13 |

**vhh-01** is the one both folding models agree on: the best averaged percentile (98.8), a tight
refold (2.7 Å), and it passes the organisers' oracle under OpenFold3. **vhh-05** has the best single
boltz2 interface score and the largest epitope overlap, but the two models disagree sharply
(0.51 against 18.83), so its evidence is one model's, not two. If only four are made, we would send
vhh-01, vhh-05, vhh-04 and vhh-02: two that the models agree on and two that cover different CDR3
clusters and the widest epitope overlap.

## Three bugs found on the way, each of which would have produced a wrong submission

1. **Designs from different runs collided.** Every run numbers its examples from zero, so
   `her2_vhh_denovo_0` exists in all three. The merged refold table keyed on that name alone, 120
   designs became 40, and the top-ranked "design" was two different molecules mixed together. Fixed
   by putting the run id in the key; `tests/test_denovo_vhh_rank.py` holds it.
2. **Every epitope overlap was 0.00.** The refold queries carry no `residue_ids`, so HER2 is
   numbered from 1 in the predictions, and nothing matched 1N8Z. Fixed with an explicit offset
   (chain C = 488, chain D = 590), checked against the chain lengths.
3. **The gate subset matched on a prefix**, so `her2_vhh_denovo_3` also let `..._39` through and
   designs that had never been gated were folded. Fixed to exact matching.

## What we did not do

- **No lab data.** Every ranking is a prediction.
- **Identity to the framework's own parent.** The official checker checks CDRH3s and identity to trastuzumab. A
  VHH on the caplacizumab framework is far from trastuzumab, but its framework is still a known
  antibody's. If the lab's database measures identity over the whole sequence against every known
  antibody, that is a risk we carry and name rather than hide.
- **One seed per design in the refold.** Replicates would give a spread per design; we have the
  spread between two models instead.

## Files

| File | What it holds |
| --- | --- |
| `denovo_submission.csv` | the submission table, checked by the official checker's `team-validate` (6/6 valid) |
| `vhh/candidates_all.csv`, `vhh/candidates_pass.csv` | all 120 designs with their numbers, and the 118 that passed the gate |
| `vhh/scores.csv` | refold metrics per model, percentiles, epitope and the combined binding axis |
| `vhh/selection.csv` | the six picks with the reason each was picked |
| `vhh/epitope_*.csv` | HER2 interface residues per predicted complex |
| `thermo.csv`, `esm2/` | the ESM-2 scores and the raw platform result |
| `fv/` | the secondary Fv designs and why they are not de novo |

