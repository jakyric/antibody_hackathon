# `moritz` branch: CDRH3 diversification and refolding screen (HER2 / trastuzumab)

Optimisation-track idea: keep trastuzumab's binding mode but change at least four CDRH3 positions (the minimum for
CDRH3 identity below 70% under immunum's CDR3 definitions), choosing the positions from data and the substitutions
with the paper's own BLOSUM-45 recipe, then check the designs by refolding the Fv-HER2 complex on three services.

## Headline results

- **Baseline binding classifiers** trained on the 524k-variant FACS screen reach AUC 0.995 (gradient-boosted trees)
  and 0.997 (the authors' CNN, re-implemented in numpy) on the authors' clonotype test split. On the paper's SPR
  designs they barely beat "fewer mutations from trastuzumab", but those designs were pre-screened by a CNN from the
  same lineage, so that comparison is range-restricted.
- **At least 4 CDRH3 changes are needed** if `abguard` uses immunum's CDR3 boundaries (13 residues IMGT, 11 Kabat):
  3 changes give 77% / 73% identity (fail), 4 give 69% / 64% (pass). The real `abguard` was not available, so
  `abguard-lite` (immunum numbering, length, cysteines, CDR3 identity) approximates it.
- **Positions**: only W99, G103, F104, Y105 of the 10 CDRH3 residues touch HER2 in 1N8Z. Choosing positions by the
  worse of the FACS and SPR "does changing this hurt" estimates gives **{2,3,4,8,9}** with contacts excluded and
  **{2,4,5,8,9}** with contacts enabled (G103 enters, G101 leaves).
- **Refolding** (Fv + HER2 domain IV, scored against the crystal): only **OpenFold3 (MSA on)** separated measured binders
  from measured non-binders (trastuzumab 3/3 seeds, both non-binders 0/3). Boltz-2 reproduces the crystal pose for
  everything; OpenDDE antibody-antigen fails trastuzumab itself. Pose recovery is weak evidence, not proof of binding.
- **Final picks** (OpenFold3 recovers the crystal pose in all seeds):

| Variant | Rank | CDRH3 | Changes (structure numbering) | Services recovering the pose |
|---|---|---|---|---|
| contacts excluded | 1 | `WAYYGFYVMD` | G100A G101Y D102Y A106V | OpenFold3, Boltz-2, OpenDDE |
| | 2 | `WDYSGFYYMD` | G100D G101Y D102S A106Y | all three |
| | 3 | `WAYDGFYYVD` | G100A G101Y A106Y M107V | all three |
| | 4 | `WYAYGFYYMD` | G100Y G101A D102Y A106Y | all three |
| | 5 | `WYSYGFYARD` | G100Y G101S D102Y M107R | OpenFold3, Boltz-2 |
| | 6 | `WYGAGFYYRD` | G100Y D102A A106Y M107R | OpenFold3, Boltz-2 |
| | 7 | `WDADGFYYGD` | G100D G101A A106Y M107G | OpenFold3, Boltz-2 |
| contacts enabled | 1 | `WDGAAFYYMD` | G100D D102A G103A A106Y | all three |
| | 2 | `WYGYSFYYMD` | G100Y D102Y G103S A106Y | all three |
| | 3 | `WYGAGFYYRD` | (same as above) | OpenFold3, Boltz-2 |

Enabling the contact residues gave fewer robust designs (3 of 12 against 7 of 12 pass OpenFold3 in every seed), but
two of the three tolerate a G103 change in the prediction. Neither variant addresses thermostability, and changing
only non-contact CDRH3 positions can at best tie trastuzumab's affinity.

## Caveats and open questions for the organisers

- No lab data: everything is computational. A recovered pose does not show binding.
- M107: all three measured M107 singles in the paper's SPR data (K, T, A) were 40 to 70 times weaker, so designs that
  change M107 carry affinity risk.
- Developability: the picks carry 4 to 6 aromatic residues against trastuzumab's 3 (polyreactivity risk).
- Confirm with the organisers: the exact `abguard` identity and mutation-count rule, and whether an affinity tie
  counts as "comparable" for the tiers.

## What is in the branch

```
baseline/            scripts, README and outputs (see baseline/README.md)
  common.py structure_tools.py            shared helpers; comparison with the 1N8Z crystal
  train_eval.py score_cdrh3.py            binding baselines and the gate
  blosum_diversity.py                     position choice + BLOSUM-45 candidates (--allow-contacts)
  refold_prepare.py refold_japanfold.py refold_nebius.sh refold_analyze.py refold_consensus.py pick_shortlist.py
  outputs/                                shortlists, metrics, consensus (…_contacts = contacts-enabled variant)
  refold_inputs*/                         HER2 domain IV sequence + alignment, JapanFold input YAMLs
nebius-runs/         Boltz-2 and MSA run records and inputs
results/             refolded structures and scores (top-ranked model per complex; PAE matrices and zips omitted)
```

Variants: `VARIANT=contacts` makes every script write `..._contacts` files and folders without touching the first run.

## Reproduce

```bash
uv venv baseline/.venv --python 3.12 && uv pip install --python baseline/.venv/bin/python -r requirements.txt
cd baseline
.venv/bin/python train_eval.py                                               # classifiers (needs data below)
.venv/bin/python blosum_diversity.py                                         # contacts excluded
VARIANT=contacts .venv/bin/python blosum_diversity.py --allow-contacts       # contacts enabled
.venv/bin/python refold_prepare.py && .venv/bin/python refold_japanfold.py --model openfold3 --msa on --tag jf_of3_msa
.venv/bin/python refold_analyze.py && .venv/bin/python refold_consensus.py && .venv/bin/python pick_shortlist.py
```

JapanFold calls read `JAPANFOLD_API_KEY` from the environment; Nebius calls need `NEBIUS_CREDS` and the
`nebius-scientific-ai` skill. No keys are stored in this repository.

## Data (not included, fetch it yourself)

| Path expected | Source |
|---|---|
| `data/her2/tz_her2_facs_3classes_524346.csv.zip` (unzip to `data/her2/facs_extract/`), `data/her2/trastuzumab_her2_spr_all_768.csv` | hackathon bucket `data/her2-trastuzumab/` |
| `data/her2/her2_fv_1n8z.cif` | hackathon bucket `rfd4-proteina/docs/Data/12-her2/her2_fv_1n8z.cif` |
| `data/mehta/{train,val,test}.csv`, `saved_weights.h5`, `shanehsazzadeh_test.csv` | Zenodo 10831512 `data.zip`, folder `clonotype_split/0.7_0.15_0.15` |
| `data/paper/spr_submitted.csv` | github.com/oxpig/Tz_her2_affinity_and_beyond, `data/spr/` |

Paper: Chinery, Hummer, Mehta et al., "Simple computational methods can outperform deep learning in designing diverse,
binder-enriched antibody libraries", bioRxiv doi:10.1101/2024.03.26.586756.
