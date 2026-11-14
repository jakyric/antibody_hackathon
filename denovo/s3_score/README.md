# S3: refolding and orthogonal binding prediction

Scores VHH + HER2 domain IV complexes with OpenFold3 (primary) and orthogonal predictors, all on JapanFold, and writes
interface metrics per design. Scores are **filters, not binding evidence**: the literature (and our own calibration below)
shows confidence scores rank binders weakly, so every gate is set from decoys and controls scored by the same code.

## Run

```bash
export JAPANFOLD_API_KEY=...                                  # never printed or logged
PY=denovo/.venv/bin/python
$PY denovo/s3_score/make_calib.py                              # controls (results/denovo/designs_control.csv) + Fv calibration set
$PY denovo/s3_score/make_vhh_calib.py                          # + real VHH-antigen pairs from Smorodina et al. (needs data/smorodina/)
$PY denovo/s3_score/run_s3.py calib --models of3,dde,esm2      # score the calibration complexes
$PY denovo/s3_score/calibrate.py                               # AUCs, per-model gates (gates.json), calibration.md
$PY denovo/s3_score/run_s3.py watch --models of3,dde,esm2      # loop: score every new designs_*.csv row, rebuild the CSVs
$PY denovo/s3_score/run_s3.py score --models af2ig             # optional extra: AF2-IG on the design (or OF3) pose
$PY denovo/s3_score/run_s3.py assemble                         # rebuild score_of3.csv / score_orth.csv from raw metrics
```

Everything is resumable: a JapanFold job is named by a hash of (model, params, target ids) and its id is cached in
`results/denovo/jobs/`, so re-running never pays twice; parsed metrics are append-only in `results/denovo/s3_raw/metrics.jsonl`;
downloaded structures sit in `results/denovo/s3_raw/<model>/<job id>/`.

## Continuous operation (what actually runs during the hackathon)

```bash
export JAPANFOLD_API_KEY=...            # in the shell that starts it; never printed
denovo/s3_score/start_watch.sh          # nohup background loop, pid in denovo/s3_score/watch.pid, log denovo/s3_score/watch.log
denovo/s3_score/stop_watch.sh           # stop now; a restart adopts every unparsed JapanFold job, so nothing is lost or paid twice
touch denovo/STOP                       # global stop file: the loop stops launching, drains in-flight jobs, exits
```

The loop wakes every 90 s, re-reads `designs_*.csv` and `rules.csv`, and launches what is missing; it stops launching at 15:40 local time or at
`denovo/STOP`, or when the JapanFold spend reaches the $25 stage budget. Each cycle rewrites `score_of3.csv` and `score_orth.csv`
atomically and logs a 6-line funnel + control canary.

* **Tier 1**: OpenFold3 (+MSA, 3 samples, seed 0) on every design whose `rules.csv` row has `denovo_vm_pass` and `canon_cys_ok` true (plus all controls);
  designs without a verdict yet wait. Order is round-robin over (framework x hotspot set), newest first, so the first batch is representative.
* **Tier 2**: designs in the top ~35% of the OF3 consensus so far, and all controls: OpenDDE-abag, ESMFold-2, a second OF3 seed (`of3s1`) and the
  lysozyme decoy for specificity.
* **AF2-IG companion** (default on, `--tier2 dde,esm2,of3s1,af2ig`): re-scores the generator's own pose (`design_cif`: target + binder coordinates, no MSA, deterministic; the nearest JapanFold
  equivalent of folding with the design as a template) for every gate-passer, and the OF3 pose for controls; metrics `af2ig_iptm`, `af2ig_interface_pae`, `af2ig_min_ipae`, `af2ig_pose_rmsd_vhh` (binder RMSD to the input after superposing on the target; the
  upstream RFD4 oracle passes at <= 2 A). Lowest priority in the job queue. Boltz-2 / RF3 exist as extras (`--models bz2,rf3`).
* `of3_recipe` documents the OF3 setup in every row: MSA on, no template. Job concurrency is capped (`--max-jobs 8`) so the shared account queue (HTTP 429 when full) stays free for the other stages.
* Ignored generators: `fwdiv` (cancelled), `boltzgen*` (stopped by the user); rows of `designs_boltzgen_metrics.csv`-style files are dropped.

## What is run

| key | JapanFold model | params | role |
|---|---|---|---|
| `of3` | `openfold3` | MSA server on, 3 diffusion samples, seed 0, `write_pae` | primary (best separator of measured trastuzumab binders in the team's earlier refolds) |
| `dde` | `opendde-abag` | MSA on, 3 samples, seed 0 | orthogonal: antibody-antigen checkpoint |
| `ddeg` | `opendde` | same | orthogonal: general complex checkpoint |
| `ddenm` | `opendde-abag` | MSA off | orthogonal: does the MSA matter |
| `esm2` | `esmfold2` | MSA off, 3 samples | orthogonal: language-model folder, folds the target unguided |
| `af2ig` | `af2ig` | `write_pae` | extra: AlphaFold2 initial guess re-scores the generator's own pose (or the OF3 pose) |
| `bz2`, `rf3` | `boltz2`, `rf3` | MSA on | extras, low weight (Boltz-2 is over-confident) |

Input per design: chain A = VHH, chain B = HER2 domain IV (119 aa, UniProt P04626 489-607). Each design is also folded against
hen egg-white lysozyme (decoy antigen) so every design carries its own specificity control. The trastuzumab Fv control is a
3-chain complex (A = VH, B = VL, C = HER2).

## Metrics (per design and model; column prefix = model key)

* `iptm`, `ptm`, `complex_plddt`, `confidence`, `iptm_sd` (spread over the 3 samples), `binder_plddt`, `cdr3_plddt`
* `min_ipae`, `mean_ipae`: cross-block PAE (VHH rows x HER2 columns and the transpose), in Angstrom
* `ipae_contacts`: mean PAE over residue pairs in contact (heavy atoms < 4.5 A); `ipae_cdr_epitope`: mean PAE between the CDR residues and the trastuzumab epitope residues
* `ipsae`: Dunbrack 2025 ipSAE, d0res variant: for each VHH residue i, mean of 1/(1+(PAE_ij/d0)^2) over HER2 residues j with PAE_ij < 10 A, d0 = 1.24 (n-15)^(1/3) - 1.8 with n = number of those residues (n >= 27, d0 >= 1); maximum over i; symmetrised by taking the larger of the two directions. PAE comes from the top-ranked sample's `_pae.npz`.
* `n_contact_pairs`, `n_vhh_res_contact`, `n_her2_res_contact`
* `epi_overlap`: share of the HER2 residues touched by the VHH that lie in trastuzumab's epitope (19 residues within 4.5 A of the Fv in 1N8Z, `targets/epitope.json`); `epi_coverage`: share of the epitope touched; `cov_hs*`: coverage of each hotspot set
* competition proxy: the predicted HER2 is superposed (CA, 109 shared residues) onto the 1N8Z HER2 chain and the VHH is moved with it; `clash25` / `clash30` = VHH heavy atoms closer than 2.5 / 3.0 A to a trastuzumab Fv atom (true contacts are >= 3 A, so these are steric impossibilities), `fv_contact40` = atoms within 4 A, `cdr_fv_frac4` = share of CDR atoms within 4 A of the Fv. A VHH that binds the trastuzumab epitope must overlap the Fv; a VHH binding elsewhere does not.
* pose stability over the 3 samples: `pose_rmsd_max` (VHH CA RMSD in the HER2 frame vs sample 0), `cdr3_cent_max` (CDR3 centroid shift), `epi_overlap_min`
* `own_vhh_rmsd_her2frame`, `own_vhh_rmsd_fold`: refold fidelity versus the generator's own design (when `design_cif` is set)
* `<model>_hewl_iptm/ipsae/min_ipae` and `<model>_spec_iptm_delta`: the same design against lysozyme and the difference
* cross-model: `x_<m>_vs_of3_vhh_rmsd` / `_cdr3_dist`, `xmodel_vhh_rmsd_mean`, `xmodel_cdr3_dist_max`, `xmodel_n_agree` (models whose CDR3 sits within 8 A of OF3's)
* `<model>_pass` (past the calibrated gate AND epitope overlap >= 0.5), `<model>_pct` (percentile rank among candidates, 1 = best), `consensus_rank_score` (mean percentile over gated models), `n_models_pass`

`score_of3.csv` holds the `of3_*` columns; `score_orth.csv` everything else plus `of3_pass`/`of3_pct`. Controls (`ctl-*`, `is_control = 1`) are
scored alongside the designs and excluded from the percentile ranks.
BoltzGen's own metrics (`designs_boltzgen_metrics.csv`) are joined as `bg_*` columns when present: informational, never a gate.

## Calibration (the point of this stage)

`make_calib.py` / `make_vhh_calib.py` build 100 complexes with known labels, scored exactly like designs:

* trastuzumab Fv variants with measured SPR binders (KD <= 10 nM) vs "no reportable KD" (12 + 12), and FACS high vs low (16 + 16), 3-6 CDRH3 changes
* crystal VHH-antigen pairs (10 cognate) vs the same VHHs against other antigens and against the shuffled pairs that AlphaFold3 scored highest (Smorodina et al. 2026, public)
* per design: the lysozyme decoy; 12 framework-only VHHs (BCII10-FGLA + random CDR3) against HER2

`calibrate.py` reports AUC with bootstrap intervals per family/model/metric (`calibration.md`), picks per model the metric with the best AUC on the
VHH families (ipSAE within 0.02 preferred), and sets the gate at the 95th percentile of the decoy values (`gates.json`). It also reports behaviour
**per generator** (RFD4 vs BoltzGen vs framework-diversified), because they may differ.

See `decisions.md` for the results and what we trust.
