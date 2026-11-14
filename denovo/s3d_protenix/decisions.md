# S3d decisions (Protenix-v2 on Nebius)
* 14:40 user asked to run protenix-v2 on Nebius in addition to ESMFold-2 / OpenDDE for the top designs. VM not used.
* Schema facts (forge.py schema protenix-v2 --full): input is ONE native Protenix JSON list (<= 256 MB, no job-count limit stated);
  `msa_mode` is const `none`; `sample_count` 1-16; `model_seeds` up to 16 unique. So: no MSA, one submission of 42 jobs.
* Parameters chosen: sample_count 3, model_seeds [7] (three diffusion samples of one seed; median over samples is the design value,
  spread is kept). More seeds were not worth the queue time.
* Selection: Pareto front 1 (38 designs) union the Pareto top-20, capped at 40 -> 38 designs; plus 4 controls. Rationale: the
  user wants the extra predictor for "our top designs"; controls calibrate it (trastuzumab Fv positive, framework-only decoys).
* Submitted 14:40:17, run `op-70cf0433a9fbb3f9`, status queued (Nebius says about 20 min for a first call).
* Gotcha: forge.py stores run records as nebius-runs/<label>.json, so an input file named like the label breaks `submit`
  (AttributeError: 'list' has no 'get'). Input file is therefore `nebius-runs/ptxtop-input.json`, label `ptx-top`.
* 14:58 FIRST RUN FAILED (op-70cf0433a9fbb3f9, one 42-job list): the model computed all 42 jobs (seed 7 completed in 404 s, last
  job `ctl-fwrand-03` succeeded) but the platform's post-check raised `confidence result count must equal the exact seed/sample
  product` (expects 1 job x 1 seed x 3 samples = 3 results per run), so no result was published (`result_available: False`).
  The schema calls the input a "list of jobs", but the validator only accepts ONE job per run. Not resumable.
* 15:00 D3: resubmitted as one run per complex (same parameters, sample_count 3, seed 7), 16 runs: 13 designs (Pareto order,
  capped at 13 to stay polite on a shared platform instead of 40 runs) + trastuzumab Fv + 2 decoys. Run labels ptx01..ptx16
  (`nebius-runs/ptxNN.json` hold the run ids), inputs `nebius-runs/ptxinNN.json`, manifest.csv maps label <-> design_id.
  `finish.py` waits, fetches into results/ptx-runs/<label> and parses after each new batch.
* README commands for the single-list route (build_input.py) are kept only as a record of what failed; use build_single.py +
  submit_all.py.
* 15:06-15:12 all 16 runs succeeded (about 6-12 min each once running; 48 structures, no PAE in the output, only summary confidence:
  pTM, ipTM, ranking_score, chain-pair ipTM, pLDDT). Nebius names artifacts by UUID; `parse_ptx.py` maps them back by sha256 from result.json.
* RESULT (n = 13 designs, 3 controls, 3 samples each, median): trastuzumab Fv control ipTM 0.68 / pTM 0.78, but its VH-HER2
  chain-pair ipTM is only 0.36 and its HER2 epitope overlap 0.16 (the crystal pose is NOT recovered without MSA). Decoys
  ctl-fwrand-01 / -02: ipTM 0.47 / 0.22. Designs: median ipTM 0.20, max 0.32; 0 of 13 above the best decoy. Spearman of ptx_iptm
  vs S3's OpenFold3 / OpenDDE-ESMFold2 consensus: -0.2 to +0.09 (n = 13, so noisy). Conclusion: Protenix-v2 without MSA on Nebius
  does not separate designs from decoys and a framework-only decoy outscores every design; treat `ptx_score` as a weak, noisy
  extra vote, not as a gate. Not evidence of binding either way.
* `ptx_score` = percentile of ptx_iptm among the 13 designs only (no PAE in the output), so it is a ranking inside this small set.
