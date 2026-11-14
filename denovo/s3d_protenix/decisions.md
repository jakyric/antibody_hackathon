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
