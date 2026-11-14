# S3d: Protenix-v2 (Nebius) as an extra orthogonal predictor

Adds `protenix-v2` on the Nebius Scientific AI platform to the OpenDDE / ESMFold-2 orthogonal set, for the top designs only.
Nothing runs on the team VM; JapanFold is not used.

## Run
```bash
PY=denovo/.venv/bin/python
$PY denovo/s3d_protenix/build_input.py          # nebius-runs/ptxtop-input.json (38 designs + 4 controls) + manifest.csv
python3 ~/.claude/skills/nebius-scientific-ai/forge.py submit protenix-v2 \
  --entry 'protenix-input|protenix-input-json/v1|nebius-runs/ptxtop-input.json|application/json' \
  --params '{"checkpoint":"protenix-v2","msa_mode":"none","sample_count":3,"model_seeds":[7]}' --label ptx-top --no-wait
nohup bash denovo/s3d_protenix/finish.sh > denovo/s3d_protenix/finish.log 2>&1 &   # waits, fetches, parses
# by hand: forge.py fetch op-70cf0433a9fbb3f9 --dest results/ptx-top && $PY denovo/s3d_protenix/parse_ptx.py
```
Run id `op-70cf0433a9fbb3f9` (recorded in `nebius-runs/ptx-top.json`; re-running the submit follows the existing run, it never pays twice).

## Input
One native Protenix JSON list, 42 jobs: chain 1 = VHH (`vhh_seq`), chain 2 = HER2 domain IV (119 aa, query of `denovo/targets/her2_domIV.a3m`),
`count: 1`, name = design_id. Controls: trastuzumab Fv (VH, VL, HER2) as positive; three framework-only decoys (`ctl-fwrand-01..03`).
Designs = Pareto front 1 union the Pareto top-20 (gate-passing population by construction), ranked by pareto_rank, maximin, composite, capped at 40.

## Output: `results/denovo/score_ptx.csv`
One row per design_id, `ptx_*` = median over the 3 samples (`ptx_iptm`, `ptx_ptm`, `ptx_ranking_score`, `ptx_chain_pair_iptm`, binder/HER2 pLDDT,
`ptx_epi_overlap`, `ptx_cdr_fv_frac4`, `ptx_clash25`, ... geometry from S3's `geom.py`), `ptx_iptm_best`, `ptx_iptm_sd`, `ptx_n_samples`,
`ptx_score` (0..1: mean percentile of iptm and, if PAE is in the output, of min interface PAE; reference population = non-control designs),
`is_control`, `stage_ts`. `report.md` has controls, designs vs decoys and Spearman vs S3's OpenFold3 / OpenDDE / ESMFold-2.

## Caveats
The Nebius profile is **no-MSA only** (schema: `msa_mode=none`), so it is a single-sequence fold of both chains; Protenix-v2 here is the
upstream v2.0.0 checkpoint. Scores are filters, not binding evidence; calibrate against the decoys and the trastuzumab control (report.md).
