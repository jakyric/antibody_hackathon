# S4 decisions log (yield / Tm / developability)

Times are CET, 14 Nov 2026. Evidence tags: **M** = measured here, **L** = literature/tool documentation, **I** = my inference.

## D1. One in-process TNP runner instead of the `TNP --file --ncores` CLI  (14:00)
* `TNP` shells out to `NanoBodyBuilder2 --n_threads 1` once per sequence on CPU: the earlier session measured about 4 min
  per sequence (M, `tools/README.md`).  For 100+ designs that is hours.
* `tnp_batch.py` loads NanoBodyBuilder2 once, models on the GPU (about 2 s per sequence on the local RTX 4070, M), strips
  hydrogens exactly as TNP does and then calls TNP's own functions (`CDR_Assigner.main`, `main_compactness`,
  `CreateAnnotation`, `assign_flag`) unchanged from `tools/TNP`.  7 reference nanobodies: 21 s end to end (M).
* Rejected: reusing OF3/RFD4 structures for the patch metrics.  TNP's patch calculator expects an IMGT-numbered heavy chain
  `H` as written by NanoBodyBuilder2; our complexes are not numbered that way and are co-folded with HER2, so the patch
  scores would also be influenced by the antigen.  A fresh apo NanoBodyBuilder2 model is the same input TNP's calibration
  used (L), and costs 2-3 s.
* Consequence (M): TNP patch metrics are model-dependent.  The same 7 sequences modelled on the VM versus locally give
  CDR3 compactness within 0.03, but PSH differs by up to 15 units (caplacizumab 102.9 vs 118.3) and PPC can flip
  (caplacizumab 0.353 vs 0.000).  Treat `tnp_psh/ppc/pnc` as +-15% noisy; use the red/amber counts, not single values.

## D2. Candidate-set tiering  (14:05, user instruction)
* Cheap for every design: `liab_*` (pure Python) and `nm_tm` (NanoMelt, about 1 s per sequence, M).
* Slow tools (TNP, monoyield) only for: controls + reference baselines + top 35% of `of3_consensus_rank_score` in
  `score_of3.csv` + S2 `designs_fwdiv.csv` variants whose parent is in that top set.  Until S3 publishes the score column
  only controls and baselines get the slow tools.  Direction of the score is assumed higher = better
  (`--score-lower-better` flips it).
* monoyield loads ESMC-6B once per fold checkpoint (5 folds x 2 models = 10 loads, about 1 min each, M), so a call costs
  about 10 min regardless of batch size.  Hence ONE batched VM call per cycle, at most one per 10 min, bf16.

## D2b. Gate and cancelled generators  (14:20, orchestrator)
* TNP and monoyield only for controls + designs with `denovo_vm_pass` true in `rules.csv` that also sit in the top 35% of the gate-passers by
  `of3_consensus_rank_score`.  `fwdiv` (framework diversification cancelled), `boltzgen`, `boltzgen-cys` (user stopped BoltzGen) are ignored.
* NanoMelt runs on the CPU only (user decision after a CUDA OOM at 248 sequences on the 8 GB laptop GPU): 8 threads, batches of 32, a failing
  batch is bisected so a single bad sequence gets `nm_error` instead of killing the batch.  TNP's NanoBodyBuilder2 uses the local GPU but peaks at 0.57 GiB
  (M); `--tnp-cpu` forces the CPU if that is still too much.
* The 248-sequence NanoMelt pass at 14:18-14:24 ran on the GPU in chunks of 24 (before the CPU decision reached me); those values are the same model
  and the same weights, so they are valid, only the device differed.
* Every cycle skips rows without `vhh_seq`: `designs_boltzgen_metrics.csv` also matches the `designs_*.csv` glob and crashed one pass (KeyError).

## D2c. monoyield paused  (14:33, user rule via orchestrator)
* The user forbids starting anything on the shared VM GPU without explicit permission.  monoyield (ESMC-6B) exists only on the VM, so it is paused.
  `my_yield`, `my_yield_sd`, `my_mdratio`, `my_mdratio_sd`, `my_dtype` are partial: 50 rows (7 reference nanobodies, 14 controls, 29 RFD4 designs of the first
  candidate batch).  Later designs are blank, NOT zero; `tools_done` shows which tools really ran.
* Guards: `props.py` forces `--where local` and drops monoyield unless `--allow-vm`; `run_vm.sh` exits 5 unless `S4_ALLOW_VM=1`.  Resume = the user says so, then
  `python3 denovo/s4_props/props.py --allow-vm --watch 240` (monoyield: one batched VM call per cycle, at most one per 10 min, bf16; about 2 min for 29 sequences).
* Measured before the pause: bf16 vs fp32 on the 7 reference nanobodies differ by at most 0.010 (yield) and 0.007 (mdratio) log10 units against an ensemble sd
  of 0.05-0.28, so bf16 is adequate; a monoyield call took 116-127 s for 21-29 sequences at bf16.
* The yield prior that remains available without the VM is the nine-rule score (`nb9_*`), which is weaker (D4b).

## D3. NanoMelt and TNP also on the VM  (14:00, user instruction)
* One micromamba env `s4` under `~/moritz/s4` (python 3.12; anarci/hmmer from bioconda, pip torch, ImmuneBuilder, the
  NanoMelt pip deps pinned to the local versions).  NanoMelt: VM and local Tm identical to 0.01 C on 7 sequences (M).
  TNP: see D1 for the model-dependence.  The VM copy is a fallback and a cross-check (`--where both`); local stays primary
  because the local GPU is idle while the VM GPU is shared.
* micromamba itself was fetched from the official `micro.mamba.pm` endpoint into `~/moritz/s4/bin` (no sudo, nothing outside
  `~/moritz`).  NanoMelt downloaded about 8 GB of ESM weights (ESM-1b 7.8 GB, ESM-2 150M 0.6 GB) into `~/moritz/s4/torch_cache`.

## D4. What the three numbers can and cannot say  (honest limits)
* NanoMelt: apparent Tm from sequence, MAE about 4-7 C on held-out nanobodies (L, paper).  The model was trained on
  natural-like nanobodies; de novo CDR3s and diversified frameworks are out of distribution.  Use it to drop clear
  outliers and to compare designs of one framework, not to claim an absolute Tm.  Licence CC BY-NC-SA 4.0
  (non-commercial): fine for the hackathon, matters for the "path to value" story.
* The hackathon judges thermostability against *measured trastuzumab*.  Trastuzumab is a Fab/IgG; a VHH Tm is not
  comparable like for like (different assay construct, no VL/CH1).  NanoMelt cannot evaluate a Fab.  Our trastuzumab-VH
  row in `baselines.csv` is a sanity check only (a VH without VL is not a stable domain, and has no VHH hallmarks).
* TNP flags were calibrated on clinical-stage nanobodies (L).  They flag developability outliers; they are not Tm
  predictors.
* monoyield: E. coli monomeric yield model (ESMC-6B + LoRA, 5-fold ensemble).  Adaptyv expresses cell-free, not in
  E. coli, so treat as a proxy for "does it express and stay monomeric".  bf16 inference is a deviation from the authors'
  fp32; see README for the measured bf16-vs-fp32 difference.
* Liability motifs are textbook heuristics.  Thresholds sit just above the range of the six well-behaved reference
  nanobodies (hydrophobic 5-mer KD >= 2.5, K+R in CDRs >= 6).  Caplacizumab (approved) carries a CDR3 DG, so a flag is a
  prompt to look, not a rejection.

## D4b. Nine-rule CDR2/CDR3 score (monoyield paper)
* Implemented exactly as in the paper's Methods ("Nanobody sequence rules") and the orchestrator's brief; unit test = parent hNbBCII10-FGLA must give
  6/9 (fails N1=V, N2=R, CDR2 pos 6 = L), which the code reproduces.  Extra: the six-rule version needs C1 = Y only (5/6 for the parent).
* Not applicable (not counted either way): the four CDR2 rules when the IMGT CDR2 is not 8 residues, the five CDR3 rules when the loop is < 5 residues.
  `nb9_applicable` tells how many rules were really evaluated; weight the score accordingly.
* Evidence level: the three added rules and the broadened C1 rule were derived on the same dataset (authors: "independent validation is needed"),
  the scaffold is a single fixed framework with CDR3 of 5-13 residues.  RFD4 CDR3s here are longer.  Treat as a weak yield prior.

## D5. Rejected / skipped
* SoluProt / NetSolP / Aggrescan3D: not attempted (time box; the three requested tools came first).  Adaptyv's own guidance
  names them, so it is a known gap in S4, listed in the README.
* Original TAP is web-only (see memory/tools README); not part of S4 (the user asked for the Therapeutic Nanobody Profiler).
* Running RFD4-Proteina or any docker image on the VM: forbidden by the user, not done.
* Operational: the shared VM's sshd intermittently resets new connections ("kex_exchange_identification: Connection reset") under load.  `run_vm.sh` and
  `vm_install.sh` retry up to 10 times on ssh exit 255 and multiplex over one ControlMaster connection.  Never edit `run_vm.sh` while a pass is running (bash
  reads scripts incrementally; that produced the transient `line 26: syntax error` at 14:16).
