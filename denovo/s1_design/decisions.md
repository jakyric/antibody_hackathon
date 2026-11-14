# S1a decisions log (JapanFold `rfd4`, VHH on HER2 domain IV)

Times are local CET, 14 Nov 2026. Evidence first, choice second, what was rejected last.

## D1. Hosted rfd4 CAN do framework-conditioned VHH design (answer to the 25-minute discovery question)
Evidence (probe jobs `probe-bcii10-v1/v3`, results in `results/denovo/jobs/`):
* `POST /v1/designs` with `protocol=rfd4-binder`, `structure` (mmCIF text), `contig`, `hotspots` (residue level) and
  `params {num_designs<=8, recipe, nsteps<=60, seed}` is accepted. The OpenAPI only documents the rfd3 contig, but the rfd4
  wrapper uses the same grammar as the RFD4 VM docs: kept ranges `B1-22` are fixed, bare ranges `13-17` are designed
  with a length drawn from the seed, `/0` starts a new chain.
* The archive contains `inputs/design.cif`, the annotated CIF the wrapper builds from our contig. Reading it shows what is
  actually conditioned: kept residues get `condition_sequence = residue` (sequence fixed), backbone atoms N,CA,C,O of ALL kept
  residues (HER2 and framework) are in ONE coordinate group (`condition_CoordinateGroup = 1`), each designed block becomes one
  placeholder residue with `expsegmin/max` = the requested length range, hotspot residues get `condition_hotspot = True`.
  There is no `condition_cath-label` column (the VM caplacizumab example has `2.60.40.x`; the hosted wrapper cannot set it).
* Designs come back with the framework residues unchanged (checked on every design in `geom_rfd4jf.csv`, `ok=1` means the
  designed chain is exactly FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4) and `motif_rmsd` 0.38-0.44 A.
* Checkpoint: not disclosed anywhere in the API. The CDR3s that come back look antibody-like (`...DY` endings, aromatic/Ser/Gly
  rich) and contact the epitope, so it is at least antibody-capable; I cannot tell antibody-2026-10-30 from the general release.
  Reported as unknown.

Consequence 1: the framework pose relative to HER2 is NOT sampled. Kept HER2 + framework backbone are one rigid group
(motif RMSD ~0.4 A). The pose is whatever we submit. We submit the framework superposed on trastuzumab's VH in 1N8Z
(CA of FR positions by IMGT number, RMSD 0.68 A for bcii10, 0.71 A for capla). That choice is what makes the designs
epitope-directed (they sit where the trastuzumab VH sits, so a VHH that binds there should compete with trastuzumab).
It also means every design shares one binding mode, so pose `p1` (10-14 degree rotation about the epitope centroid plus
1.5 A shift away, rejected if any framework CA is closer than 4 A to HER2) is included as a second mode and compared in
the feedback table below.

Consequence 2: HER2 has a gap (581-590 unresolved in 1N8Z). The API refuses two kept stretches of one chain with nothing
between them (400 with a clear message) and rfd4 itself fails with a bond-length error if the numbering is contiguous across the
gap (`probe-bcii10-v1`). Solution: HER2 as two chains, A 489-580 and C 591-607, like the crystal file (C and D).

## D2. Frameworks
* `bcii10` = hNbBCII10-FGLA, RCSB 3EAK chain A ("NbBCII10 humanized (FGLA mutant)", 1.95 A, verified title via the RCSB API).
  Chosen because RFantibody and Germinal both fixed it after expression testing (research report), it has a human-like FR2
  (F37 G44 L45 A47, so no camelid hallmark liability for solubility but also lower hydrophilicity), and a real crystal structure.
* `capla` = caplacizumab VHH framework from the RFD4 VM example `Data/03-vhh/pmvk_caplacizumab.cif` (CDRs removed). Chosen as the
  second, structurally different, clinically validated VHH framework with camelid FR2 (R44? E/R hallmark), so the pipeline
  compares a humanised and a camelid-hallmark framework.
* Not used: any framework for which I had no structure I could verify from RCSB in the time box.
* IMGT/CDR boundaries computed with `immunum` (baseline/.venv), confidence 0.81 (bcii10) and 0.79 (capla, FR pieces numbered with
  stub CDRs because the example has no CDR residues). CDR = IMGT 24-40 / 56-65 / 105-117 (extended, as the RFD4 examples).

## D3. Campaign design
* Grid = framework (2) x hotspot set (3: hs1_core, hs2_broad, hs3_alt_patch from `targets/epitope.json`) x CDR3 length
  range (c3a 12-16, c3b 16-21) x pose (p0, p1) x repeats; one job = 8 designs (API max), seed unique per job.
* `recipe=high/front_loaded`, `nsteps=60` (the API says fewer steps is a different, worse model; high/front_loaded is the
  upstream setting for antibody CDRs per the model blurb and VM docs). Not varied.
* CDR1 13-17 (bcii10) / 12-15 (capla), CDR2 8-10 / 9-11: native lengths +-.
* Cost measured: $0.0034 for a 2-design job (47 billable s) so about $0.002 per design. The bottleneck is the queue, not money.
* VM route (Foundry docker `antibody-2026-10-30`) is NOT used: the user ruled it out (coordinator message), and it is not needed.

## D4. Known weaknesses (stated honestly)
* Every design keeps the framework 100% identical to its parent by construction. Measured on the first 96 bcii10 designs
  (`summary.py`): 65-71% of each design's residues are inherited framework, and global-alignment identity to the parent 3EAK
  sequence is 74-82% (mean 0.79); ALL designs are at or above 70% to their own parent. The 70% rule is about the NEAREST known
  antibody, and germline-like VHH frameworks (FR1 `QVQLVESGGGLVQPGGSLRLSC`, FR4 `WGQGTLVTVS`) are shared by many VHHs, so the
  nearest-neighbour identity can be higher still. This stage cannot fix that (framework is held fixed); S2 must diversify it
  (about 10+ framework substitutions, hallmarks/core/Cys held) and S3 must re-fold the diversified sequences.
* Epitope = trastuzumab's contact residues in the 1N8Z crystal (19 residues). A design that binds there will compete with
  trastuzumab by construction; whether it binds at all is the question for S3 and then the assay.
* rfd4 returns no confidence score. Geometry columns in `geom_rfd4jf.csv` (contacts, clashes, hotspot contacts, epitope
  precision) are only a cheap pre-filter, not evidence of binding.
* One rigid binding mode per pose: all `p0` designs sit exactly where the trastuzumab VH sits. If a VHH cannot bind in that mode
  (no VL, different loop geometry), no CDR design fixes it. `p1` is a weak hedge against that.

## D5. Wave 2: what changed after the first 96 designs were scored (evidence from S2, via the coordinator)
* S2's VM `abguard` run on wave-1 bcii10 designs: de novo class for 53%, and the pass rate tracks the CDR3 block length (12-13 aa: 0%, 14: 21%,
  15-17: 50-62%, 18-21: 80-100%; total VHH >= 128 aa: ~90%). The gate computes identity to the trastuzumab VH as 1 - Levenshtein / max(Fv length, 120)
  (`denovo/s2_rules/gate_drivers.md`), so for bcii10 the 70% line sits at a VHH length of ~125. The caplacizumab-like framework passes at any length
  (63/64 in wave 1). Decision: new waves use only the long CDR3 block, total 17-22 aa, for both frameworks, and capla gets real weight.
* The monoyield paper (Fig. 5, same hNbBCII10-FGLA scaffold): CDR3 packs on the invariant C104-A105-A106 and F42, which drives monomeric yield. In wave 1 my designed
  block started at IMGT 105, so 48% of designs changed A106 (AR/AK/AV/AL). Decision: wave 2 keeps C104-A105-A106 FIXED (added to FR3 in the contig, the designed block
  follows); `caa_fixed=1` in `params_json`. "Keep total length": the `cdr3` column (A105 A106 + designed residues) is 17-22 aa, so the contig's designed range is 15-20 and wave-1 /
  wave-2 `cdr3` strings are comparable. (If the intent was 17-22 designed residues AFTER the fixed AA, the total would be 19-24; say so and re-run `c3n` with `CDR3_N = (17, 22)`.)
* capla: the example structure has no CDR3 residues, so A105/A106 cannot be kept from "its own" structure. They are grafted from bcii10 after superposing the four residues that end FR3
  (RMSD 0.09 A over 16 backbone atoms, peptide bond C104-A105 1.31 A vs 1.33 ideal) and written as ALA 76, 77; FR4 renumbered +3. An engineering assumption, flagged `caa_grafted` in `params_json`.
* Wave 1 (uniform grid, 24 cells) was stopped at 216 designs; its in-flight jobs were collected, nothing deleted. All of its designs stay in the CSV (`caa_fixed` absent = wave 1).

## D6. How the pipeline feeds back into design (`feedback.py`; allocation and reason for every wave in `plan_log.jsonl`)
* Success criterion (coordinator): `denovo_vm_pass` in `rules.csv`. Quality among successes: S3's OF3 results. I did NOT use S3's `consensus_rank_score`: no rfd4 design clears S3's gates
  (95th percentile of decoys; OF3 ipSAE >= 0.45, OpenDDE min_ipAE <= 0.43), so it is the constant 0.464 for 93% of designs. The quality term is instead the mean percentile (among our scored designs)
  of OF3 `iptm`, `-min_ipae` and trastuzumab-epitope overlap. S3's calibration gives these metrics AUC of only about 0.65 on real VHH-antigen pairs, so the factor is shrunk (4 pseudo-observations),
  clamped to 0.5-2, and used only after 40 passing designs were scored.
* Honest reading: OF3 gives our designs iptm 0.06-0.38 (median 0.12) and ipSAE ~0 (max 0.014), the same level as the framework-only decoy controls (median iptm 0.14) and far below the
  trastuzumab Fv control (0.555). OF3 therefore does NOT say these VHHs bind HER2; it only mildly orders them. The settings table is steering evidence, not binding evidence.
* Allocation: weight = shrunk VM pass rate x OF3 quality factor, floor 0.25, zero for clear losers (none reached the 10-scored-pass threshold with a significant deficit). First OF3-informed plan (14:45):
  capla hs1_core p0 got 4 jobs, capla hs3 p0 3, bcii10 hs1_core / hs2_broad cells 1 each (quality 0.34-0.40 vs 0.49 overall).

### Settings table at the end of the campaign (cells = framework x hotspot set x CDR3 label x pose; `denovo/s1_design/feedback.csv`, `cell_table.md`)
CDR3 labels: `c3a`/`c3b` wave 1 (block starts at IMGT 105), `c3n` designed 15-20 after the fixed CAA, `c3m` (capla 11-17, bcii10 14-18), `c3m2` (capla 10-16, bcii10 13-16; the final range, IMGT block 12-18 / 15-18). "VM de novo pass" = `denovo_vm_pass` in S2's rules.csv.
OF3 columns are S3's `of3_iptm` / `of3_epi_overlap` means; quality = mean percentile among our scored designs (iptm high, min_ipae low, epitope overlap high).

| framework | hotspots | CDR3 | pose | designs | VM de novo pass | OF3 iptm (mean) | OF3 epitope overlap | quality among passes (n) |
|---|---|---|---|---|---|---|---|---|
| bcii10 | hs1_core | c3a | p0 | 16 | 2/9 | 0.16 | 0.21 | 0.57 (2) |
| bcii10 | hs1_core | c3a | p1 | 16 | 2/9 | 0.14 | 0.25 | 0.68 (2) |
| bcii10 | hs1_core | c3b | p0 | 16 | 9/9 | 0.10 | 0.08 | 0.29 (9) |
| bcii10 | hs1_core | c3b | p1 | 16 | 6/9 | 0.15 | 0.28 | 0.44 (6) |
| bcii10 | hs1_core | c3m | p0 | 7 | 4/7 | 0.13 | 0.11 | 0.36 (4) |
| bcii10 | hs1_core | c3m | p1 | 10 | 6/10 | 0.11 | 0.15 | 0.32 (6) |
| bcii10 | hs1_core | c3m2 | p0 | 40 | 28/40 | 0.10 | 0.12 | 0.32 (16) |
| bcii10 | hs1_core | c3m2 | p1 | 40 | 22/40 | 0.14 | 0.18 | 0.45 (10) |
| bcii10 | hs1_core | c3m2L | p0 | 8 | 7/8 | - | - | - (0) |
| bcii10 | hs1_core | c3m2L | p1 | 8 | 7/8 | - | - | - (0) |
| bcii10 | hs1_core | c3n | p0 | 103 | 33/40 | 0.12 | 0.13 | 0.38 (30) |
| bcii10 | hs1_core | c3n | p1 | 94 | 35/38 | 0.13 | 0.19 | 0.37 (35) |
| bcii10 | hs2_broad | c3a | p0 | 16 | 3/9 | 0.15 | 0.29 | 0.67 (3) |
| bcii10 | hs2_broad | c3a | p1 | 16 | 1/9 | 0.16 | 0.26 | 0.56 (1) |
| bcii10 | hs2_broad | c3b | p0 | 16 | 8/9 | 0.15 | 0.21 | 0.57 (8) |
| bcii10 | hs2_broad | c3b | p1 | 16 | 6/9 | 0.12 | 0.16 | 0.47 (6) |
| bcii10 | hs2_broad | c3m | p0 | 10 | 7/10 | 0.11 | 0.22 | 0.37 (7) |
| bcii10 | hs2_broad | c3m | p1 | 13 | 12/13 | 0.10 | 0.28 | 0.42 (12) |
| bcii10 | hs2_broad | c3m2 | p0 | 48 | 33/48 | 0.16 | 0.23 | 0.53 (16) |
| bcii10 | hs2_broad | c3m2 | p1 | 32 | 13/32 | 0.17 | 0.24 | 0.55 (10) |
| bcii10 | hs2_broad | c3m2L | p0 | 8 | 5/8 | - | - | - (0) |
| bcii10 | hs2_broad | c3n | p0 | 94 | 36/38 | 0.12 | 0.21 | 0.43 (32) |
| bcii10 | hs2_broad | c3n | p1 | 73 | 27/31 | 0.11 | 0.23 | 0.39 (27) |
| bcii10 | hs3_alt_patch | c3a | p0 | 16 | 5/9 | 0.16 | 0.30 | 0.63 (5) |
| bcii10 | hs3_alt_patch | c3a | p1 | 16 | 2/9 | 0.16 | 0.25 | 0.46 (2) |
| bcii10 | hs3_alt_patch | c3b | p0 | 16 | 6/9 | 0.18 | 0.24 | 0.56 (6) |
| bcii10 | hs3_alt_patch | c3b | p1 | 16 | 8/9 | 0.17 | 0.28 | 0.58 (8) |
| bcii10 | hs3_alt_patch | c3m | p0 | 8 | 7/8 | 0.14 | 0.22 | 0.52 (7) |
| bcii10 | hs3_alt_patch | c3m | p1 | 8 | 5/8 | 0.15 | 0.39 | 0.57 (5) |
| bcii10 | hs3_alt_patch | c3m2 | p0 | 40 | 28/40 | 0.17 | 0.22 | 0.59 (7) |
| bcii10 | hs3_alt_patch | c3m2 | p1 | 40 | 26/40 | 0.20 | 0.29 | 0.68 (9) |
| bcii10 | hs3_alt_patch | c3m2L | p0 | 8 | 6/8 | - | - | - (0) |
| bcii10 | hs3_alt_patch | c3m2L | p1 | 8 | - | - | - | - (0) |
| bcii10 | hs3_alt_patch | c3n | p0 | 150 | 50/59 | 0.16 | 0.22 | 0.53 (43) |
| bcii10 | hs3_alt_patch | c3n | p1 | 125 | 38/48 | 0.19 | 0.27 | 0.58 (29) |
| capla | hs1_core | c3a | p0 | 16 | 9/9 | 0.12 | 0.24 | 0.50 (9) |
| capla | hs1_core | c3a | p1 | 16 | 9/9 | 0.10 | 0.22 | 0.34 (9) |
| capla | hs1_core | c3b | p0 | 16 | 9/9 | 0.14 | 0.11 | 0.39 (9) |
| capla | hs1_core | c3b | p1 | 8 | 8/8 | 0.09 | 0.27 | 0.30 (8) |
| capla | hs1_core | c3m | p0 | 6 | 6/6 | 0.13 | 0.20 | 0.51 (6) |
| capla | hs1_core | c3m | p1 | 8 | 8/8 | 0.09 | 0.28 | 0.37 (8) |
| capla | hs1_core | c3m2 | p0 | 56 | 55/56 | 0.14 | 0.29 | 0.51 (16) |
| capla | hs1_core | c3m2 | p1 | 48 | 47/48 | 0.11 | 0.20 | 0.36 (16) |
| capla | hs1_core | c3m2L | p0 | 8 | 8/8 | - | - | - (0) |
| capla | hs1_core | c3n | p0 | 159 | 54/54 | 0.14 | 0.25 | 0.53 (47) |
| capla | hs1_core | c3n | p1 | 90 | 34/34 | 0.11 | 0.23 | 0.39 (23) |
| capla | hs2_broad | c3a | p0 | 8 | 7/8 | 0.16 | 0.15 | 0.52 (7) |
| capla | hs2_broad | c3a | p1 | 8 | 7/8 | 0.16 | 0.26 | 0.55 (7) |
| capla | hs2_broad | c3b | p0 | 8 | 8/8 | 0.13 | 0.30 | 0.51 (8) |
| capla | hs2_broad | c3b | p1 | 8 | 8/8 | 0.12 | 0.25 | 0.47 (8) |
| capla | hs2_broad | c3m | p0 | 7 | 6/7 | 0.19 | 0.12 | 0.53 (4) |
| capla | hs2_broad | c3m | p1 | 6 | 6/6 | 0.17 | 0.19 | 0.57 (6) |
| capla | hs2_broad | c3m2 | p0 | 56 | 52/56 | 0.16 | 0.18 | 0.55 (22) |
| capla | hs2_broad | c3m2 | p1 | 56 | 46/48 | 0.23 | 0.23 | 0.62 (14) |
| capla | hs2_broad | c3m2L | p0 | 8 | 8/8 | - | - | - (0) |
| capla | hs2_broad | c3m2L | p1 | 8 | - | - | - | - (0) |
| capla | hs2_broad | c3n | p0 | 137 | 53/53 | 0.14 | 0.21 | 0.56 (47) |
| capla | hs2_broad | c3n | p1 | 119 | 42/42 | 0.14 | 0.22 | 0.49 (31) |
| capla | hs3_alt_patch | c3a | p0 | 8 | 8/8 | 0.11 | 0.11 | 0.38 (8) |
| capla | hs3_alt_patch | c3a | p1 | 8 | 8/8 | 0.20 | 0.34 | 0.72 (8) |
| capla | hs3_alt_patch | c3b | p0 | 8 | 8/8 | 0.12 | 0.32 | 0.58 (8) |
| capla | hs3_alt_patch | c3b | p1 | 8 | 8/8 | 0.12 | 0.17 | 0.43 (8) |
| capla | hs3_alt_patch | c3m | p0 | 6 | 6/6 | 0.16 | 0.38 | 0.68 (6) |
| capla | hs3_alt_patch | c3m | p1 | 7 | 7/7 | 0.17 | 0.25 | 0.59 (7) |
| capla | hs3_alt_patch | c3m2 | p0 | 48 | 42/48 | 0.17 | 0.26 | 0.64 (15) |
| capla | hs3_alt_patch | c3m2 | p1 | 72 | 71/72 | 0.20 | 0.25 | 0.63 (16) |
| capla | hs3_alt_patch | c3m2L | p1 | 8 | 8/8 | - | - | - (0) |
| capla | hs3_alt_patch | c3n | p0 | 136 | 45/46 | 0.18 | 0.26 | 0.57 (36) |
| capla | hs3_alt_patch | c3n | p1 | 130 | 46/46 | 0.15 | 0.25 | 0.54 (46) |

Reading it (weak evidence, small n per cell, all OF3 numbers near decoy level): capla passes the VM de novo gate in 98% of designs at any length, bcii10 only in 71% (63% inside the final 12-18 block); `p0` (trastuzumab VH pose) gives
higher OF3 quality than `p1` in most cells; the capla cells (quality 0.5-0.6 among passes) sit above the bcii10 cells (0.25-0.55), with bcii10 hs2_broad p1 the lowest. The planner followed this: the 14 logged plans move from uniform to
capla-heavy (last plan: capla 15 of 24 jobs, bcii10 9). None of these differences is large enough to call a winner; the one robust effect is the gate pass rate.


## D7. Incidents
* `HTTP 402 ... runs still in progress can cost up to $90.72`: the account's in-flight runs reserve their maximum cost, so with S3's 20-40-target OpenFold3 batches in flight a new design job is refused until some
  finish. Not an empty wallet (balance was $90). Fixed by retrying with back-off (a real out-of-credit 402 has different text and still fails the job).
* Head-of-line blocking in my first runner: as many threads as in-flight jobs, so jobs that had already `succeeded` waited behind queued ones for ~15 minutes. Fixed with 40 threads, a semaphore of 12 for NEW
  submissions (`--max-inflight`, the coordinator's cap) and `recover_jobs()`, which collects any submitted job that never reached the CSV (restarting the runner is safe).
* **Duplicate designs (found 14:55, my bug).** A job with seed S makes its 8 designs with seeds S..S+7, and I gave consecutive repeats consecutive seeds (S, S+1, ...), so each new repeat re-made 7 of the
  previous job's 8 designs: at 14:55 only 422 of 1360 rows had a unique `vhh_seq`. Effects: every wasted row was scored by S2/S3 and would have taken a ranking slot. Fix: (1) `load_rows()` keeps the first row per
  `vhh_seq` and moves the rest to `duplicates_removed.csv` (the CSV shrank from 1360 to ~430 rows; design ids of removed rows disappear, their CIFs stay in `results/denovo/designs/`); (2) new jobs are named
  `s1a3-*` with seeds `200000 + cell*1000 + rep*8` (disjoint); (3) `run_wave` drops any new row whose sequence is already present. Cost of the bug: about $1 and ~25 minutes of queue time.
* `pkill -f` with the runner's command line also killed my own shell once; harmless.

## D8. Loop-aware batch (user decision 15:20): modelled HER2 including the unresolved 581-590 loop
* Evidence (coordinator, `denovo/validation/gap_check.py`): the 10-residue gap sits inside the epitope (epitope residue 591 is 2.6 A from the gap centre) and 22% of crystal-input designs have VHH atoms within 3 A of the loop anchors.
* Input: `denovo/targets/her2_domIV_modeled.cif` (AlphaFold DB model of domain IV superposed on the 1N8Z HER2, ONE contiguous chain A, 119 aa, resi k = HER2 488+k) plus the SAME placed framework chain B as before (same crystal frame, so p0/p1 stay valid).
  Contig `A1-119,/0,B1-22,...` (kept sequence + backbone, one rigid group), hotspot n -> `A(n-488)`. The API accepted it (the earlier refusal was only for two kept stretches in one chain).
* Smoke (1 job, 8 designs, 0.007 USD): 2 chains, framework intact, motif_rmsd 0.47-0.62 A, 144-276 contacts, 5-6 of 6 hotspots touched. Epitope precision is lower (0.50-0.65 vs 0.8-0.9) because the designs also contact the modelled loop residues, which are not in the 19-residue epitope list.
* Batch: 9 jobs = 72 designs (both frameworks, three hotspot sets, both poses, c3m2 ranges, `caa_fixed=1`, `her2_input: modeled_loop` in `params_json`), generator `rfd4jf-loop`, ids `<fw>-rfd4jfloop-s<seed>-<k>`, seeds 500000+, job names `s1l-*`.
  Integration rule from the coordinator: the output CIF has 2 chains, so the CIFs are in `results/denovo/designs_loop/` and `design_cif` is EMPTY in the CSV (the path is in `params_json.design_cif_loop`); rows share the dedup-by-sequence rule.
* Caveat: the modelled loop (pLDDT ~64) is a guess; "VHH does not clash with the modelled loop" is weak evidence only.

## D9. End state (15:25)
* `results/denovo/designs_rfd4jf.csv`: 1480 rows, 1480 unique `vhh_seq` (capla 782, bcii10 698; 1408 `rfd4jf`, 72 `rfd4jf-loop`). 1092 rows have an IMGT-block CDR3 (incl. AA) of 12-18 aa, the leaderboard window; 388 rows from earlier waves have 19-22 and were left in the CSV
  (the leaderboard drops them). Per unique cell see the table above; per framework x hotspot x pose: 94-149 designs in each of the 12 cells.
* Removed or never written: 986 duplicate-sequence rows (`duplicates_removed.csv`, the seed-overlap bug), 147 designs outside 12-18 collected after the cap changed (`out_of_range_cdr3.csv`).
* Spend on rfd4: 3.72 USD by the job ledger (about 0.0025 USD per kept design including the wasted duplicates); all 1480 designs, 0 failed jobs except 7 early 402s that were never submitted (retried by name later).
* Budget was never the limit: the shared queue (fleet at 260/260 chips) and my own seed bug were.
