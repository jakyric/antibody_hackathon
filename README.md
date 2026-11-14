# HER2 / trastuzumab: CDR-H3 and Vernier-zone redesign pipeline

Branch `jakob` of the Oslo HER2 hackathon repository. It turns the exploratory analysis of the hackathon into one
reproducible pipeline that proposes **novel trastuzumab variants** (CDR-H3 changes plus compensating Vernier-zone
swaps), folds them, and ranks them with several structure-based metrics. The earlier exploratory scripts
(`research_*.py`, old panels) were removed; they stay in the git history.

> Everything here is computational. No score below predicts affinity; the folds only say whether a design docks like
> a known binder and keeps the crystal contacts.

## What the pipeline does

```
best tested CDR-H3 (WGSQGFYAYD, SPR 0.023 nM; trastuzumab 0.114 nM)
   │ 1 enumerate   BLOSUM-45 recipe, 3-4 changes at the 6 non-contact positions (2,3,4,8,9,10)
   │               no Cys · liabilities ≤ WT · ≥3 Hamming from all 525,063 known sequences
   │ 2 panel       diverse shortlist (≥3 apart) × Vernier genotypes (BLOSUM-45 ≥ +1, only where a Vernier
   │               residue touches a changed H3 position: H94 H93 H47 L36 L46 L49)
   │ 3 validate    organisers' abguard (edit-distance novelty, liabilities, Cys parity)
   │ 4 fold        JapanFold OpenDDE antibody-antigen, + 4 controls (WT, base loop, 2 measured non-binders)
   │ 5 contacts    GetContacts on every fold: Fv–HER2 interface and CDR-H3-vs-rest, compared with 1N8Z
   │ 5b pose       RMSD to the crystal after fitting the Fv framework: HER2 epitope RMSD, HER2 RMSD, CDR-H3 RMSD
   │ 6 score       ipSAE and iLIS from PAE for every fold (score.scope: filtered limits it to designs passing the contacts)
   │ 7 rank        epitope RMSD → contact filter → salt-bridge rule → rank by iLIS (ipSAE tie-break, ≤3 per CDR-H3) → JSON / FASTA / CSV
   └ 8 figures, 9 deck, 10 plots   optional: PyMOL figures, a PowerPoint summary, score distributions with trastuzumab marked
```

Stages are decoupled by files in `paths.out_dir`, so any stage can be re-run alone.

## Quick start

```bash
git clone https://github.com/jakyric/antibody_hackathon && cd antibody_hackathon && git checkout jakob
python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
cp config/pipeline.example.yaml config/pipeline.yaml        # edit paths and hosts; the file is git-ignored
export JAPANFOLD_API_KEY=...                                # https://japanfold.aiand.com/account/token (never commit it)
python run_pipeline.py --list
python run_pipeline.py --stages enumerate panel             # data-only stages, no network
python run_pipeline.py                                      # stages 1-7
python run_pipeline.py --stages figures plots deck          # optional
python tests/test_basic.py                                  # data-free unit tests
```

### External pieces (not in this repository)

| Needed by | What | How |
|---|---|---|
| enumerate, panel | hackathon data: `trastuzumab_her2_spr_all_768.csv`, `tz_her2_facs_3classes_524346.csv.zip` | set `paths.data_dir` |
| validate | the organisers' `abguard` build (it calls itself a demo; the jury's may differ) | `validate.mode: ssh \| local \| skip` |
| fold, score | a JapanFold account and `JAPANFOLD_API_KEY` | set `fold.job_id` to re-use an earlier job instead of folding again |
| contacts | [GetContacts](https://github.com/getcontacts/getcontacts) and `vmd-python` | `conda env create -p gc/env -f environment-contacts.yml`; `git clone https://github.com/getcontacts/getcontacts gc/getcontacts` |
| contacts, figures | PyMOL on the PATH (adds hydrogens, renders figures) | `contacts.pymol` |
| deck | `python-pptx` and a `.pptx` template | `paths.template_pptx` |

## Method notes

* **Start from the best tested loop, not wild type.** `WGSQGFYAYD` (G3S, D4Q, M9Y) is the tightest SPR-measured CDR-H3
  (0.023 nM, one replicate; trastuzumab 0.114 nM, eight). Mutations are counted from it; every design is still ≥3 edits
  from all known loops, as `abguard` requires.
* **H3 restriction** follows Moritz's `moritz` branch: positions that touch HER2 (1, 5, 6, 7) are fixed, substitutions
  use the paper's BLOSUM-45 recipe, P(b|a) ∝ f_b·exp(0.25·BLOSUM45(a,b)) with CDR-H3 background frequencies, top 6 per
  position. Here the recipe is applied to the base-loop residue, and D10 is allowed on request.
* **Vernier zone.** The Vernier residues that touch CDR-H3 in 1N8Z are H94 (R98), H93, H47, L36, L46, L49. A residue is
  varied only if it touches a position that changed from the base loop, one swap at a time, BLOSUM-45 ≥ +1, no Cys/Pro.
  The datasets carry no Vernier variation, so nothing in the data scores these choices.
* **Fold model.** In this project's controls only OpenDDE antibody-antigen put both measured binders above both
  non-binders (ipTM 0.947 and 0.687 vs 0.655 and 0.622). OpenFold3 with MSA was flat (~0.5) and ESMFold2-fast scored
  the non-binders above trastuzumab, so OpenDDE alone is used.
* **GetContacts, two views, residue-pair level, vs the crystal.** *Interface*: Fv vs HER2 (27 crystal pairs, 6 with
  CDR-H3). *Intra*: CDR-H3 (99–108) vs the rest of the Fv (22 pairs), which covers the Vernier and VL contacts and the
  R98–D108 salt bridge. HER2 residues 581–590 are unresolved in the crystal and ignored when counting new contacts.
  The base loop's own prediction sets the pass level (85% of interface pairs, 1 H3–HER2 pair missing).
* **Same epitope (pose RMSD).** Every prediction is superposed on the 1N8Z crystal with the C-alpha atoms of the Fv
  framework (VH without CDR-H3 + VL). `epitope_rmsd` is then the C-alpha RMSD of the 19 HER2 residues that touch the
  crystal Fv; `her2_rmsd` covers all resolved HER2 residues and `h3_rmsd` the loop. The filter is epitope RMSD ≤ 3 Å
  (the best tested loop's prediction is 1.3 Å; trastuzumab's own prediction is 30 Å away because it does not dock).
  Every design with ipTM ≥ 0.9 is within 3 Å (median 1.2 Å), so in this run docking and binding the crystal epitope
  coincide.
* **ipSAE / iLIS** are re-implementations (Dunbrack d0res form, PAE < 10 Å; iLIS = √(LIS·cLIS), PAE < 12 Å, C-beta
  contacts < 8 Å). They have not been checked against the authors' code. In the controls both put the wild-type
  prediction in the non-binder range, so they read as "docks like the base loop", not as affinity.
* **Salt-bridge rule (hard filter on the final list).** The R98 (H94) - position-10 salt bridge is the prominent
  CDR-H3-to-framework contact. A design is kept if GetContacts finds that salt bridge in its fold; if it does not, the
  design is kept only when **both** partners are mutated (an H94 swap *and* a changed position 10); a design with just
  one mutated partner and no bridge is excluded. Keeping both partners mutated is a compensation hypothesis that no
  data here can confirm.
* **Plots.** `fig_distributions.png` shows iLIS, ipSAE and the contact counts (Fv-HER2 and CDR-H3-vs-rest residue
  pairs) for all folded designs; `fig_ipsae_ilis.png` the ipSAE-iLIS plane. Trastuzumab is marked as the prediction of
  its sequence (same model and pipeline as the designs) and, for contact counts, as the crystal value (27 and 22 pairs).
* **Ridge ranking aid.** `ridge_sum` is an additive ridge model on log10 K_D of the 180 SPR designs (alpha 10); it is
  informative but weak on new combinations (Spearman ~0.36) and is not used for selection.

## Results of the 2026-11-14 run (`results/2026-11-14/`)

| Step | Count |
|---|---|
| raw designs (3–4 BLOSUM-45 changes at positions 2,3,4,8,9,10) | 23,760 |
| after no-Cys, liabilities ≤ WT, ≥3 from all known sequences | 5,247 |
| shortlist × Vernier genotypes | 36 H3 → 528 sequences |
| pass abguard | 517 (35 H3) |
| folded (+4 controls) | 521 |
| docked (ipTM ≥ 0.9) and within 3 Å of the crystal epitope | 201 |
| keep the crystal contacts | 38 designs, 16 H3 |
| salt-bridge rule | 24 designs, 11 H3 |
| final list | 20 designs, 11 distinct H3; 17 with the R98–pos 10 salt bridge, 3 with both partners mutated |

Controls (single seed):

| Control | ipTM | ipSAE | iLIS | Fv–HER2 pairs kept | CDR-H3-vs-rest pairs kept | epitope RMSD (Å) |
|---|---|---|---|---|---|---|
| base loop `WGSQGFYAYD` | 0.947 | 0.938 | 0.647 | 0.85 | 0.91 | 1.3 |
| trastuzumab | 0.687 | 0.021 | 0.103 | 0.15 | 0.82 | 29.8 |
| non-binder 1 | 0.655 | 0.016 | 0.035 | 0.15 | 0.77 | 11.2 |
| non-binder 2 | 0.622 | 0.000 | 0.006 | 0.00 | 0.68 | 33.5 |

Top 10 (full list with all scores and VH/VL sequences: `top_candidates.json`, `.fasta`, `top_mutations.csv`):

| # | CDR-H3 | mutations vs WT | Vernier | iLIS | ipSAE | epitope RMSD (Å) | R98–pos 10 |
|---|---|---|---|---|---|---|---|
| 1 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L46L>M | 0.743 | 0.940 | 1.44 | salt bridge present |
| 2 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L49Y>F | 0.734 | 0.948 | 1.51 | salt bridge present |
| 3 | `WSSQGFYGRD` | G2S G3S D4Q A8G M9R | H47W>F | 0.724 | 0.926 | 1.18 | salt bridge present |
| 4 | `WGSSGFYAWY` | G3S D4S M9W D10Y | H94R>K | 0.723 | 0.941 | 0.63 | both partners mutated |
| 5 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L49Y>H | 0.718 | 0.949 | 1.49 | salt bridge present |
| 6 | `WGYQGFYGDD` | G3Y D4Q A8G M9D | H93S>T | 0.680 | 0.943 | 1.26 | salt bridge present |
| 7 | `WGDQGFYYYY` | G3D D4Q A8Y M9Y D10Y | H94R>Q | 0.678 | 0.932 | 0.93 | both partners mutated |
| 8 | `WGYQGFYGDD` | G3Y D4Q A8G M9D | L46L>F | 0.676 | 0.933 | 1.49 | salt bridge present |
| 9 | `WDSQGFYAWR` | G2D G3S D4Q M9W D10R | H94R>Q | 0.675 | 0.924 | 1.39 | both partners mutated |
| 10 | `WGSAGFYGWD` | G3S D4A A8G M9W | L49Y>H | 0.670 | 0.942 | 1.3 | salt bridge present |

What the run says about the Vernier zone: single swaps barely move ipTM on average (median change −0.006 to +0.003 per site), but for 9 of the 20 final designs with a swap the same CDR-H3 with WT Vernier does **not** dock (ipTM below 0.9) and the swap restores docking.

### Caveats

* One seed per complex; the top of the list is within noise. ipTM jumps between ≈0.6 (not docked) and ≈0.95 (docked).
* No metric here predicts affinity. In the SPR data no loop with a non-acidic residue at position 10 beat wild type, and
  in the FACS screen Ser at position 10 is the least enriched residue. The salt-bridge rule keeps 3 designs without the
  bridge (D10Y with H94 K or Q, D10R with H94Q), whose partners are both mutated; that compensation is a hypothesis.
* `abguard` is the organisers' demo build; the jury's checker may differ. Its CDR-3 string includes the H94 anchor, so a
  swap at H94 can make an otherwise failing loop pass.
* The base loop has one SPR replicate.

## Repository layout

```
run_pipeline.py            stage runner
pipeline/                  s1_enumerate … s5b_pose … s10_plots, japanfold.py (client), common.py, pymol_prep.py, pymol_figures.py
config/pipeline.example.yaml
reference/her2_fv_1n8z.cif trastuzumab Fv + HER2 domain IV (PDB 1N8Z)
results/2026-11-14/        outputs of the run above (no structures or PAE matrices)
tests/test_basic.py
```

## Credits

Position choice and BLOSUM-45 recipe: Moritz Ertelt's `moritz` branch, after Chinery et al., *Simple computational
methods can outperform deep learning in designing diverse, binder-enriched antibody libraries* (bioRxiv
10.1101/2024.03.26.586756). ipSAE: Dunbrack. iLIS: Kim et al. Contacts: GetContacts (Stanford).
