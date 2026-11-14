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
   │ 6 score       ipSAE and iLIS from PAE (only for designs that pass the contact pre-filter)
   │ 7 rank        filter → rank by iLIS (ipSAE tie-break, ≤3 per CDR-H3) → JSON / FASTA / CSV
   └ 8 figures, 9 deck   optional: PyMOL figures and a PowerPoint summary
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
python run_pipeline.py --stages figures deck                # optional
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
* **ipSAE / iLIS** are re-implementations (Dunbrack d0res form, PAE < 10 Å; iLIS = √(LIS·cLIS), PAE < 12 Å, C-beta
  contacts < 8 Å). They have not been checked against the authors' code. In the controls both put the wild-type
  prediction in the non-binder range, so they read as "docks like the base loop", not as affinity.
* **Ridge ranking aid.** `ridge_sum` is an additive ridge model on log10 K_D of the 180 SPR designs (alpha 10); it is
  informative but weak on new combinations (Spearman ~0.36) and is not used for selection.

## Results of the 2026-11-14 run (`results/2026-11-14/`)

| Step | Count |
|---|---|
| raw designs (3–4 BLOSUM-45 changes at positions 2,3,4,8,9,10) | 23,760 |
| after no-Cys, liabilities ≤ WT, ≥3 from all known sequences | 5,247 |
| shortlist × Vernier genotypes | 36 H3 → 528 sequences |
| pass abguard | 517 (35 H3) |
| folded (+4 controls) / docked (ipTM ≥ 0.9) | 521 / 201 |
| keep the crystal contacts | 38 designs, 16 H3 |
| final list | 20 designs, 11 distinct H3; 11 keep the R98–D108 salt bridge |

Controls (single seed):

| Control | ipTM | ipSAE | iLIS | Fv–HER2 pairs kept | CDR-H3-vs-rest pairs kept |
|---|---|---|---|---|---|
| base loop `WGSQGFYAYD` | 0.947 | 0.938 | 0.647 | 0.85 | 0.91 |
| trastuzumab | 0.687 | 0.021 | 0.103 | 0.15 | 0.82 |
| non-binder 1 | 0.655 | 0.016 | 0.035 | 0.15 | 0.77 |
| non-binder 2 | 0.622 | 0.000 | 0.006 | 0.00 | 0.68 |

Top 10 (full list with all scores and VH/VL sequences: `top_candidates.json`, `.fasta`, `top_mutations.csv`):

| # | CDR-H3 | mutations vs WT | Vernier | iLIS | ipSAE | R98–pos10 salt bridge |
|---|---|---|---|---|---|---|
| 1 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L46L>M | 0.743 | 0.940 | yes |
| 2 | `WGSSGFYAWY` | G3S D4S M9W D10Y | H93S>T | 0.737 | 0.944 | no |
| 3 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L49Y>F | 0.734 | 0.948 | yes |
| 4 | `WGSSGFYAWY` | G3S D4S M9W D10Y | – | 0.731 | 0.941 | no |
| 5 | `WGSSGFYAWY` | G3S D4S M9W D10Y | H93S>A | 0.729 | 0.938 | no |
| 6 | `WSSQGFYGRD` | G2S G3S D4Q A8G M9R | H47W>F | 0.724 | 0.926 | yes |
| 7 | `WDSQGFYAWR` | G2D G3S D4Q M9W D10R | H93S>N | 0.719 | 0.934 | no |
| 8 | `WDYQGFYYYD` | G2D G3Y D4Q A8Y M9Y | L49Y>H | 0.718 | 0.949 | yes |
| 9 | `WYYQGFYAYN` | G2Y G3Y D4Q M9Y D10N | L46L>I | 0.712 | 0.927 | no |
| 10 | `WGRQGFYADS` | G3R D4Q M9D D10S | H47W>F | 0.710 | 0.936 | no |

What the run says about the Vernier zone: single swaps barely move ipTM on average (median change −0.006 to +0.003 per
site), but for 9 of the 19 top designs with a swap the same CDR-H3 with WT Vernier does **not** dock (ipTM 0.63–0.69)
and the swap restores docking (0.94–0.95).

### Caveats

* One seed per complex; the top of the list is within noise. ipTM jumps between ≈0.6 (not docked) and ≈0.95 (docked).
* No metric here predicts affinity. In the SPR data no loop with a non-acidic residue at position 10 beat wild type, and
  in the FACS screen Ser at position 10 is the least enriched residue; 9 of the 20 designs change D10 and lose the
  R98–D108 salt bridge. Prefer the 11 that keep it unless you want to test that risk.
* `abguard` is the organisers' demo build; the jury's checker may differ. Its CDR-3 string includes the H94 anchor, so a
  swap at H94 can make an otherwise failing loop pass.
* The base loop has one SPR replicate.

## Repository layout

```
run_pipeline.py            stage runner
pipeline/                  s1_enumerate … s9_deck, japanfold.py (client), common.py, pymol_prep.py, pymol_figures.py
config/pipeline.example.yaml
reference/her2_fv_1n8z.cif trastuzumab Fv + HER2 domain IV (PDB 1N8Z)
results/2026-11-14/        outputs of the run above (no structures or PAE matrices)
tests/test_basic.py
```

## Credits

Position choice and BLOSUM-45 recipe: Moritz Ertelt's `moritz` branch, after Chinery et al., *Simple computational
methods can outperform deep learning in designing diverse, binder-enriched antibody libraries* (bioRxiv
10.1101/2024.03.26.586756). ipSAE: Dunbrack. iLIS: Kim et al. Contacts: GetContacts (Stanford).
