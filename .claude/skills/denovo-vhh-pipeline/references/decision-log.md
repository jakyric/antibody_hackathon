# Decision log (chronological, condensed; evidence in each stage's decisions.md)

1. **RFD4-Proteina on JapanFold, framework fixed, CDR1-3 designed, hotspots from the 1N8Z contacts.** The hosted wrapper keeps kept residues' sequence and one rigid coordinate group, so the framework pose is what we submit
   (superposed on trastuzumab's VH): epitope direction is built in, only one binding mode is sampled (plus a rotated pose). Checkpoint is not disclosed. Cost about $0.003 per design.
2. **Frameworks:** hNbBCII10-FGLA (PDB 3EAK; used by RFantibody and Germinal) and a caplacizumab-like framework (camelid hallmarks). Both with structures.
3. **No RFD4 on the shared VM** (user decision): design stays on JapanFold.
4. **Framework diversification built, tested, cancelled.** First estimate (about 75% identity to the parent framework, so ~10 substitutions) came from an unverified literature summary; the organisers then said fully new CDRs land under 70%;
   measured: the VM abguard gate is a length effect (`gate-analysis.md`). The fallback (`s2_rules/diversify.py`: natural-residue substitutions, protected positions, lite and deep tiers) stays in the repo, not used.
4b. **Generation steered to longer CDR3 and C104-A105-A106 kept fixed:** CDR3 <= 14 residues passed the gate 0-21%, >= 18 residues 80-100% (BCII10); the monoyield paper ties CDR3 packing against that FR3 motif and F42 to monomeric yield.
5. **BoltzGen added, then stopped (user):** its hosted nanobody protocol designs with Cys excluded, so no canonical disulfide; Cys-restored variants were made, then the generator was dropped to focus on RFD4.
6. **Canonical Cys23/Cys104 is a hard gate** (zero-Cys designs lose Tm: NanoMelt median 60.9 C vs 68.2 C restored).
7. **Gate before scoring:** only designs with abguard class `de_novo` and canonical Cys are scored (halves the cost).
8. **Scorers:** OpenFold3 (JapanFold) as primary because it was the only refolder that separated measured trastuzumab binders in earlier team work; OpenDDE-abag, ESMFold-2, second seed and AF2-IG as orthogonal tier for the top third;
   Protenix-v2 (Nebius, no MSA) for the top designs. Gates fixed a priori at the 95th percentile of decoys. Calibration shows weak discrimination on real VHH-antigen pairs: use for ordering only.
9. **Properties:** NanoMelt (CPU; the 8 GB laptop GPU ran out of memory), TNP flags, monoyield nanobody models (VM GPU, bf16, paused unless the user allows it), sequence liabilities, nine CDR2/CDR3 rules (monoyield paper, Fig. 5h; CDR2 rules only for 8-residue CDR2).
10. **Final ranking by Pareto dominance** (bind, epitope, Tm, yield, developability) with weakest-objective tie-break and diversity, next to the weighted composite; missing evidence is pessimistic (0 on the objective, 0.3 prior in the composite).
11. **Developability objective from liabilities (all designs), TNP only sharpens it** (an early version required TNP and was therefore dropped for lack of coverage).
12. **Upstream's refolding oracle is templated and single-sequence** (RFD4-Proteina repo: `oracle_ppi.template_query`); JapanFold/Nebius OpenFold3 cannot template, so our OpenFold3 numbers are not comparable to upstream thresholds. The templated test was later run on the shared VM with the user's permission (decision 13).
13. **Pose check and ranking (16:27): no scorer is dropped.** The templated VM OpenFold3 does not recover the designed pose better than ESMFold-2 or OpenDDE (pose medians 41 / 37 / 44 A; OpenFold3 misfolds the VHH without an MSA), so the ranking keeps the VM OpenFold3 + ESMFold-2 binding axis with OpenDDE as a gate. The binding scores are filters of unknown validity.
14. **Cropped HER2 re-run (539-607):** to limit random placement of the VHH, the VM OpenFold3 and OpenDDE are re-run on the cropped domain for the gate passers, the best-binding eligible designs and the controls (`denovo/vm_templated/make_crop.py`, `run_vm_crop.sh`, `collect_crop.py`, `denovo/s3_score/run_crop.py`). Results are in `results_denovo/crop/`; using them for ranking is a separate, recorded decision.

