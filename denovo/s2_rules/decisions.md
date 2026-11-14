# S2 decisions log (rules + framework diversification)

Written 14 Nov 2026 during the hackathon. Evidence level is stated for every claim: **measured** (we ran it), **read** (in the abguard source/README),
**inferred** (my reasoning, not verified).

## 1. What the VM abguard really checks (read + measured)

Build used: `~/.cache/uv/archive-v0/jsRWGsUwPztt-JwD/bin/abguard` on the team VM, `abguard 0.1.0`, immunum 2.0.0 (the build of 2026-11-11 16:17 is the
newest; the other two builds in the uv cache are byte-identical in the `abguard/` package). `team-validate FILE --all-liabilities -o out.json` returns one row per design with
`errors`, `is_valid`, LAP liabilities and the novelty columns; with the bundled reference (default) the novelty gate is on.

| check | exactly what it compares | evidence |
|---|---|---|
| format / numbering | immunum IMGT numbering, chain label H for a `vhh`, confidence >= 0.1, exactly one domain | read (`config.py`, `schemas.py`) |
| length | each polypeptide <= 250 aa | read |
| unpaired Cys | **even** Cys count over the whole sequence AND even within the numbered domain. It does NOT require the canonical Cys23-Cys104 pair: a sequence with **zero** Cys passes | read; measured: all 8 first BoltzGen designs have 0 Cys and pass (`abguard_pass = True`) |
| track (de novo vs optimisation) | whole numbered Fv Levenshtein identity **to trastuzumab VH only** (`1 - d/max(len)`), < 0.70 -> `de_novo`, >= 0.70 -> `novel` (= trastuzumab variant) | read (`novelty.py`) + measured, our local replica reproduces it to 4 decimals |
| prior-art gate | CDRH3 only: Levenshtein distance of the IMGT CDR3 to **every one of 610,109 reference CDRH3s** shipped in the wheel (sources: 524k FACS-screen variants, 38k patents (USPTO), 38k PLAbDab, 5.9k PDB, 1.1k PLAbDab-nano, 751 SPR designs, 571 proteinbase_public, 108 swissprot, 34 Thera-SAbDab...) and to trastuzumab's own `SRWGGDGFYAMDY`; needs >= 2 mutations if class `de_novo`, >= 3 if class `novel` | read + measured |
| whole-sequence identity to *other* known antibodies | **not done** by this build: the package contains only the CDRH3 index, no sequences | read: reference_data holds one parquet of CDRH3s |

Measured on stand-ins (`results` in `cache/rules_standins.csv`): NbBCII10-FGLA wild type (3EAK) -> rejected, 0 mutations from `pdb:3dwt_A`; trastuzumab VH as a `vhh` -> class `novel`,
rejected (0 mutations from trastuzumab); BCII10 framework + invented CDR1-3 -> **identity to trastuzumab VH 0.711-0.732, class `novel`** (the humanised BCII10 framework sits on the
0.70 line; single residues flip the class), accepted; an extra Cys -> `unpaired_cysteine`.

**The gap.** The task PDF (Figure 2) says the type-aware comparison is VHH against *known VHHs and heavy chains* and "Identity to nearest known antibody under 70 %?" decides the
track. This build only does the trastuzumab half. If the jury's build matches the PDF, a design that passes the VM abguard can still be classed optimisation (or rejected) because of a
known VHH we cannot see. `rules.csv` therefore carries both `denovo_vm_pass` (what we can run) and `denovo_strict_pass` (and < 70 % to our panel). **Open question for the organisers:**
does the jury's abguard compare the whole Fv with a database of known VHHs / heavy chains, and which one?

## 2. The known-antibody panel (measured) and how complete it is

183,567 distinct heavy variable domains (Fv = FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4 as numbered by immunum, deduplicated on that string):
6,363 from every PDB chain with a heavy domain (wwPDB seqres of 14 Nov 2026, 92k unique chains of 100-300 aa screened), 111 reviewed human IGHV germlines (UniProt),
1,349 PLAbDab-nano VHHs, 175,744 PLAbDab heavy chains (GenBank/patent-derived, paired + unpaired). Trastuzumab VH is always included. 6,742 entries are VHH-like (single domain
with camelid hallmarks, or from the nano set).

Honest incompleteness: no OAS (billions of repertoire sequences), no full patent sequence listings (abguard's 38k USPTO CDRH3s prove they exist; we only have PLAbDab's share), no camelid
germline segments (UniProt returned none for llama/alpaca/camel IGHV), no IMGT/GENE-DB, and the HER2 VHHs of the literature that are not in PDB/GenBank (Nikkhoi 2023, Hussack 2018, 5F7 ...):
attempted via Europe PMC full text, none of the papers prints sequences in the text (figures/supplements), so none were added and none were invented. `cache/extra_known.fasta` is read
by `build_panel.py` if a teammate adds sequences with provenance. **Direction of the error: our identity numbers can only be too low** (a missing known VHH can only raise the true max).
That is why the deep tier aims at 0.64-0.68 instead of 0.69.

## 3. What the first real designs look like (measured, 14:12 CET)

`python denovo/s2_rules/rules_summary.py` (regenerated for the final report). RFD4/JapanFold designs on the BCII10 framework (96 designs): all pass the local hard rules and abguard,
but **51/96 are class `de_novo`** (the rest sit at >= 0.70 to trastuzumab VH), and **0/96 are below 0.70 to the nearest known antibody** (median 0.804, range 0.76-0.87; the nearest
entries are BCII10 relatives such as 3EAK, 9vq6_A and PLAbDab nanobodies). The CDR1/CDR2 that RFD4 redesigns lower the identity from 0.90 (CDR3-only change) to 0.80, the framework
residues are the rest. **BoltzGen nanobody-anything designs have no cysteines at all** (C23 -> E, C104 -> A/E/I/W, W118 -> F in the first batch; abguard accepts them because zero is even) ->
`canon_cys_ok = False`, `basic_pass = False`; several also have a CDRH3 within 1-2 edits of a patent CDRH3 (`insufficient_cdrh3_novelty`).

## 4. Framework diversification: method and the choices behind it (NOT RUN in the final pipeline)

> Decision 14:15 CET (user + organisers): the framework stays intact. Everything in this section was built and tested (stand-ins, 3 real RFD4 designs in scratch, 38 variants emitted by the watcher between 14:08 and 14:15, now parked in `results/denovo/parked_fwdiv.csv`) but its variants were never refolded or profiled. Kept as a documented fallback.

Goal: reach < 70 % to anything known without breaking the VHH. Alternatives I considered:

* **Hosted ProteinMPNN with framework designable / CDRs fixed (what the brief suggested).** The hosted adapter (`infer_proteinmpnn_native`, v_48_020) has no fixed-position or tied-position
  input, only chains. I ran it on the 3EAK backbone (96 samples, T = 0.3, Nebius op-a8d5dcf3562d5b00): mean sequence recovery 0.56, sequences look unlike VHHs (`GITLTESGGGTVAPGG...`;
  it has no antibody prior). Decision: use it only as a *vote* (weight 0.004 per unit of vote frequency) on the framework positions, never as the source of substitutions. **rejected as generator.**
* **AntiFold** keeps the framework and redesigns CDRs: wrong direction. **rejected.**
* **Panel-conditioned natural substitutions with greedy identity descent (chosen).** Every allowed substitution is a residue that >= 1.5 % of known heavy domains carry at that IMGT
  position (--min-freq; 3 % left only 38/77 positions with any alternative, 1.5 % gives 47/77 and 106 substitutions, so the deep tier needs it). Each step picks the substitution that most
  lowers the mean of the top-10 nearest-known identities (+ top-1), with small bonuses for naturalness (log-odds), the MPNN vote and residues of a different human germline. Protected and
  buried positions as listed in README. Variants of one parent differ by random position masks and sampled choices.
* **ESM-2 650M masked-marginal check** (Nebius, `score-variants`): run on the stand-in variants (op-b8733dd4b66a2a19): summed log-odds -25 to -72, **mean -1.1 to -2.9 per substitution**
  (n = 18 variants, 17-27 substitutions each). There is no calibration (no natural-homolog baseline was run), so it is *informational*: it ranks variants of one parent (fewer/more
  conservative substitutions score better) but I set no pass/fail threshold on it. `esm_filter.py` is ready for the real variants.

**Two tiers** (coordinator request): `lite` = the minimum natural substitutions (<= 8) that bring trastuzumab-VH identity to <= 0.66; `deep` = the descent above. Measured on real
designs (`before_after.md`, 3 parents x 2 lite + 3 deep): lite needs 4-8 substitutions (0.685-0.724 -> 0.653-0.659 vs trastuzumab) and leaves the nearest-known identity at 0.75-0.78;
deep needs 16-25 substitutions (id. to parent 0.80-0.86) to reach 0.63-0.68. Roughly 1.2 substitutions per percentage point. **Cost is unknown until S3/S4 report**: OpenFold3 retention
(`fwdiv_eval.py`) and NanoMelt Tm per tier are the evidence. The deep tier is an unverified hedge against an unseen database.

**Cys repair (`-cr`)**: BoltzGen designs lack the canonical disulfide. The repair restores C23/C104 (and W118) and nothing else, one variant per design, so S3/S4 can say whether the
Cys-free design was merely mis-written. It is a rescue, not a claim that the BoltzGen design binds.

## 5. Before / after table

Snapshot, real RFD4 designs `bcii10-rfd4jf-s1100-000/-001/-002` (scratch run, `cache/dv_real3_rep.csv`); the live table for everything the watcher emits is `before_after.md`
(`python denovo/s2_rules/make_table.py denovo/s2_rules/fwdiv_report.csv > denovo/s2_rules/before_after.md`; the watcher's `fwdiv_report.csv` holds the Cys-repair rows' companions only if it ran long enough, see `parked_fwdiv.csv`).

| variant | tier | subs | id. to parent | vs trastuzumab VH before -> after | vs nearest known before -> after (nearest id) | <=0.66 trast | <0.70 known |
|---|---|---|---|---|---|---|---|
| bcii10-rfd4jf-s1100-000-fl1 | lite | 8 | 0.935 | 0.718 -> 0.653 | 0.815 -> 0.75 (9vq6_A) | yes | no |
| bcii10-rfd4jf-s1100-000-fl2 | lite | 8 | 0.935 | 0.718 -> 0.653 | 0.815 -> 0.75 (9vq6_A) | yes | no |
| bcii10-rfd4jf-s1100-000-fd1 | deep | 17 | 0.863 | 0.718 -> 0.605 | 0.815 -> 0.677 (9vq6_A) | yes | yes |
| bcii10-rfd4jf-s1100-000-fd2 | deep | 20 | 0.839 | 0.718 -> 0.589 | 0.815 -> 0.653 (9vq6_A) | yes | yes |
| bcii10-rfd4jf-s1100-000-fd3 | deep | 22 | 0.823 | 0.718 -> 0.573 | 0.815 -> 0.637 (9vq6_A) | yes | yes |
| bcii10-rfd4jf-s1100-001-fl1 | lite | 8 | 0.935 | 0.724 -> 0.659 | 0.837 -> 0.78 (pl:AOG97027.1) | yes | no |
| bcii10-rfd4jf-s1100-001-fl2 | lite | 8 | 0.935 | 0.724 -> 0.659 | 0.837 -> 0.772 (pl:QBE07652.1) | yes | no |
| bcii10-rfd4jf-s1100-001-fd1 | deep | 21 | 0.829 | 0.724 -> 0.577 | 0.837 -> 0.675 (pl:AML32071.1) | yes | yes |
| bcii10-rfd4jf-s1100-001-fd2 | deep | 23 | 0.813 | 0.724 -> 0.553 | 0.837 -> 0.659 (pl:AOG97027.1) | yes | yes |
| bcii10-rfd4jf-s1100-001-fd3 | deep | 25 | 0.797 | 0.724 -> 0.545 | 0.837 -> 0.634 (pl:QBE07652.1) | yes | yes |
| bcii10-rfd4jf-s1100-002-fl1 | lite | 4 | 0.968 | 0.685 -> 0.653 | 0.805 -> 0.773 (3eak_A) | yes | no |
| bcii10-rfd4jf-s1100-002-fl2 | lite | 4 | 0.968 | 0.685 -> 0.653 | 0.805 -> 0.773 (3eak_A) | yes | no |
| bcii10-rfd4jf-s1100-002-fd1 | deep | 16 | 0.871 | 0.685 -> 0.589 | 0.805 -> 0.68 (3eak_A) | yes | yes |
| bcii10-rfd4jf-s1100-002-fd2 | deep | 19 | 0.847 | 0.685 -> 0.556 | 0.805 -> 0.656 (3eak_A) | yes | yes |
| bcii10-rfd4jf-s1100-002-fd3 | deep | 22 | 0.823 | 0.685 -> 0.548 | 0.805 -> 0.633 (3eak_A) | yes | yes |

Stand-in parents (BCII10 and caplacizumab-like frameworks with invented CDR1-3, used before real designs existed):

| variant | tier | subs | id. to parent | vs trastuzumab VH before -> after | vs nearest known before -> after (nearest id) | <=0.66 trast | <0.70 known |
|---|---|---|---|---|---|---|---|
| si3-bcii10-a-fl1 | lite | 8 | 0.937 | 0.717 -> 0.654 | 0.85 -> 0.795 (pl:QEP15795.1) | yes | no |
| si3-bcii10-a-fl2 | lite | 8 | 0.937 | 0.717 -> 0.654 | 0.85 -> 0.795 (pl:QEP15795.1) | yes | no |
| si3-bcii10-a-fd1 | deep | 22 | 0.827 | 0.717 -> 0.575 | 0.85 -> 0.677 (pl:ACK08510.1) | yes | yes |
| si3-bcii10-a-fd2 | deep | 25 | 0.803 | 0.717 -> 0.528 | 0.85 -> 0.654 (pl:ACK08510.1) | yes | yes |
| si3-bcii10-a-fd3 | deep | 27 | 0.787 | 0.717 -> 0.528 | 0.85 -> 0.638 (pl:ACK08510.1) | yes | yes |
| si3-bcii10-b-fl1 | lite | 8 | 0.937 | 0.732 -> 0.669 | 0.835 -> 0.787 (pl:QKY76618.1) | no | no |
| si3-bcii10-b-fl2 | lite | 8 | 0.937 | 0.732 -> 0.669 | 0.835 -> 0.787 (pl:QKY76618.1) | no | no |
| si3-bcii10-b-fd1 | deep | 20 | 0.843 | 0.732 -> 0.583 | 0.835 -> 0.679 (9g48_C) | yes | yes |
| si3-bcii10-b-fd2 | deep | 24 | 0.811 | 0.732 -> 0.567 | 0.835 -> 0.654 (pl:AOC15334.1) | yes | yes |
| si3-bcii10-b-fd3 | deep | 25 | 0.803 | 0.732 -> 0.551 | 0.835 -> 0.638 (9ena_B) | yes | yes |
| si3-capla-like-fl1 | lite | 3 | 0.976 | 0.683 -> 0.659 | 0.856 -> 0.856 (pl:CAQ53179.1) | yes | no |
| si3-capla-like-fl2 | lite | 3 | 0.976 | 0.683 -> 0.659 | 0.856 -> 0.848 (pl:CAQ53179.1) | yes | no |
| si3-capla-like-fd1 | deep | 23 | 0.813 | 0.683 -> 0.537 | 0.856 -> 0.672 (pl:CAQ53179.1) | yes | yes |
| si3-capla-like-fd2 | deep | 26 | 0.789 | 0.683 -> 0.496 | 0.856 -> 0.656 (pl:CAQ53179.1) | yes | yes |
| si3-capla-like-fd3 | deep | 28 | 0.772 | 0.683 -> 0.504 | 0.856 -> 0.672 (pl:CAQ53179.1) | yes | yes |

## 5b. Gate behaviour with the framework intact

See [`gate_drivers.md`](gate_drivers.md): identity to trastuzumab VH = 1 - Levenshtein/max(Fv length, 120); for BCII10-framework designs it is ~87.2/L, so the class is a length effect (line at Fv length ~125) plus about +-0.02 from chance CDR matches; the caplacizumab-like framework is 98 % `de_novo` at any length.

## 6. Risks I could not remove

1. The jury's identity rule is unconfirmed (section 1). If it is the PDF's, the lite tier is not enough and the deep tier is verified only against our panel.
2. 20-28 framework substitutions can cost expression, Tm or binding; nothing in S2 measures that (S3 and S4 do).
3. immunum is the only authority on CDR boundaries here (same as abguard), including its weak numbering of some designs (`imm_conf` < 0.5).
4. BoltzGen's missing Cys: if the generator's own structure model was fine without the disulfide, our `-cr` repair may change the fold; check S3's pose RMSD of `-cr` vs parent.
5. `nearest_known` is a Levenshtein identity on the numbered Fv (abguard's formula); other jury definitions (identity on aligned positions only, CDR-excluded) give different numbers.
