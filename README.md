# Evidence-first antibody design

We built a reproducible antibody-design pipeline that prioritises experimental evidence over model
confidence.

Our central finding was that **highly accurate computational scores did not reliably predict experimental
binding affinity**. A model we trained reached AUROC 0.983 on its own label and still lost to plain edit
distance against measured KD. Rather than optimise a misleading metric, we built a selection process on
things we could verify: novelty against the full provided dataset, chemical liabilities from the
organisers' own reference, structural consistency across folding models and seeds, and explicit
uncertainty on every claim.

We deliver **ten trastuzumab-derived scFv/paired candidates** for experimental evaluation, including a
controlled investigation of the VL R66G thermostability mutation.

Our goal is not to claim that AI has already designed a superior antibody. It is to help scientists make
better decisions about which antibodies to test.

```bash
make setup && make test     # 97 tests, no data needed
make demo                   # a complete run: target in -> ranked candidates out
```

**Nothing in this repository is a laboratory measurement of our designs.** Adaptyv measures affinity,
epitope preservation and thermostability weeks after the event. Every claim is labelled `measured`,
`predicted`, `literature` or `unknown`, and we never write "binds".

---

## What changed because we tested our assumptions

The pipeline's main output is not the candidate list. It is a set of conclusions we did not want, each of
which changed what we submitted.

| We assumed | We tested it | What we did |
| --- | --- | --- |
| Our binder model is good — AUROC 0.983 | Against measured SPR it scores +0.16, worse than edit distance's +0.23 | **Removed binding from the ranking entirely** |
| StaB-ddG can order our designs | +0.38 overall, but **+0.14 [−0.26, +0.53]** at the 4–5 mutations our designs sit at | **Did not run it on the candidates** |
| R66G is a free thermostability win | Its own paper reports a 2.4× affinity loss; an independent group installed that residue *to gain* affinity | **Submitted matched pairs instead of stabilising everything** |
| Our novelty check covers the dataset | It assumed all 525,063 CDR-H3s are 10-mers. **21 are not**, and those were silently skipped | **Rewrote it to measure edit distance exhaustively** |
| Our cysteine filter works | It counted across the whole molecule; the rule is per domain | **Fixed; an odd count in each chain is now caught** |
| The repo builds for a judge | `make setup` **failed** on a clean Debian box — venv without pip | **Fixed; verified from a fresh clone** |
| Two folding seeds are enough | A third seed flipped the sign of our R66G interface result | **Report it as a hint, not a finding** |

Each row is a filter that returned a clean answer to the wrong question. That is the failure mode this
pipeline is built to catch.

---

This README is organised around the four pitch criteria, so each can be checked directly:

| Weight | Criterion | Where it is answered |
| ---: | --- | --- |
| **40%** | Judgment — clear reasoning behind the molecules we chose | [§1](#1-judgment-40) |
| **20%** | Execution — a reproducible pipeline in a GitHub repo | [§2](#2-execution-20) |
| **20%** | User fit — who actually benefits | [§3](#3-user-fit-20) |
| **20%** | Path to value — a credible route to real use | [§4](#4-path-to-value-20) |

This is one person's branch of the team's private repository. It holds methods, tests, measured results and
the candidates. It holds no credentials and does not redistribute the organisers' task text or data files.

---

## 1. Judgment (40%)

### 1.1 The honest starting point: what we cannot do

We spent the day establishing what our own methods *cannot* support, because a ranking built on a scorer
that does not work is worse than no ranking.

**Beating trastuzumab is hard.** Of 182 provided variants with a measured KD, **6 (3.3%)** beat
trastuzumab's median of 0.135 nM. The median variant is ~590× worse. These are published, professionally
designed variants from a study doing exactly this task.

**The novelty rule excludes the zone that works.** Two separate constraints apply, and they are not the
same number:

| Constraint | Minimum mutations | Where it comes from |
| --- | ---: | --- |
| CDR-H3 identity to trastuzumab under 70% | **4** | 3 mutations leaves 70% on the 10-mer (77% on abguard's 13-mer), which is not *under* 70% |
| Distance from every provided sequence | **3** | `abguard`'s `required_cdrh3_mutations: 3` |

Our designs carry **4 mutations from trastuzumab** and sit at **3 from the nearest provided sequence**, so
they satisfy both. Binder quality falls steeply with distance:

| CDR-H3 mutations | n | fraction "high" |
| --- | --- | --- |
| 2 | 600 | 97.5% |
| 3 | 1,490 | 90.0% |
| **4** | 6,392 | **74.6%** |
| 7 (median) | 160,128 | 33.6% |

All our designs therefore sit at exactly 4 mutations — as close to the original as the rules permit.

**The FACS label does not transfer to affinity.** We trained a net to AUROC 0.983 on FACS class, then tested
it against SPR: Spearman **+0.16 [−0.03, +0.33]**, worse than edit distance's **+0.23**. On the 27 sequences
measured both ways, "medium" binds *tighter* than "high", and the FACS class predicts affinity at AUROC
0.61. The net saturated a label that is itself barely informative about affinity.

**Consequence, and this is the central judgment call: we do not rank these ten on binding.** All ten sit at
4 mutations, so edit distance — the only scorer that transferred — is constant across them and orders
nothing. We say so in the table rather than inventing a column.

### 1.2 StaB-ddG: tested, and rejected for ranking

Before spending compute, we validated StaB-ddG on 60 provided variants whose KD is already measured
(`results/stabddg/ANALYSIS.md`):

| Comparison | Spearman vs log₁₀KD |
| --- | ---: |
| StaB-ddG, all 60 | +0.38 [+0.14, +0.59] |
| mutation count alone | **+0.54** |
| StaB-ddG at 4–5 mutations (n=25) | **+0.14 [−0.26, +0.53]** |

It does not beat counting mutations, and **at the mutation count our designs actually sit at, its interval
includes zero**. It is also a *binding* ddG, not a folding free energy, and every variant in the run mutates
the heavy chain only — so it says nothing about thermostability and nothing about our light-chain mutation.
We did not run it on the candidates.

### 1.3 The one real trade: VL R66G

Six of the ten rows carry a single light-chain mutation, VL R66G. It is the only thermostability claim in
the submission, and we interrogated it rather than inheriting it.

**What the source paper reports** (Okazaki et al. 2023, *ACS Omega* 8(25):22945, trastuzumab scFv):

| | wild type | VL-R66G |
| --- | --- | --- |
| Tm (DSF) | 61.8 ± 0.2 °C | **66.7 ± 0.2 °C** (+4.9) |
| KD (SPR) | 2.0 ± 0.4 nM | **4.8 ± 0.7 nM** |

The paper describes the affinity as "unchanged". Its own table shows a **2.4× loss with non-overlapping
error bars**, no n and no statistical test. We do not repeat the word "unchanged".

**Independent evidence points the same way.** Carter et al. 1992 (*PNAS* 89:4285) engineered this residue in
the opposite direction: huMAb4D5-3 (VL66 **Gly**) KD 4.4 nM vs huMAb4D5-5 (VL66 **Arg**) KD 1.1 nM, differing
only at VL66. Trastuzumab has Arg66 on purpose.

**Our structural analysis challenges the assumption that R66G is functionally neutral.** The mutation is
justified on the grounds that position 66 "was not located at the antigen-binding region". Arg66 indeed
makes no direct contact with HER2 — but it is positioned to stabilise a loop that does. Measured by us in
crystal structure 1N8Z
(`results/r66g/R66_GEOMETRY.md`, `scripts/r66_contacts.py`):

| | Distance |
| --- | ---: |
| Arg66 side chain → CDR-L1 Asn30 backbone O | **2.54 Å** (consistent with a hydrogen bond) |
| Arg66 side chain → CDR-L1 Asp28 | 3.84 Å |
| Arg66 side chain → nearest HER2 atom | 5.45 Å |
| for scale, CDR-L1 → HER2 | 2.99 Å |

Arg66 is not in the interface. It sits against CDR-L1 at a distance consistent with a hydrogen
bond, and CDR-L1 contacts HER2 at 2.99 Å. Distance alone does not prove the bond, and this does
not show the original paper's geometry to be wrong — but it does give a **plausible mechanism**
for the affinity loss those authors measured and then described as absent.

**Our own folding agrees in direction** (`results/r66g/R66G.md`). We folded four matched pairs —
trastuzumab and three designs, with and without R66G — with Boltz-2 at three seeds (42, 77, 123), same template:

- **Epitope is preserved.** Median pair overlap is 1.00 for trastuzumab, tz-opt-02 and tz-opt-18, against
  a seed-to-seed noise floor of 0.97 measured on the same sequences. tz-opt-02 holds 1.00 at all three
  seeds; tz-opt-18 holds it at two of three. **tz-opt-07 is the exception** — 0.86 / 0.94 / 0.90, below
  the noise floor at every seed. R66G moves that design's predicted epitope and does not move the others',
  which is a per-design result, not a property of the mutation.
- **Interface confidence leans down, but only leans.** Antibody–antigen ipTM is lower with R66G in
  **8 of 12** comparisons, median **−0.012**. Per design, almost all of that sits inside the seed-to-seed
  spread of a single sequence (median 0.014). This is a directional hint consistent with the published
  affinity loss — it is **not** independent confirmation, and we do not present it as one.

  The third seed is why that sentence is hedged. On two seeds the drop for trastuzumab looked like it
  cleared the noise floor; the third seed flipped its sign. Two seeds would have let us overstate this.

**Why this changes the submission.** The brief counts a property only if it beats the measured baseline —
"ties and regressions are hollow". R66G *may* improve thermostability — no design of ours has had its
Tm measured — while putting the affinity vertex at
risk, so it is a reasonable bet for A-tier and a poor one for S-tier. We therefore spend part of the
submission on **matched pairs** rather than stabilising everything.

**Exactly what the ten rows are**, since the pair argument depends on it:

| | Count | Which |
| --- | ---: | --- |
| Stabilised (carry R66G) | 6 | `tz-opt-02/07/18/01/15/03-stab` |
| Unmodified parents | 4 | `tz-opt-02/07/18/01` |
| **Complete matched pairs** | **4** | `tz-opt-02`, `tz-opt-07`, `tz-opt-18`, `tz-opt-01` |
| Stabilised with no parent in the list | 2 | `tz-opt-15-stab`, `tz-opt-03-stab` |

Ten unique rows allow at most five complete pairs; we submit four, and spend the two remaining slots on
additional loops rather than a fifth control. **The two unpaired stabilised designs cannot support a
controlled read of R66G** — if only those are tested, we learn nothing about the mutation. That is the
cost of the choice, and it is why the pairs are ranked above them.

### 1.4 The candidates

`results/CANDIDATES.md` is the single table: ten rows, the three judged axes, the uncertainty on each, and
what each claim rests on. `results/SUBMISSION.csv` is the submitted file.

**One row, in full, so the format is clear** — `tz-opt-02-stab`, our first-ranked candidate:

| | |
| --- | --- |
| **Why it is in the list** | Best developability profile of all twenty: it drops trastuzumab's DG isomerisation motif and the CDR-H3 tryptophan, giving 7 liability motifs against trastuzumab's 8, with one fewer at high severity. It is the only design strong on both axes we could actually measure. |
| **CDR-H3** | `MGPEPFYAMD` — 4 mutations from trastuzumab, 3 from the nearest of 525,063 provided sequences (`CGAEGFYAMD`) |
| **Affinity** | `unknown`. No scorer we have orders designs at 4 mutations. R66G additionally costs affinity in the literature (2.0 → 4.8 nM) |
| **Epitope** | `predicted` 0.89 vs trastuzumab; R66G pair overlap **1.00 at all three seeds**, noise floor 0.97 |
| **Thermostability** | `literature` +4.9 °C, on a different molecule, one lab, never replicated |
| **What the lab should find out** | Whether R66G's Tm gain transfers to a redesigned loop, and what it costs in SPR. Its parent `tz-opt-02` is in the list for exactly this comparison. |
| **What could go wrong** | The affinity loss exceeds the Tm gain in value; the Tm gain does not transfer out of Okazaki's construct; or the design simply fails to express |

The other nine rows are in the same form in `results/CANDIDATES.md`.

---

## 2. Execution (20%)

### 2.1 Reproduce it

```bash
make setup     # venv + numpy/pytest (bootstraps pip where ensurepip is missing)
make test      # 97 tests
make verify    # the benchmark: does our scorer recover answers we already know?
```

`make test` and `make verify` were run **from a clean clone of this branch, on a machine that had never seen
the project** — 97 passed, 6 skipped, benchmark exit 0, reproducing the headline
`edit_distance Spearman +0.226` in the legal zone. That run found a real bug: on a stock Debian/Ubuntu
python, `python3 -m venv` creates an environment *without* pip, and every later step died with "No module
named pip". `make setup` now bootstraps pip itself, so the repository builds on a machine with no
`python3-venv` package.

### 2.1b A complete run: target in, ranked candidates out

`make demo` is the whole method in one command, on a fresh target definition:

```bash
make demo        # ~1 minute
```

| Step | What happens | Typical numbers |
| ---: | --- | --- |
| 1 | Propose CDR-H3 variants from a position-weight matrix fitted to the FACS library | 5,000 unique |
| 2 | **Entry gate** — novelty against all 525,063 provided sequences | 3,354 pass (67%); 1,646 too close |
| 3 | Drop what the binder model rates unlikely (top 30% kept) | 1,006 left |
| 4 | Rank on the only signal that transferred to measured SPR | 907 clear the full gate |
| 5 | Select a diverse set, no near-duplicates of each other | 10 chosen |
| out | `results/demo/designs_demo.csv` — ranked, with a written rationale per design | |

It writes to `results/demo/`, so it never touches the submitted designs. To point it at a different
target, replace the reference sequences in `src/abx/sequences.py` and the assay tables in `data/`; the
gate, the ranking and the reporting are target-agnostic.

Then check any candidate table against the entry rules:

```bash
python scripts/requirement_check.py --submission results/SUBMISSION.csv
```

### 2.1c Getting the data

`make demo` and `make verify` need two public tables in `data/`. `make test` does **not** — it passes on a
clean clone with no data at all.

| File | What it is |
| --- | --- |
| `data/tz_her2_facs_3classes_524346.csv` | 524,346 CDR-H3 variants with a FACS binding class |
| `data/trastuzumab_her2_spr_all_768.csv` | 768 variants, 190 with a measured KD |

Both come from Chinery et al., bioRxiv [10.1101/2024.03.26.586756](https://doi.org/10.1101/2024.03.26.586756),
distributed by the organisers at `s3://apheris-oslo-hackathon/data/`. We do not redistribute them here. If
a file is missing, every script that needs it says which file and where to put it rather than failing
obscurely.

### 2.2 One command that says whether we pass

```bash
python scripts/requirement_check.py --md results/REQUIREMENTS.md
```

One row per design, one column per entry requirement, and a verdict that is **PASS**, **FAIL** or
**UNKNOWN** — never a blank, and never a PASS for something that cannot be tested here. Current state
(`results/REQUIREMENTS.md`): **all ten designs PASS every locally checkable requirement.**

| Requirement | Status |
| --- | --- |
| Format is VHH / scFv / Fab | PASS, 10/10 |
| Total length ≤ 250 aa | PASS, 10/10 (227 aa; 23 to spare) |
| No unpaired cysteines, **per domain** | PASS, 10/10 |
| CDR-H3 identity to trastuzumab under 70% | PASS, 10/10 (60% on the 10-mer; 4 mutations) |
| ≥3 mutations from every provided sequence | PASS, 10/10 (distance 3 against all 525,063) |
| `immunum` numbering confidence | **UNKNOWN** — only `abguard` can answer |
| Identity to known antibodies (de novo rule) | **UNKNOWN** — needs Adaptyv's database |

**UNKNOWN is never reported as a pass.** The two UNKNOWN rows cannot be satisfied on any machine without
the organisers' tooling, so UNKNOWN warns and FAIL blocks; the command exits 1 only on a FAIL.
`--strict-unknown` makes UNKNOWN blocking too, for anyone who wants the stricter behaviour.

**The margin is printed as a number, because we have none.** Every design sits at exactly 3 mutations from
the nearest provided sequence and the requirement is 3 — a margin of zero, by construction, because the
generator was built to stay as close to trastuzumab as the rules permit. A design that only just passes can
fail if the organisers' tool counts even slightly differently, so this is shown rather than hidden behind a
tick.

**Where the sources disagree we take the strictest reading.** The slides say "≤70%" and "at least 2
mutations"; the brief says "under 70%" and "strictly more than 2"; `abguard` enforces 3. We apply under-70%
and ≥3, which satisfies every reading.

**The check is tested against the thing it must catch.** `tests/test_requirement_check.py` asserts that
trastuzumab itself FAILS both novelty rules, that a design one mutation from a provided sequence FAILS,
and that an odd cysteine count in *each* chain is caught even though the molecule-wide total is even.

### 2.3 The entry gate, checked against the whole dataset

Every submitted row passes the organisers' own `abguard`, run on the exact committed file
(`results/gate/final_validation.txt`, which records the md5 of the file it validated):

```
1. the ten alone               -> EXIT 0   (10/10 valid)
2. the ten plus trastuzumab    -> EXIT 1   (trastuzumab correctly rejected)
```

**Trastuzumab is always included as a negative control.** A gate that passes everything is
indistinguishable from a gate that is switched off.

The novelty rule is measured against *all* 525,063 provided CDR-H3s, not a subset
(`scripts/verify_novelty_full.py`, `results/gate/novelty_full.csv`). This check exists because our abguard
reference was built from a 4,579-sequence shortlist, on the argument that anything more than 4 mutations
away cannot change a verdict. The argument is sound, but it assumed every provided sequence is a 10-mer.
They are not: **20 are 9 residues and one is 16**, and a Hamming comparison silently scored those as
unreachable. Measured exhaustively — Hamming where lengths match, Levenshtein where they do not — every
design sits at distance **3** (the requirement) and trastuzumab sits at **0**.

### 2.4 Every number is traceable

`tools/record_run.py` wraps each run and records the command, the commit, the seed and input checksums into
a run folder and a ledger. `tools/scan_repo.py` scans for credentials and oversized files before anything is
published. Sequence columns are stripped from shared result tables by `scripts/share_results.py`, and a test
asserts that none survive.

### 2.5 What we got wrong, and how we caught it

- A disulfide filter required an atom that **glycine does not have**, so glycine positions were silently
  skipped.
- A cysteine filter counted cysteines across the **whole molecule** when the rule applies per domain.
- `novelty-report` does **not** reflect `--reference-db`; only `team-validate` enforces it. Reading a
  novelty report as evidence about prior art nearly produced a false all-clear.

The first two returned clean answers to the wrong question. That is why the gate reports `unknown` for
anything it cannot verify locally, never a pass, and why every gate run includes a control that must fail.

---

## 3. User fit (20%)

**Our users are antibody engineers and experimental scientists who must decide which candidates deserve
expensive laboratory testing.** That is the decision this pipeline is built around.

Existing AI tools can produce thousands of candidates and highly confident predictions, but those scores
do not always predict experimental success. We demonstrated that on our own work: a model scoring AUROC
0.983 on its own label was beaten by edit distance against measured affinity, and a published ddG model
lost its signal entirely in the regime our designs occupy.

Our pipeline helps researchers avoid false confidence. It validates scorers against experimental evidence
*in the regime where they will be used*, filters invalid and chemically risky designs, reports a noise
floor beside every prediction, and states what it cannot know. The result is not just a ranked list — it
is an **auditable experimental plan**: every candidate carries what we expect, how uncertain it is, what
the claim rests on, and what the assay would have to show to refute it.

**Why matched pairs are the right deliverable for this user.** A lab that tests ten unrelated molecules
gets ten numbers. A lab that tests matched pairs learns what one mutation *does*, with the loop held
constant. We hand over a controlled experiment, not a wish list.

**What we explicitly do not claim.**

- We have **no measured affinity, epitope or thermostability** for any design, and we do not claim any
  design binds.
- We do **not** rank the ten on binding, because no scorer we tested can order designs at 4 mutations.
- We deliver an **scFv / paired variable fragment**, not the clinical IgG. Properties do not transfer
  between formats automatically, and the R66G literature is scFv-only.
- A higher Tm is a **developability proxy, not a formulation result**. It does not by itself show that a
  formulation survives a broken cold chain or supports high-concentration subcutaneous dosing; those need
  separate aggregation, solubility and viscosity work that we have not done.
- We do **not** claim the method saves time or money. We have not compared it against how candidates are
  normally selected, and until we do, that claim would be exactly the kind of unvalidated confidence this
  pipeline exists to catch.

---

## 4. Path to value (20%)

**For the molecules.** The submission is built from matched pairs, so the assay is interpretable whichever
way it falls: if the stabilised members gain Tm without losing SPR, R66G transfers to a redesigned loop and
is worth carrying forward; if they lose affinity as Okazaki's and Carter's data predict, the parents are the
candidates and the mutation is dropped. An all-stabilised submission would not distinguish these.

**For the method, which is the more transferable product.** Three pieces generalise beyond HER2:

1. **Validate before you rank.** Every scorer is tested against data where the answer is already measured,
   *in the regime where it will be used*. StaB-ddG passed the easy test (+0.38 overall) and failed the one
   that mattered (+0.14 at fixed mutation count). Checking only the first would have put a meaningless
   column in the submission.
2. **Noise floors before effects.** Every predicted quantity is reported against the spread of the same
   model on the same sequence at another seed. The epitope number moves by up to 0.105 between seeds — the
   size of the entire spread between our designs — so it cannot rank them, and we say so.
3. **Controls that must fail.** Trastuzumab in every gate run; planted near-duplicates in the reference
   check. A filter is only evidence once you have watched it reject something.

### The route from here

| Phase | What we do | The value |
| --- | --- | --- |
| **Now** | Generate candidates and quality-assure them | Invalid and weakly-justified designs never reach the bench |
| **The lab** | Adaptyv measures binding, epitope and stability | We find out what actually works |
| **Next version** | Use those measurements to score our own selection process | Fewer mis-prioritised experiments, measured rather than asserted |
| **Beyond** | Run the method on other antibodies and targets | Tests whether the approach generalises |

**The key point: the method improves even if none of our ten beats trastuzumab.** We have handed the lab a
falsifiable hypothesis about R66G. If it raises Tm across different CDR-H3 backgrounds, we have found
something portable. If it costs affinity as Okazaki's and Carter's data predict, we have identified a real
limitation and the parents become the candidates. Both outcomes are worth the assay; an all-stabilised
submission would have produced neither.

**The measurable decision point.** The next version of this pipeline is only better if, scored against the
lab results, its ranking beats our current one — which is, by our own admission, no ranking at all on the
binding axis. That is a low bar and an honest one, and it is the first number we will report.

---

## Layout

| Path | What it is |
| --- | --- |
| `results/REQUIREMENTS.md` | **Pass / fail / unknown per design, per entry requirement, with margin** |
| `results/CANDIDATES.md` | **The ten candidates, three judged axes, uncertainty and basis** |
| `results/SUBMISSION.csv`, `SUBMISSION.md` | The submitted file, and why each row is there |
| `results/r66g/` | The R66G investigation: crystal geometry, matched-pair folding, seed noise |
| `results/stabddg/ANALYSIS.md` | Why StaB-ddG was not used for ranking |
| `results/gate/` | abguard runs, the exhaustive novelty check, the controls |
| `results/ranked_submission.csv`, `RANKING.md` | All twenty designs with sequences |
| `results/EVIDENCE.md` | Per candidate, which claims are measured and which expected |
| `FINDINGS.md` | The full record of what we measured, including what did not work |
| `src/abx/data.py` | The FACS and SPR tables, with the CDR-H3 conventions made explicit |
| `src/abx/gate.py` | Entry rules, run locally; reports `unknown` rather than a pass |
| `src/abx/benchmark.py` | Does a scorer rank known answers? All rows, and the legal region |
| `src/abx/splits.py` | Splits that do not leak through near-neighbours |
| `src/abx/model.py`, `rank.py`, `pipeline.py` | The binder model and the ranking |
| `src/abx/metrics.py` | Spearman, AUROC, top-k precision, bootstrap intervals, no scipy |
| `src/abx/thermo.py` | Disulfide geometry, building the C-beta glycine does not have |
| `scripts/requirement_check.py` | The one command: every entry requirement, per design, with margin |
| `scripts/verify_novelty_full.py` | Exhaustive novelty against all 525,063 provided CDR-H3s |
| `scripts/r66_contacts.py` | Arg66 geometry measured in 1N8Z |
| `scripts/r66g_analyse.py` | Matched-pair epitope and interface analysis against seed noise |
| `scripts/stabddg_analyse.py` | The StaB-ddG validation |
| `scripts/epitope.py`, `fold_analyse.py`, `fold_queries.py` | Folding panel, epitope and confidence analysis |
| `tools/record_run.py`, `tools/scan_repo.py` | Run provenance; secret and size scanning |
| `tests/` | Tests for all of the above |

Data source for the HER2 tables: Chinery et al., bioRxiv
[10.1101/2024.03.26.586756](https://doi.org/10.1101/2024.03.26.586756).
