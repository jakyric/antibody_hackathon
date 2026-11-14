# krishagen — HER2 / trastuzumab optimisation

Ten optimisation-track candidates for a better trastuzumab, and the tested pipeline that produced them.

This README is organised around the four pitch criteria, so each one can be checked directly:

| Weight | Criterion | Where it is answered |
| ---: | --- | --- |
| **40%** | Judgment — clear reasoning behind the molecules we chose | [§1](#1-judgment-40) |
| **20%** | Execution — a reproducible pipeline in a GitHub repo | [§2](#2-execution-20) |
| **20%** | User fit — what the design is for, and who needs it | [§3](#3-user-fit-20) |
| **20%** | Path to value — a credible route to real use | [§4](#4-path-to-value-20) |

One rule governs every number below: **nothing here is a laboratory measurement of our designs.** Adaptyv
measures affinity, epitope preservation and thermostability weeks after the event. Every claim is labelled
`measured`, `predicted`, `literature` or `unknown`, and we never write "binds".

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

**The novelty rule excludes the zone that works.** `abguard` enforces ≥3 CDR-H3 mutations from trastuzumab
and from every provided sequence. Binder quality falls steeply with distance:

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

**The structural premise of the mutation is wrong.** It is justified on the grounds that position 66 "was
not located at the antigen-binding region". Measured by us in crystal structure 1N8Z
(`results/r66g/R66_GEOMETRY.md`, `scripts/r66_contacts.py`):

| | Distance |
| --- | ---: |
| Arg66 side chain → CDR-L1 Asn30 backbone O | **2.54 Å** (hydrogen bond) |
| Arg66 side chain → CDR-L1 Asp28 | 3.84 Å |
| Arg66 side chain → nearest HER2 atom | 5.45 Å |
| for scale, CDR-L1 → HER2 | 2.99 Å |

Arg66 is not in the interface, but it hydrogen-bonds a loop that is.

**Our own folding agrees in direction** (`results/r66g/R66G.md`). We folded four matched pairs —
trastuzumab and three designs, with and without R66G — with Boltz-2 at multiple seeds, same template:

- **Epitope is preserved.** Pair overlap is 1.00 for trastuzumab, tz-opt-02 and tz-opt-18, against a
  seed-to-seed noise floor of 0.97 measured on the same sequences. Only tz-opt-07 shifts (0.90).
- **Interface confidence drops slightly.** Antibody–antigen ipTM is lower with R66G in **6 of 8**
  comparisons; for trastuzumab itself the drop clears the noise floor at both seeds.

**Why this changes the submission.** The brief counts a property only if it beats the measured baseline —
"ties and regressions are hollow". R66G buys the thermostability vertex and puts the affinity vertex at
risk, so it is a reasonable bet for A-tier and a poor one for S-tier. We therefore submit **matched pairs**:
each stabilised design alongside its unmodified parent, so the lab result separates the loop from the
mutation instead of confounding them.

### 1.4 The candidates

`results/CANDIDATES.md` is the single table: ten rows, the three judged axes, the uncertainty on each, and
what each claim rests on. `results/SUBMISSION.csv` is the submitted file.

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

`make verify` needs the two public Chinery tables in `data/` (`tz_her2_facs_3classes_524346.csv`,
`trastuzumab_her2_spr_all_768.csv`). They are not redistributed here. `make test` passes without them.

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

**The molecule.** Trastuzumab treats HER2-positive breast cancer. We optimise it in the format the challenge
specifies, and the three properties we target — affinity, retained epitope, thermostability — are the ones
that decide whether a binder survives manufacture, shipping and storage, not only a binding assay.

**Who needs the thermostability axis.** A higher Tm is a developability proxy: it is what lets a formulation
tolerate a cold chain that breaks, and it is the usual route to high-concentration subcutaneous dosing. That
is the argument for spending any of a limited candidate budget on it.

**Why matched pairs are the right deliverable here.** A lab that tests ten unrelated molecules gets ten
numbers. A lab that tests matched pairs learns what one mutation *does*, on several backgrounds, with the
loop held constant. Given that the published justification for R66G is contradicted by the crystal
structure, the controlled comparison is the more useful thing to hand over.

**What we explicitly do not claim.** We have no measured affinity, epitope or thermostability for any
design. We do not claim any design binds. We do not rank the ten on binding, because no scorer we tested
can.

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

**Honest next step.** The highest-value experiment is not another model. It is the matched-pair lab result,
which settles a published claim that rests on a structural premise we have shown to be false and that no
independent group has replicated.

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
