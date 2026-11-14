# What we measured, and what it means for the submission

Everything below is computed from the provided data. Every claim can be
re-run: `python scripts/status.py`

Norwegian version: [FUNN.md](FUNN.md)

---

## 1. Beating trastuzumab is brutally hard

Of the 182 variants with a measured KD in the SPR table:

| | |
| --- | --- |
| beat trastuzumab's median (0.135 nM) | **6 (3.3 %)** |
| beat trastuzumab's best measurement (0.029 nM) | **1 (0.5 %)** |
| median variant | 79.5 nM — about **590× worse** |

These are published, professionally designed variants from a study doing
exactly this task. We should not promise the jury that we win the
affinity axis. We should show that we know how hard it is.

## 2. The novelty filter excludes precisely the zone that works

Under 70 % identity on a 10-mer means **at least 4 mutations**
(7 of 10 matches is 70 %, which is not *under*; 6 matches is 60 %).

But binder quality falls steeply with distance:

| CDRH3 mutations | n | fraction "high" |
| --- | --- | --- |
| 1 | 50 | 96.0 % |
| 2 | 600 | 97.5 % |
| 3 | 1,490 | 90.0 % |
| **4** ← first legal | 6,392 | **74.6 %** |
| 5 | 28,510 | 59.4 % |
| 7 (median in set) | 160,128 | 33.6 % |
| 10 | 19,739 | 15.8 % |

We are not allowed to stay close, and going far is expensive. All 20 of
our designs therefore sit at **exactly 4 mutations** — as close to the
original as the rules permit.

## 3. The label did not transfer to affinity

This is the most important thing we found today, and it changed how we
rank.

We trained a small net on the 524,346 FACS-measured variants. On a split
where no test sequence lies closer than 3 mutations to anything in
training:

> **AUROC 0.983** vs edit distance's 0.693

That looks like a win. Then we tested it against SPR — a different lab,
a different measurement principle:

> **Spearman +0.16, CI [−0.03, +0.33]** in the legal zone.
> Edit distance gets **+0.23** on the same data.

The explanation is in the 27 sequences measured *both* ways:

| FACS class | n | median KD |
| --- | --- | --- |
| high | 16 | **7.20 nM** |
| medium | 6 | **4.20 nM** |
| low | 5 | 15.37 nM |

"Medium" binds tighter than "high". The FACS class predicts measured
affinity with **AUROC 0.61** — barely above chance. Within "high", KD
spans 0.09 to 174 nM, a factor of 2000.

**The net has saturated its label.** It hits FACS class at 0.98, while
the label itself hits affinity at 0.61.

Be precise about what this is and is not. AUROC 0.61 is **not** a
mathematical ceiling. A noisy label is not an information-theoretic
wall: a model can in principle pick up features that correlate with
affinity beyond what the label itself captures, and we cannot rule out
that a different architecture or training setup would manage it.

What we can say is what we measured: **our** FACS-trained net did not
transfer. Spearman -0.04 in the legal zone, against edit distance's
+0.23. We found no transfer, and we have no reason to assume it exists
until someone demonstrates it.

FACS measures whether something displays on a cell surface and binds at
all. That is not the same as how tightly it binds.

### What we did about it

We do **not** rank the submission by the net. We use it for what it can
do — screen out probable non-binders — and rank the survivors by what
actually transferred to SPR: proximity to trastuzumab, plus conservative
chemical substitutions.

That is a weaker ranking. Spearman +0.23 is not much. We say so in the
pitch rather than hiding it behind a prettier number.

> **CORRECTED 14 November — both halves of that ranking were wrong.**
> See section 11 and `results/audit/`.
>
> **Proximity to trastuzumab cannot order our designs.** Every one of the
> twenty sits at exactly 4 CDRH3 mutations, so edit distance is *constant*
> across them. The one signal that transferred contributes nothing to the
> order of our own list. Restricted to measured variants at exactly 4
> mutations (n = 34), no scorer we have beats chance.
>
> **Conservative substitutions point the wrong way.** `conservative_fraction`
> carried the largest weight in the formula and has a *negative* point
> estimate against measured affinity in every subset tested: −0.133 (all),
> −0.162 (legal), −0.159 (d=4), −0.277 (d=5). No interval excludes zero, so
> this is not established — but it was never validated before it went in.
>
> The corrected number for the net in the legal zone is **+0.156, CI
> [−0.027, +0.326]**, not −0.04. The −0.04 was the d=5 slice quoted as if it
> were the whole legal zone.

## 4. We leaked first, and caught it ourselves

Our first split grouped on the first five residues. It gave AUROC 0.995.

Then we measured the distance from test to train: **194 of 300 test
sequences sat one mutation from a training sequence.** Median distance 1.
The library is a near-exhaustive mutagenesis — all ten positions carry
all twenty amino acids — so prefix grouping separates nothing.

Replaced with a distance split that removes everything within 2
mutations of the test set. Training shrank from 334,000 to 111,000. The
number fell from 0.995 to 0.983, and now it means something.

Worth saying out loud: **a near-perfect score on a biological assay is a
warning sign, not a victory.**

## 5. A third of a naive proposal set is illegal

Of 20,000 generated variants, **6,639 (33 %)** sat too close to a
provided sequence and would be rejected by the official checker.

The task description warns about this in its own RFD4 weights table: the
HER2 fine-tuned weights reproduce known binders, and `her2-2026-11-07`
on the conservative setting passes the novelty filter **0 %** of the time.

This is why the gate sits *before* the ranking in our pipeline, not
after. Anyone generating candidates should run them through it before
spending GPU time:

```
python scripts/check_novelty.py
```

## 6. The pitch is a technical decision, not a side task

Base rate for beating trastuzumab on affinity: 3.3 %. Designs tested
depends on pitch placement:

| Designs tested | P(at least one beats affinity) | |
| --- | --- | --- |
| 4 | **13 %** | every team gets this |
| 8 | 24 % | 3rd place |
| 12 | 33 % | 2nd place |
| 16 | **42 %** | 1st place |

S-tier needs all three axes at once:

| | 4 designs | 16 designs |
| --- | --- | --- |
| thermostability left to chance | 5.8 % | 21 % |
| thermostability engineered | 9.2 % | **32 %** |

We cannot move the 3.3 % base rate much. We can quadruple the number of
tickets. **Winning the pitch is the single largest lever we have on the
experimental outcome.**

---

## The two novelty rules, which are separate

For the optimisation track, the official checker applies **two** independent checks.
A design failing either is rejected — never reclassified into the other
track.

| | Rule | For a 10-mer |
| --- | --- | --- |
| **A** | CDRH3 identity to *trastuzumab* under 70 % | ≥ 4 mutations |
| **B** | More than 2 point mutations from *every provided data point* | ≥ 3 from all 525,063 |

Provided data = the union of both CSVs: 524,346 FACS + 754 unique SPR,
overlapping in 37, giving exactly **525,063** — matching the data README.

Note trastuzumab itself is in the SPR table, so rule B alone already
forces ≥3 from the original. Rule A forces ≥4.

**Our 20 designs: 20/20 on both.**

### One ambiguity worth asking the organisers about

The task text contradicts itself on rule B:

- one place reads as *at least* two mutations from any provided antibody → ≥2
- another reads as *strictly more than* two point mutations → ≥3

We implement the strict reading (≥3). If the real rule is ≥2, we are
leaving candidates on the table.

---

## A convention that can silently break everything

The dataset calls CDRH3 a **10-mer**. The heavy chain is built as

```
...YYCSR + <10-mer> + Y + WGQGT...
```

The trailing Y is **not** included, even though Kabat and IMGT count it.
Trastuzumab's loop is `WGGDGFYAMD` here, not `WGGDGFYAMDY`.

The novelty filter is measured against these tables and must be computed
in their convention. Mix the two and every distance is off by one
residue, and a design that looks novel gets rejected.

Tested in `test_datasettets_cdrh3_er_en_rest_kortere`.

## 7. Co-folding confidence is the fourth ranking signal to lose to edit distance

We put every candidate ranking signal through the same test: take variants
whose KD was **measured**, score them, and see whether the score orders them.
A signal that cannot order known answers cannot order our designs.

The panel is 24 SPR variants spread evenly across the measured range
(0.023 nM to 13.9 uM) plus trastuzumab, all with a 10-residue CDRH3 so that
they share one template. Each one is co-folded with HER2 domain IV as a
four-chain complex, same template, same seed, same settings for every model.

Boltz-2, n = 25, Spearman against -log10(KD), 95 % bootstrap CI:

| Signal | Spearman | 95 % CI | AUROC (tighter half) |
| --- | --- | --- | --- |
| interface ipTM, antibody-antigen | +0.41 | [0.00, 0.72] | 0.61 |
| interface ipTM, VH-antigen only | +0.41 | [0.00, 0.73] | 0.61 |
| ipTM (whole complex) | +0.43 | [0.02, 0.72] | 0.63 |
| pTM | +0.44 | [0.02, 0.73] | 0.63 |
| ranking score | +0.30 | [-0.07, 0.63] | 0.61 |
| pLDDT | +0.18 | [-0.22, 0.51] | 0.59 |
| **edit distance to trastuzumab** | **+0.51** | **[0.16, 0.78]** | **0.72** |
| our FACS-trained net | +0.22 | [-0.23, 0.57] | 0.59 |

OpenFold3 (OpenBind weights), same 25 variants, same template, same seed:

| Signal | Spearman | 95 % CI | AUROC |
| --- | --- | --- | --- |
| interface ipTM, antibody-antigen | +0.24 | [-0.18, 0.61] | 0.58 |
| ipTM (whole complex) | +0.15 | [-0.23, 0.51] | 0.56 |
| pLDDT | -0.02 | [-0.41, 0.40] | 0.50 |

Read this carefully, because it cuts both ways.

- Boltz-2's interface confidence **does** carry affinity signal. +0.41 is not
  nothing, and it is better than the net we trained on half a million FACS
  labels. A folding model that never saw this assay beats a model that was
  fitted to it. That is finding 3 again, from a different direction: training
  on the FACS label bought us nothing for affinity, while a model that never
  saw the label carries some.
- But on the same 25 variants, **edit distance still wins**, and its interval
  is the only one clearly clear of zero. This is the fourth ranking method we
  have tried and the fourth to lose to counting mutations.
- pLDDT is close to useless here (+0.18). Whole-structure confidence is not
  interface quality. If anyone ranks on pLDDT, stop them.

- OpenFold3 is weaker than Boltz-2 on every metric and its intervals all
  contain zero. On this task, with this template, the two models are not
  interchangeable.

The honest summary: **co-folding confidence is a weak second opinion, not a
ranking function.** We use it to flag designs the folding model dislikes, and
we keep ranking on distance plus conservative substitution.

### And the two models do not agree with each other

The strongest argument against ranking on co-folding is not any single
correlation. It is that the two models, given identical inputs, rank our
20 designs almost independently:

| Metric | Boltz-2 vs OpenFold3, Spearman, n = 20 |
| --- | --- |
| interface ipTM (ab-ag) | +0.18 |
| pLDDT | +0.34 |

The absolute scales are not comparable either. Trastuzumab's own complex
scores interface ipTM **0.94 in Boltz-2 and 0.21 in OpenFold3** — the same
four chains, the same template, the same crystal structure behind it. A
threshold carried from one model to the other is meaningless, and so is any
"consensus" built by averaging them.

This is the practical rule: **compare within one model, never across, and
treat agreement between models as the evidence, not either score alone.**
Here there is very little agreement, so there is very little evidence.

## 8. Boltz-2 does not think any of our designs beats trastuzumab

Same run, same template, same seed. Trastuzumab is in the panel as the
baseline, scored by exactly the same pipeline as the designs:

| Model | Metric | Trastuzumab | Designs, median | Designs, best | Above baseline |
| --- | --- | --- | --- | --- | --- |
| Boltz-2 | interface ipTM (ab-ag) | 0.941 | 0.916 | 0.936 | **0 of 20** |
| Boltz-2 | ipTM | 0.960 | 0.943 | 0.956 | **0 of 20** |
| Boltz-2 | pLDDT | 97.2 | 96.6 | 97.6 | 4 of 20 |
| OpenFold3 | interface ipTM (ab-ag) | 0.206 | 0.121 | 0.302 | 2 of 20 |
| OpenFold3 | ipTM | 0.599 | 0.554 | 0.632 | 2 of 20 |
| OpenFold3 | pLDDT | 87.7 | 86.8 | 88.9 | 3 of 20 |

Boltz-2 puts no design above trastuzumab on either interface metric.
OpenFold3 puts two above. Only one of those two is near the top of the
Boltz-2 ordering; the other sits mid-table there. Across the top five the
models share three designs out of five — some signal, not a ranking.

Given section 7 this is **not** a prediction that the designs are worse. The
metric that produced this table orders known answers at rho 0.41 with an
interval touching zero, so it cannot resolve a 0.025 ipTM difference. What it
does tell us is that nothing in our set stands out as obviously better than
the clinical antibody, which is the expected result from finding 1: only
3.3 % of measured variants beat trastuzumab, and the novelty rule keeps us
out of the zone where the good ones live.

It is also a warning about the pitch. If we claim a design is better, this is
the first number a judge with the same tooling will compute, and it will not
support us. The claim we can support is about thermostability and epitope.

## 9. Our disulfide scan could not see glycine, so it missed the known answer

> **DO NOT ACT ON THE RECOMMENDATION IN THIS SECTION.**
> The VH44–VL100 disulfide it proposes is **rejected by the official checker** in both
> submission formats, with `errors: ['unpaired_cysteine']`. A design carrying
> it is never sent to the lab. The literature does not support it either:
> Weatherill et al. 2012 (*PEDS* 25(7):321–329) tested it across six scFvs
> and found "no or little increase in thermal stability".
>
> Full result and the replacement in `results/disulfide/` and section 11.
> The bug analysis below is still correct and worth reading — a C-beta the
> scan could not see is a real lesson. The conclusion it led to is not.

A bug worth writing down, because it was invisible and the fix changed our
recommendation.

The scan looks for positions where a new disulfide between VL and VH is
geometrically possible, measured in the real structure. It required a
C-beta atom at both positions. **Glycine has no C-beta.** Every glycine was
therefore skipped in silence — and the most-used engineered VH-VL disulfide
in the literature (Kabat VH44 + VL100) sits on a glycine in the heavy chain.

Two fixes:

1. Build the C-beta where it would be, from N, CA and C. Tested against every
   residue in 1N8Z that has a real one (938 of them): median error 0.08 A,
   90th percentile 0.19 A, worst 0.33 A.
2. The acceptance window was too tight. The known pair sits at C-beta-C-beta
   5.18 A; our cutoff was 4.8 A. The literature's own geometric screens use a
   CA-CA window of roughly 4.2-6.6 A, and side chains rotate when both
   positions become cysteine, so the wild-type C-beta distance is a cruder
   measure than it looks. Widened to 5.5 A.

With both fixes the scan returns three candidates instead of one, and the
known pair is among them.

**This changes our recommendation.** We now propose the literature pair:

| | Our geometric pick | The literature pair |
| --- | --- | --- |
| Light chain | Phe98 -> Cys | **Gln100 -> Cys** |
| Heavy chain | Leu45 -> Cys | **Gly44 -> Cys** |
| C-beta-C-beta | 4.38 A | 5.18 A |
| CA-CA | 6.37 A | 5.46 A |
| Built in a lab | no | **yes, on 4D5** |

4D5 is the antibody trastuzumab was humanised from, and this exact bond has
been made on it (Woern & Plueckthun, *Biochemistry* 1999). Reported: a clear
gain in resistance to denaturation, binding essentially unchanged, and a
sharp drop in expression yield. The same pair is first on the ranked list in
the Pastan/Brinkmann disulfide-stabilised Fv patents, and our own pick
(VH45 + VL98) is eighth on that same list -- so the scan was not wrong, just
not reading the literature.

Numbering note that cost us an hour: these are Kabat positions. For
trastuzumab's light chain Kabat and sequential numbering coincide, so Kabat
VL100 is sequential 100. For the heavy chain they coincide up to position 52,
so Kabat VH44 is sequential 44. In 1N8Z both are the plain residue numbers.
Do not copy "Gly100Cys" from a paper: in trastuzumab's humanised FR4 that
position is **Gln**, not Gly.

Both positions are framework, so neither novelty rule is touched, and two new
cysteines keep the parity even.

---

## 10. The way out: ask the models a question that has a known answer

Section 7 left us with three folding models that disagree. The usual answers
to that — average them, or trust the most confident — are both guesses. There
is a better one available here, because this target has a crystal structure.

Trastuzumab's epitope on HER2 is known from 1N8Z: 19 HER2 residues with a
heavy atom within 4.5 A of the antibody. We can ask each model to predict the
trastuzumab complex and then simply check how much of that epitope it
recovers. No confidence score involved — just which residues ended up in
contact.

| Model | Trastuzumab's predicted epitope vs the crystal (Jaccard) |
| --- | --- |
| **Boltz-2** | **1.00** |
| Boltz-1 (no template) | 0.29 |
| Protenix v2 | 0.26 |
| OpenFold3 OpenBind | 0.10 |

The contact cutoff does not decide this. At 4.0, 4.5 and 5.0 A, Boltz-2
scores 0.94, 1.00 and 1.00, while no other model gets above 0.30 at any of
the three. ESMFold2, which folds from sequence alone with no template, sits
at 0.06-0.07: it does not assemble the complex at all and cannot be used for
this target.

Boltz-2 reproduces the crystal epitope exactly. The other three do not come
close. **That is the tie-break we were missing:** we are not choosing between
three equally plausible opinions, we are choosing the one model that can get
the one case we can check right.

Caveat, and it matters: every templated model was given the same template,
which already places HER2 and the framework. So this measures whether a model
stays consistent with a structure it was handed, not whether it could find the
epitope unaided. It is still a real discriminator — three models were handed
the same help and lost the epitope anyway — but it is not evidence that
Boltz-2 would find a novel epitope.

### And it gives us a computed number on the second scored axis

Epitope preservation is one of the three things the lab measures. We do not
measure it. What follows is a **prediction**: contact overlap computed from
Boltz-2's predicted structures, one model, one seed.

| | Predicted overlap with trastuzumab's epitope (Boltz-2, single seed) |
| --- | --- |
| Our 20 designs, median | **0.95** |
| Our 20 designs, worst | **0.85** |

Every design is predicted to keep essentially the whole epitope. That is the
expected result for four point mutations in CDRH3 on an unchanged framework.

**It replaces an assumption with a model's opinion, not with a measurement.**
The actual epitope result is the lab's competitive binding assay against
trastuzumab, and it arrives after the event. Three limits on the number above,
all of which matter:

- **One seed.** We have no replicate, so we cannot say how much of the 0.95
  is signal and how much is the sampling noise of a single run. An overlap
  computed from one structure has no error bar.
- **Templated.** Boltz-2 was handed a template that already places HER2 and
  the framework. A high overlap partly reflects that it kept what it was
  given.
- **Self-referential.** The reference epitope is Boltz-2's own prediction of
  the trastuzumab complex, which happens to reproduce the crystal exactly.
  We are comparing the model against itself on the designs.

Taken together: this number supports the claim that our mutations are
unlikely to move the binding site, which is what four substitutions in one
loop on an unchanged framework should do. It is not evidence that the epitope
is preserved, and it must not be written up as one.

The measure behaves the way a measure of a *different* axis should: across
the calibration variants it has no relationship with measured affinity
(Spearman -0.02), and all of them sit between 0.86 and 1.00. Variants of
trastuzumab bind the same patch whether they bind it tightly or weakly.
Epitope and affinity are separate questions, and this number answers the
first one only.

---

## 11. The gate is verified, against the provided data, with controls

The official checker has now run on all twenty designs. **All twenty pass**, in both
`paired` and `scfv` format.

The bundled reference does **not** contain the 525,063 provided CDRH3s — for
every design it reports the nearest known loop as trastuzumab itself. So the
rule that actually constrains us was untested. We built a reference from the
provided data (every sequence within 4 mutations of a design, plus 400 random
controls: 4,579 rows) and re-ran.

> Every design sits at **exactly 3** mutations from the nearest provided
> sequence. The requirement is 3. **Zero margin, uniformly.**

All twenty still pass. And because a gate that passes everything is
indistinguishable from a gate that is switched off, we planted three designs
one mutation from a sequence *in* the reference. In the same run, `--drop`
removed exactly those three and kept all twenty:

```
[novelty] control_d1_0_0: CDRH3 SRCDTTACTLWIY is 1 mutation(s) from
          SRADTTACTLWIY in provided_ref; a novel design needs at least 3
```

Three things this run settled that the written rules did not:

- **The enforced rule is 3 mutations, not 4.** `required_cdrh3_mutations: 3`.
  We had read the 70 % identity clause as implying 4. It does not; the
  identity threshold decides the *track*. Our designs carry 4 anyway, so they
  satisfy both readings — but the margin is ours, not required.
- **CDRH3 is 13 residues** in the official checker's numbering (`SRWGGDGFYAMDY`), not the
  10-mer core `src/abx/` works on. Never quote an identity percentage without
  saying which length it is over.
- **`novelty-report` does not report the reference comparison.** With
  `--reference-db` it still prints `cdrh3_nearest_source: "trastuzumab"`,
  even for a control sitting 1 mutation from a reference sequence. Only
  `team-validate` enforces the provided-data check. This nearly produced a
  false alarm.

## 12. Developability: a ranking signal that is not a model we fitted

The official checker ships a 13-motif liability reference. Run against it, trastuzumab
itself carries **8 liabilities, 3 of them high severity**, including a `DG`
isomerisation motif inside CDRH3. A licensed drug manufactured since 1998
carries that motif. Liabilities are a risk budget, not a disqualification,
and a raw count means nothing — the delta against trastuzumab does.

| | |
| --- | ---: |
| designs **worse** than trastuzumab on high-severity motifs | **0** |
| the same | 16 |
| **better** | **4** |

The four that improve — `tz-opt-02`, `18`, `07`, `15` — all drop the `DG`
motif. But **nine designs introduce `FragH:DP`**, an acid-labile Asp-Pro bond
that trastuzumab does not carry. Several look neutral on the count because
they traded `DG` for `DP`. That trade is not neutral: `DG` is demonstrably
tolerated in a licensed product; `DP` is risk we would be introducing with no
such evidence. Those nine are excluded from the submission.

## 13. Thermostability: the disulfide is dead, and something better replaces it

The VH44–VL100 plan from section 9 fails twice over.

**The gate rejects it.** Four designs were built with and without the
substitutions (positions located with `immunum`, not assumed: Kabat VH44 = G,
VL100 = Q in all four). Every variant carrying the bond fails:

| Format | Parents | With disulfide | Failure |
| --- | ---: | ---: | --- |
| `paired` | 4 pass | **4 fail** | `heavy_cysteines_are_paired` |
| `scfv` | 4 pass | **4 fail** | `heavy_domain_cysteines_are_paired` |

scFv does not rescue it. The official checker pairs cysteines **per domain**, not per
chain, so the VH domain still holds three. **This generalises: any engineered
disulfide placing an odd number of cysteines in one variable domain is
rejected.** An interdomain staple is by construction exactly that.

**The literature does not support it either.** Weatherill et al. 2012 tested
it across six scFvs in both orientations — the only systematic generality
test — and found "no or little increase in thermal stability", though it was
the best of six positions. The nearest 4D5-relevant measurement (Wörn &
Plückthun 1999, *Biochemistry* 38(27):8739–8750) bought ~0.4 M GdnHCl of
chemical stability for a **29-fold** collapse in expression yield; Zhao et al.
2010, using exactly G44C/Q100C, saw a 4-fold drop.

Two misattributions circulate and should not be repeated. The widely quoted
"Tm 60→69 °C" belongs to **Young et al. 1995**, antibody Se155-4, positions
A-L57/D-H106 — a different bond on a different antibody. And Okazaki 2023 is
often cited as disulfide work; it contains none.

### The replacement: VL Arg66Gly

> Okazaki et al. 2023, *ACS Omega* 8(25):22945–22954, PMID 37396255.
> Trastuzumab scFv by DSF: **61.8 ± 0.2 °C → 66.7 ± 0.2 °C**, ΔTm **+4.9 °C**,
> binding unchanged.

Measured on our scaffold, on the judged axis, by direct measurement. It
touches only light-chain framework: CDRH3 is untouched so novelty is
unaffected, and it adds no cysteine. The gate run confirms liability counts
are **identical** with and without it, and all variants pass.

It is one paper and one mutation, and we have not measured it. So every
stabilised candidate is submitted alongside its unmodified parent.

## 14. The submission

The candidate list, its composition and the per-design rationale are kept
out of this repository by the team's sharing rules. What is shareable is
the reasoning that produced it, in sections 11-13 above: the gate is
verified with planted controls, nine designs are excluded for introducing a
liability trastuzumab does not carry, and the thermostability axis rests on
VL R66G rather than the disulfide.

---

## What we have not done

**Thermostability, beyond the one mutation.** Section 13 replaced the
disulfide with VL R66G, which is gate-legal and has a measured ΔTm on this
scaffold. What is still missing is a **second, independent stabilising
change**, so that the thermostability axis does not rest on a single
published number from a single paper. We have not folded the R66G variants
either — worth doing to confirm the CDR loops do not move, though R66 is
framework and distant from the paratope, so the prior is that they do not.

**Interface PAE.** Section 7 ranks on ipTM and pLDDT, read from the summary
scores each model writes. Interface PAE needs the detailed score files, which
are CBOR, and the reader for them is not installed in our environment. It is
the metric the organisers' own pass criteria use, so it should be added.

**Calibration against the shuffled pairs.** Section 7 calibrates against
measured KD, which is the stronger test, but we have not measured the false
positive rate on the Smorodina true negatives at the thresholds we use.

**Replicates.** Everything in sections 7 and 8 is one seed and one diffusion
sample per design. We therefore do not know the run-to-run spread, and we
should not separate two designs on a difference we have not shown is larger
than the noise. A second seed is the cheapest thing left to run.

**Epitope, on more than one model.** Section 10 measures epitope overlap, but
only Boltz-2 passed the crystal check, so the number rests on a single model.
It also rests on a template that already placed HER2. An apo prediction, or a
model that recovers the epitope without the template, would make it stronger.

~~**The official checker itself.**~~ **Done** — section 11. All twenty pass, against a
reference built from the provided data, with planted controls proving the
filter is on.

---

## What we think should happen next

### Thermostability: resolved differently than planned

The four-step disulfide plan that stood here is obsolete — step 1 produces a
molecule the gate rejects (section 13). It has been replaced by VL R66G,
which is already built, gate-checked and in the submission.

What remains on this axis is one thing: **a second, independent stabilising
change**, so the axis does not rest on a single published ΔTm. That is a
literature and design question, not a compute question.

The expression-yield judgement that was step 4 no longer applies — R66G has
no reported yield penalty, which was part of why it won over the disulfide.

### The two questions this document said would decide the de novo track

Both are answered in the task description, and neither needed to be asked:

- **How the de novo novelty rule is computed.** Not whole-sequence. The official checker
  does a *type-aware* comparison: a VHH against known VHHs and heavy chains,
  an scFv or Fab split into paired chains and compared as paired antibodies.
  The de novo threshold is under 70 % identity to the nearest known antibody,
  plus at least 2 CDRH3 mutations from trastuzumab.
- **Whether the four guaranteed slots are per team or per track.** **Per
  team.** "Each team gets 4 designs whatever happens." Placement adds 12/8/4,
  and each of 6 judges nominates a team for 2 more, to a maximum of 16.

So splitting across tracks does halve our guaranteed shots. The prize fund is
split equally between tracks, which pulls the other way.

### Both tracks, or only optimisation?

Everything we have built so far is optimisation. Our recommendation is to
**keep optimisation as the primary track, and enter de novo only if two
things are confirmed by the organisers:**

- **How the de novo novelty rule is computed.** "Under 70 % identity with a
  known antibody" is straightforward if it is measured on the CDRs. If it is
  measured on the whole sequence, then nothing we can build in the time left
  qualifies, because every real antibody framework is far more than 70 %
  identical to known antibodies. This single answer decides the track.
- **Whether the four guaranteed lab slots are per team or per track.** If
  they are per team, splitting two and two halves our shots in a track where
  only 3.3 % of measured variants beat the baseline (section 1). If they are
  per track, entering both is close to free.

If both answers come back favourably, the argument for de novo is not
affinity. It is that the tier system gives credit for **any one** axis, and a
fresh scaffold has a genuine chance on thermostability where a trastuzumab
variant is fighting an already-good clinical antibody. That makes de novo a
cheap ticket on the axis we are strongest on, not a second attempt at the
axis we are weakest on.

We have queued a de novo experiment that redesigns all six CDRs against HER2
on an Fv framework, to find out whether such designs pass our own gate at
all. It is there to answer the question, not to be submitted.

