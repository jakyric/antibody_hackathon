# What we measured, and what it means for the submission

Everything below is computed from the provided data with the code in this
branch. Each section says which script produces it, so every number can be
re-run rather than taken on trust.

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
provided sequence and would be rejected by `abguard`.

The task description warns about this in its own RFD4 weights table: the
HER2 fine-tuned weights reproduce known binders, and `her2-2026-11-07`
on the conservative setting passes the novelty filter **0 %** of the time.

This is why the gate sits *before* the ranking in our pipeline, not
after. Anyone generating candidates should run them through it before
spending GPU time:

```
python scripts/check_novelty.py   # not in this branch: it reads the provided data
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

For the optimisation track, `abguard` applies **two** independent checks.
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

### And it gives us a number on the second scored axis

Epitope preservation is one of the three things the lab measures, and until
now we had only an assumption about it. With Boltz-2 we have a measurement:

| | Overlap with trastuzumab's epitope, same model, same seed |
| --- | --- |
| Our 20 designs, median | **0.95** |
| Our 20 designs, worst | **0.85** |

Every design keeps essentially the whole epitope. That is the expected result
for four point mutations in CDRH3 on an unchanged framework, and now it is
measured rather than assumed.

The measure behaves the way a measure of a *different* axis should: across
the calibration variants it has no relationship with measured affinity
(Spearman -0.02), and all of them sit between 0.86 and 1.00. Variants of
trastuzumab bind the same patch whether they bind it tightly or weakly.
Epitope and affinity are separate questions, and this number answers the
first one only.

---

## What we have not done

**Thermostability, beyond the one pair.** Section 9 gives a disulfide with
published lab support, and it costs nothing in novelty because it is in the
framework. What is still missing is everything else on that axis: a second,
independent stabilising change; a fold of the disulfide variant to check that
the bond is modelled closed and the CDR loops do not move; and a decision on
whether the known drop in expression yield matters for a lab that has to make
the protein before it can test it. This remains the cheapest of the three
axes and the one where we are most likely to beat the baseline.

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

**`abguard` itself.** Our gate is a local mirror, tested against known
answers, but theirs is the authority. It must run before submission.

---

## What we think should happen next

### Thermostability needs an owner, and the job is now small

It is the only one of the three axes the novelty rules do not constrain, the
one with the clearest published precedent, and therefore the one where we are
most likely to actually beat the baseline. Section 9 turned it from an open
question into a short list:

1. Build the VH Gly44Cys + VL Gln100Cys variant of our top designs.
2. Co-fold it against HER2 and check two things: that the two new cysteines
   come out within bonding distance of each other, and that the CDR loops do
   not move relative to the same design without the bond. Same model, same
   seed, parent as the control.
3. A second, independent opinion on stability (StaB-ddG, ESM-2 likelihood).
4. Decide whether the documented drop in expression yield is acceptable. A
   protein the lab cannot make cannot be measured.

Steps 1-3 are a couple of hours of compute. **Step 4 is a judgement call and
needs the person with wet-lab antibody experience, not the person running the
pipeline.** It should not be the same person as the design and ranking thread:
the tooling does not overlap, and both threads are time-boxed by the same
deadline.

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
