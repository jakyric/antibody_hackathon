# The ranked list: 20 designs on two axes

`results/ranked_submission.csv`, produced by `scripts/rerank.py`.

Teams submit a ranked list of up to 20 designs. Every team gets its top 4
tested whatever happens; pitch placement and judge nominations add up to 12
more. **The order is what matters, not the length.**

## The two axes

| Axis | What it is | How good it is |
| --- | --- | --- |
| **Binding** | closeness to trastuzumab, the original ranking | weak: Spearman +0.23 against measured SPR in the legal zone |
| **Developability** | liability motifs vs trastuzumab, from the official checker's LAP reference | clean, but not a prediction of any judged property |

Neither is strong. The binding proxy is the only signal we found that
transferred to measured affinity at all, and it transferred weakly. The
developability axis is sound, but it speaks to chemical degradation and
expression risk, not to the three properties Adaptyv measures.

## How they are combined

**Rank aggregation**: each design's position on each axis, added. No weights,
no units, no invented precision. A design that is 2nd on binding and 1st on
developability (sum 3) outranks one that is 1st and 5th (sum 6).

Designs that introduce a high-severity motif trastuzumab does not carry are
demoted below everything else, whatever their sum.

### What we tried first and rejected

Sorting strictly by developability tier. It promoted designs from binding rank
15 and 18 straight into the top four on the strength of a liability count.

That is over-trusting the proxy. Liability motifs are not a thermostability
prediction, and thermostability is the judged axis. Trading away the only
signal we have for a judged property, in favour of one that speaks to a
different concern entirely, is not a trade the evidence supports.

Rank aggregation keeps both axes visible and lets neither dominate.

## The result

| New | Design | Bind | Dev | Sum | Move | CDRH3 |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| **1** | **tz-opt-02** | 2 | 1 | 3 | +1 | MGPEPFYAMD |
| **2** | **tz-opt-01** | 1 | 5 | 6 | −1 | APGDGFQIMD |
| **3** | **tz-opt-03** | 3 | 6 | 9 | — | LGGDGAYMVD |
| **4** | **tz-opt-07** | 7 | 3 | 10 | +3 | FGGWPWYAMD |
| 5 | tz-opt-04 | 4 | 7 | 11 | −1 | WPGDGAYKME |
| 6 | tz-opt-05 | 5 | 10 | 15 | −1 | WGGDGWYNAE |
| 7 | tz-opt-11 | 11 | 8 | 19 | +4 | WGPDGAYNWD |
| 8 | tz-opt-15 | 15 | 4 | 19 | +7 | FGGDMWYWMD |
| 9 | tz-opt-18 | 18 | 2 | 20 | +9 | WGGWPAYALD |
| 10 | tz-opt-20 | 20 | 9 | 29 | +10 | LGADGVYMMD |
| 11 | tz-opt-19 | 19 | 11 | 30 | +8 | WGPDGADWMD |
| 12–20 | the `FragH:DP` group | | | | | demoted |

**tz-opt-02 leads because it is the only design strong on both axes**: second
on the binding proxy, and the cleanest developability profile of the twenty
(7 liabilities against trastuzumab's 8, one fewer high-severity, nothing new
introduced). It drops trastuzumab's `DG` isomerisation motif and its CDRH3
tryptophan, and pays only a second methionine.

The top four changed by one member against the binding-only ranking:
`tz-opt-04` out, `tz-opt-07` in.

## What this ranking does not claim

These designs are **expected** to bind. Nothing here is a measurement. The
affinity signal behind the binding axis is weak enough that the difference
between rank 1 and rank 5 should not be read as meaningful; the difference
between the top group and the demoted `FragH:DP` group is the part we would
defend.

Beating trastuzumab on affinity is rare: of 182 published, professionally
designed variants with measured KD, 6 beat its median and 1 beat its best. We
are not predicting that any of these will. We are ordering twenty candidates
by how good a bet each looks on the evidence we have.

## Reproduce

```
abguard team-validate <table> --all-liabilities -o liabilities.json
python scripts/liability_analyse.py --liabilities liabilities.json \
       --out liability_delta.csv
python scripts/rerank.py --designs results/designs.csv \
       --liabilities liability_delta.csv --out results/ranked_submission.csv
```

Entry-gate evidence for all 20 is in `results/gate/`; the developability
analysis is in `results/developability/`.

---

> **SUPERSEDED, 14 November 2026.** The rank aggregation below assumed the
> binding axis could order the twenty. The audit in `results/audit/` showed it
> cannot: all twenty sit at exactly 4 CDRH3 mutations, so edit distance is
> constant across them and orders nothing, and inside that stratum no scorer
> beats chance. The term that actually dominated the old formula,
> `conservative_fraction`, has a negative point estimate against measured
> affinity in every subset tested.
>
> The live list is **`results/SUBMISSION.md`** and `results/SUBMISSION.csv`.
> This file is kept because the reasoning that was wrong is part of the record.
