# The final ten: what we submit and why

`results/SUBMISSION.csv`. Built by `scripts/build_final.py`, validated with
`abguard team-validate --reference-db` against a reference built from the
provided HER2 data. All ten pass. Trastuzumab, added as a control, fails.

## The three facts that determined this list

**1. We cannot rank these on affinity, and we will not pretend to.**
All twenty designs sit at exactly 4 CDRH3 mutations from trastuzumab. Edit
distance — the only scorer that tracks measured affinity (rho +0.226 in the
legal zone, interval clear of zero) — is therefore *constant* across them and
orders nothing. Inside that stratum (n = 34 measured variants) no scorer we
have beats chance. Full working in `results/audit/`.

**2. Developability can rank them.** Measured against trastuzumab with the
organisers' own liability reference. Nine of twenty introduce `FragH:DP`, a
high-severity acid-labile motif trastuzumab does not carry; they are excluded.
Details in `results/developability/`.

**3. Only a stabilised design can reach S-tier.** S-tier needs all three
properties improved. Four substitutions in one CDR loop give no reason to
expect a higher melting temperature, so an unmodified design tops out at
A-tier and realistically B. We need a stabilising change to compete at all on
the third axis.

## The stabilising mutation: VL Arg66Gly

> Okazaki et al. 2023, *ACS Omega* 8(25):22945–22954, PMID 37396255.
> Trastuzumab scFv by DSF: **61.8 ± 0.2 °C → 66.7 ± 0.2 °C, ΔTm +4.9 °C**,
> binding unchanged.

Why this one and not the obvious candidate:

| | VH44–VL100 disulfide | VL R66G |
| --- | --- | --- |
| passes the gate | **no** — `unpaired_cysteine`, both formats | **yes**, cleanly |
| measured on trastuzumab scFv | no | **yes** |
| effect on Tm | "no or little increase" across six scFvs | **+4.9 °C measured** |
| expression risk | 29-fold and 4-fold yield drops reported | none reported |
| adds cysteines | two | none |

The disulfide is dead twice over: the official checker rejects it outright (see
`results/disulfide/`), and Weatherill et al. 2012 (*PEDS* 25(7):321–329)
tested it across six scFvs and found no meaningful thermal gain. The nearest
4D5-relevant measurement (Wörn & Plückthun 1999) bought ~0.4 M GdnHCl of
chemical stability for a **29-fold** collapse in expression yield.

R66G touches only light-chain framework. It does not alter CDRH3, so our
novelty margin is untouched, and it adds no cysteine. The gate run confirms
liability counts are identical with and without it.

**It is still one paper and one mutation.** We have not measured it. Every
stabilised candidate is paired with its unmodified parent below so the lab
result can tell them apart.

## The list

| # | ID | CDRH3 | Stab | Why this one |
| ---: | --- | --- | :-: | --- |
| 1 | `tz-opt-02-stab` | MGPEPFYAMD | ✓ | best developability of all twenty: drops both trastuzumab's `DG` isomerisation motif and its CDRH3 tryptophan, adds nothing new. Stabilised, so all three axes are in play. |
| 2 | `tz-opt-18-stab` | WGGWPAYALD | ✓ | second-best developability, a clearly distinct loop from #1. |
| 3 | `tz-opt-07-stab` | FGGWPWYAMD | ✓ | removes the `DG` motif, loop distinct from both above. |
| 4 | `tz-opt-15-stab` | FGGDMWYWMD | ✓ | fourth design improving on trastuzumab's liability profile. |
| 5 | `tz-opt-01-stab` | APGDGFQIMD | ✓ | top loop under the old binding proxy — kept as a hedge in case that proxy carries signal we could not demonstrate. |
| 6 | `tz-opt-03-stab` | LGGDGAYMVD | ✓ | second loop from the old proxy, same hedge, different sequence. |
| 7 | `tz-opt-02` | MGPEPFYAMD | — | **control** for #1: isolates what R66G did. |
| 8 | `tz-opt-18` | WGGWPAYALD | — | **control** for #2. |
| 9 | `tz-opt-07` | FGGWPWYAMD | — | **control** for #3. |
| 10 | `tz-opt-01` | APGDGFQIMD | — | **control** for #5. |

Six distinct CDRH3 loops. Pairwise distance between them: **min 3, median 6,
max 8** — these are different bets, not one bet written six ways.

## Why the stabilised versions take the guaranteed slots

Every team gets its **top four** tested whatever happens. Those four are all
stabilised and all carry different loops: four independent shots at improving
all three properties.

The controls sit at 7–10 deliberately. If we place well enough for them to be
made, we learn what R66G did on our own molecules — which is worth having
whichever way it goes. If we do not place, we lose nothing, because the
guaranteed four were never the controls.

## What we are not claiming

These designs are **expected** to bind. Nothing here is a measurement.

- **Affinity:** we have no validated way to order these ten. Of 182 published,
  professionally designed variants with measured KD, 6 beat trastuzumab's
  median and 1 beat its best. We are not predicting that any of ours will.
- **Epitope:** the 0.95 median overlap in `FINDINGS.md` is a *prediction* from
  one Boltz-2 run at one seed, against a reference the same model produced.
  It supports the claim that four substitutions in one loop should not move
  the binding site. It is not a measurement, and the real answer is Adaptyv's
  competitive binding assay.
- **Thermostability:** +4.9 °C is Okazaki's number on their molecule, not ours.

## Reproduce

```
python scripts/build_final.py --designs results/designs.csv --out SUBMISSION.csv
abguard team-validate SUBMISSION.csv --reference-db <provided reference>
```

Gate evidence `results/gate/` · liabilities `results/developability/` ·
ranking audit `results/audit/` · disulfide `results/disulfide/` ·
stabilisation gate run `results/stabilise/`.
