# The final ten: what we submit and why

`results/SUBMISSION.csv`. Six loops, each as a stabilised design, with the top
four also submitted unmodified as controls. Built by
`scripts/build_final_clean.py` from `results/pool/final_six.csv`.

**Checks, all green:**

| Check | Result |
| --- | --- |
| the official checker's team validation, reference built from the provided HER2 data | the ten alone **exit 0**; with trastuzumab added, 10 of 11 pass and the control is the one that fails |
| the official checker's own bundled reference, 610,109 known CDRH3s | every design **≥3 edits**; trastuzumab collides at 0 |
| exhaustive edit distance to all 525,063 provided loops | every design at **exactly 3**; control at 0 |
| our per-requirement check, every design and every rule | **10 pass, 0 fail, 0 unknown** |
| test suite | 201 passed |

Earlier lists are kept as `SUBMISSION_v1.*` and `SUBMISSION_v2.*`.

## The list

| # | ID | CDRH3 | Stab | Hypothesis | Contacts changed | Epitope | Liabilities |
| ---: | --- | --- | :-: | --- | --- | ---: | --- |
| 1 | `tz-clean-02-stab` | RGGEGFYANE | ✓ | tolerance | W99 | **1.00** | −3 / −1 high |
| 2 | `tz-grip-08-stab` | WGSEGFYNRD | ✓ | **contact** | **none** | **1.00** | −2 / −1 high |
| 3 | `tz-clean-01-stab` | YGGDYFYGKD | ✓ | tolerance | W99, G103 | 0.95 | −3 / −1 high |
| 4 | `tz-grip-02-stab` | WGSKGFYGPD | ✓ | **contact** | **none** | 0.95 | −2 / −1 high |
| 5 | `tz-clean-03-stab` | FGGGGNYAKD | ✓ | tolerance | W99, F104 | **1.00** | −3 / −1 high |
| 6 | `tz-clean-07-stab` | HGGDRFYKRD | ✓ | tolerance | W99, G103 | 0.95 | −3 / −1 high |
| 7–10 | the unmodified parents of #1–#4 | | — | | | | **controls** |

Pairwise edit distance between the six loops: min 3, median 6, max 6.

## The decision this list encodes

Four of CDR-H3's ten positions touch HER2. Measured from 1N8Z, heavy atom to
heavy atom, and the set is identical at 4.5 Å and 5.0 Å, so it is not an
artefact of the cutoff:

| | W99 | G100 | G101 | D102 | **G103** | **F104** | **Y105** | A106 | M107 | D108 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Å to HER2 | **4.05** | 6.99 | 5.89 | 5.08 | **2.42** | **4.27** | **3.32** | 6.83 | 9.36 | 9.28 |

Trastuzumab's own Trp-oxidation liability **is W99** — a contact. Removing all
three of its CDR-H3 liabilities therefore costs a residue that touches the
antigen. That is a real choice, and the evidence is split.

**Against touching the contacts:** the geometry, and the fact that no other
entry we looked at gives up W99.

**For touching them:** the 524,346 measured FACS variants, which say every one
of the four contact positions tolerates substitution, several preferring it:

| Position | Trastuzumab | Its high-fraction | Best alternative |
| --- | :-: | ---: | --- |
| W99 | W | 0.532 (n=177,081) | N 0.708, Y 0.656 (n=62,854) |
| G103 | G | 0.625 | **R 0.935** (n=42,561) |
| F104 | F | 0.614 | G 0.785 |
| Y105 | Y | 0.539 | W 0.674 |

A marginal from a binding screen is not an affinity measurement, and a crystal
contact is not an energy. Neither side wins on the evidence we have.

**So the list is split deliberately.** Four designs follow the tolerance
reading; two follow the contact reading. Both appear in the top four, which is
the only part of the ranking we control — every team's top four is tested
whatever happens.

The reason to split rather than pick is not timidity. Before this, all six of
our designs changed W99. That is not six bets — it is one bet written six
times, and if W99 turns out to be load-bearing, every design fails together
for the same reason. The split removes that single point of failure, and the
lab result then tells us which reading was right, which is worth having.

## What the three candidate sets actually measured

| | liability delta | contacts changed | FACS proxy rank | epitope (median, min) | ipTM |
| --- | --- | --- | --- | --- | --- |
| v1, sampled then inspected | −1 / −1 at best | 1–2 | — | 0.947, 0.85 | 0.9162 |
| v2, liability-clean | **−3 / −1** | 1–2, always W99 | **25–248** | **1.00**, 0.947 | 0.9118 |
| contact-preserving | −2 / −1 | **none** | 169–2,971 | 0.947, 0.90 | **0.9207** |

Freezing the contacts did **not** improve the predicted epitope and cost a
liability and a great deal of proxy score. Its only advantage is ipTM, by
0.009 — and our own calibration found that boltz2's ipTM does not rank the 182
measured variants, so that is close to no advantage at all. That is why the
contact reading gets two slots and not six.

## The stabilising mutation: VL Arg66Gly

Unchanged, on all six.

> Okazaki et al. 2023, *ACS Omega* 8(25):22945–22954, PMID 37396255.
> Trastuzumab scFv by DSF: **61.8 ± 0.2 °C → 66.7 ± 0.2 °C, ΔTm +4.9 °C**,
> binding unchanged.

S-tier needs all three judged properties improved, and four substitutions in
one loop give no reason to expect a higher melting temperature. R66G touches
only light-chain framework, leaves CDR-H3 and our novelty margin alone, and
adds no cysteine. It is one paper and one mutation, which is why the top four
are also submitted unmodified: the lab can then say what R66G did.

## What we are not claiming

- **Not that any of these binds.** Every number here is a prediction or a
  sequence property. The lab results come after the event.
- **Not that removed liability motifs are thermostability.** Isomerisation and
  oxidation are chemical degradation over time; the lab measures a melting
  temperature. Different things, and we do not conflate them.
- **Not that the proxy ranks binding.** It has no demonstrated skill inside a
  fixed mutation count, and every design here sits at exactly four.
- **The novelty margin is zero.** Three edits against a requirement of three,
  for every design, measured by edit distance against both the provided data
  and the official checker's full 610,109-sequence reference.

