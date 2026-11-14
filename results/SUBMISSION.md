# The final ten: what we submit and why

`results/SUBMISSION.csv`. Built by `scripts/build_final_clean.py` from
`results/pool/ranked.csv`, validated with the official checker's `team-validate
--reference-db` against a reference built from the provided HER2 data.
**The ten alone exit 0. With trastuzumab added as a control, 10 of 11 pass
and the control is the one that fails.**

The previous list is kept as `results/SUBMISSION_v1.md` and
`results/SUBMISSION_v1.csv`. What changed, and why, is at the bottom.

## The list

| # | ID | CDRH3 | Stab | Epitope | Liabilities vs trastuzumab |
| ---: | --- | --- | :-: | ---: | --- |
| 1 | `tz-clean-02-stab` | RGGEGFYANE | ✓ | **1.00** | −3 total, −1 high |
| 2 | `tz-clean-01-stab` | YGGDYFYGKD | ✓ | 0.95 | −3 total, −1 high |
| 3 | `tz-clean-03-stab` | FGGGGNYAKD | ✓ | **1.00** | −3 total, −1 high |
| 4 | `tz-clean-07-stab` | HGGDRFYKRD | ✓ | 0.95 | −3 total, −1 high |
| 5 | `tz-clean-06-stab` | YGGAGFHAVD | ✓ | **1.00** | −3 total, −1 high |
| 6 | `tz-clean-04-stab` | HTGDEFYAID | ✓ | 0.947 | −3 total, −1 high |
| 7 | `tz-clean-02` | RGGEGFYANE | — | | **control** for #1 |
| 8 | `tz-clean-01` | YGGDYFYGKD | — | | **control** for #2 |
| 9 | `tz-clean-03` | FGGGGNYAKD | — | | **control** for #3 |
| 10 | `tz-clean-07` | HGGDRFYKRD | — | | **control** for #4 |

Six distinct loops, pairwise edit distance **min 4, median 5, max 6** — six
different bets, not one written six times. Every design is exactly 4
substitutions from trastuzumab's CDR-H3 and exactly 3 edits from the nearest
of the 525,063 provided loops.

## The three things that determined this list

**1. Every design removes all three of trastuzumab's own CDR-H3 liabilities,
and adds none.**

Trastuzumab's loop `WGGDGFYAMD` carries `Isom:DG` (high severity), `MetOx:M`
and `TrpOx:W`. Removing all three costs three substitutions; the novelty rule
demands four; so exactly one position is free. That space is small enough to
enumerate completely — 691,064 loops — rather than sample.

The official checker's own liability report puts every one of the ten at **5 liabilities
(2 high) against trastuzumab's 8 (3 high)**. That is the floor reachable by
changing CDR-H3 alone, and all ten reach it. Expression failure means no
measurement on any axis, so this is the one sequence property worth
optimising before the lab sees anything.

**2. Removing those three residues is what the measured data wanted anyway.**

This is the part that could have gone badly: trastuzumab's CDR-H3 tryptophan
is a paratope residue, and trading it for a developability motif could cost
binding. The 524,346 FACS-measured variants say otherwise. Fraction scoring
`high`, per residue at each position:

| Position | Trastuzumab | Its score | Better options in the data |
| --- | :-: | ---: | --- |
| 1 | W | 0.532 (n=177,081) | Y 0.656 (n=62,854), N 0.708 (n=917) |
| 4 | D | 0.363 (n=30,934) | G 0.827, P 0.702, A 0.696 |
| 5 | G | 0.625 (n=123,584) | R 0.935 (n=42,561) |
| 9 | M | 0.499 (n=4,737) | V 0.608, L 0.584 (n=72,039) |

All four motif-bearing positions are ones the data says are improvable. The
liability-clean constraint and the binding evidence point the same way. That
is luck rather than design, but it is measured luck, and it is why this swap
is not a trade.

*Caveat: these are marginal per-position fractions from a FACS assay, which
says whether a variant binds, not how tightly.*

**3. Epitope decided the order, because nothing else could.**

All ten carry an identical liability profile, so developability — which
ordered the old list — is now a constant and orders nothing. The remaining
axes are the predicted HER2 contact set against trastuzumab's, and the FACS
proxy. Combined by rank aggregation (`scripts/clean_rank.py`), ties averaged.

The epitope floor holds: trastuzumab reproduces its own 1N8Z crystal epitope
at Jaccard **1.00** in this run, so a design's overlap means something. Two
candidates were dropped outright for drifting off it (0.75 and 0.90) rather
than ranked low — an optimisation-track design that no longer competes for
trastuzumab's site has lost that axis, not scored badly on it.

## The stabilising mutation: VL Arg66Gly

Unchanged from the previous list, and the argument is unchanged.

> Okazaki et al. 2023, *ACS Omega* 8(25):22945–22954, PMID 37396255.
> Trastuzumab scFv by DSF: **61.8 ± 0.2 °C → 66.7 ± 0.2 °C, ΔTm +4.9 °C**,
> binding unchanged.

S-tier needs all three judged properties improved. Four substitutions in one
CDR loop give no reason to expect a higher melting temperature, so an
unmodified design tops out at A-tier. R66G touches only light-chain
framework, leaves CDR-H3 and therefore our novelty margin untouched, and adds
no cysteine. It is still one paper and one mutation, which is why four of the
ten are the unmodified parents of the top four: the lab result can then say
what R66G did.

## What we are not claiming

- **Not that these bind better than the old list.** On the FACS proxy they
  score higher (9.2–11.6 against 6.3–9.4), but that proxy has no demonstrated
  skill inside a fixed mutation count, and every design here sits at exactly
  four. Read it as reassurance that we did not trade binding away for
  developability, not as evidence that we bought any.
- **Not that any of them binds at all.** Every number here is a prediction or
  a sequence property. "Expected", never "binds". The lab results come after
  the event.
- **The novelty margin is still zero.** Three edits against a requirement of
  three, for every design, as before — now measured by edit distance rather
  than by position, which is the metric the official checker actually uses.

## What changed from the previous list, and why

A comparison set of published designs, scored through our own pipeline
(`results/external/README.md`), showed one real gap and it was not in our
favour: nine of ten of theirs reduced the high-severity liability count,
against four of our twenty, and nine of our twenty introduced `FragH:DP`, an
acid-labile motif trastuzumab does not carry and none of theirs did.

That was a search problem, not a scoring problem. Our twenty were sampled
from a position-weight matrix and inspected for liabilities afterwards;
nothing in the draw knew what a liability was. `src/abx/pool.py` enumerates
the clean space instead. The new list is better than the old on every axis we
can measure and worse on none:

| | Old ten | New ten |
| --- | --- | --- |
| liability delta vs trastuzumab | −1 total, −1 high at best; most of the pool 0 or worse | **−3 total, −1 high, all ten** |
| designs introducing `FragH:DP` | 9 of the 20 it was drawn from | **none** |
| epitope overlap with trastuzumab | median 0.947, min 0.85 | **median 1.00, min 0.947** |
| FACS proxy score | 6.3 – 9.4 | **9.2 – 11.6** |
| pairwise distance between loops | min 3 | **min 4** |
| novelty margin | 3 edits, zero margin | 3 edits, zero margin |

