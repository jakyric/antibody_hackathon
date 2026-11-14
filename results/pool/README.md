# The liability-clean pool: how the new designs were found

Every file here was produced by the scripts listed below. The submission built from it is
`../SUBMISSION.csv`; the reasoning is in `../SUBMISSION.md`.

## Why a new search at all

Scoring a published comparison set through our own pipeline
(`../external/README.md`) showed the one place we were behind: their pool
reduced trastuzumab's high-severity liability count in nine designs of ten,
ours in four of twenty, and nine of ours introduced `FragH:DP`.

Our twenty were drawn from a position-weight matrix and checked for
liabilities afterwards. Nothing in the draw knew what a liability was, so the
clean loops were never enumerated — not rejected, never looked at.

## The search

| Step | Script | Result |
| --- | --- | --- |
| enumerate the clean space | `src/abx/pool.py` | **691,064** loops at exactly 4 substitutions, carrying none of trastuzumab's CDR-H3 motifs and adding none |
| score and cut | `scripts/clean_pool.py` | top 3,000 by the FACS model |
| novelty, exhaustive | same | **855** at ≥3 edits from all 525,063 provided loops — every one at exactly 3 |
| re-filter on the widened motif list | `scripts/clean_select.py` | 499 |
| spread | same | 8 picks, farthest-point from the proxy's best |
| fold + epitope | `scripts/clean_fold.sh`, `scripts/epitope.py` | 2 dropped below the epitope floor |
| rank | `scripts/clean_rank.py` | `ranked.csv`, the final six |
| build + validate | `scripts/build_final_clean.py`, the official checker | `SUBMISSION_v2.csv`, exit 0 |

The enumeration is complete, not sampled: no seed, and the same list every
time. That is why it is worth doing — the clean loops either exist or they do
not, and now we know.

## The motif list was wrong once, and the official checker caught it

`FORBIDDEN_PAIRS` in `src/abx/pool.py` started as the obvious motifs. A probe
of eight enumerated loops through the official checker's all-liabilities mode showed three of
them still carrying an `Isom` motif in CDR-H3: the tool reads **`DT`** as
isomerisation, which our list did not. The list is now wider than what
the official checker flags rather than narrower, and the finalists are still validated by
the official checker itself.

A motif list we wrote is a way to narrow a search. It is never a way to
declare a design clean.

## Controls that ran, and what they said

| Control | Result |
| --- | --- |
| trastuzumab through the novelty check | distance **0** — it collides with the provided data, as it must |
| trastuzumab through the official checker's team validation | **rejected**, `insufficient_cdrh3_novelty`; a gate it passed would be broken |
| trastuzumab folded in the same run | `iptm_ab_ag` **0.9405**, identical to the two earlier runs to four decimals |
| trastuzumab's predicted epitope vs its 1N8Z crystal epitope | Jaccard **1.00** — the floor under every epitope number here |

## Files

| File | What it is |
| --- | --- |
| `clean_pool.csv` | the 855 novel clean loops, with proxy score and nearest provided sequence |
| `selected.csv`, `selected8.csv` | the farthest-point picks (6 and 8) |
| `epitope/epitope.csv` | predicted HER2 contact overlap, all eight plus the parent |
| `fold_compare.csv` | ipTM and pLDDT for the new loops and the old twenty, one scale |
| `ranked.csv` | the final six, with the rank aggregation |
| `SUBMISSION_v2.csv` | the ten submitted rows |
| `validation_v2.json`, `liability_delta_v2.csv` | the official checker's verdict and the liability deltas |
| `novelty_full_v2.csv` | every submitted loop against all 525,063, by edit distance |
| `provided_ref_v2.csv` | the reference the official checker run was built from |

The two `.duckdb` reference databases are rebuilt from `provided_ref_v2.csv`
by the official checker's reference ingest and build steps and
are not kept here.

