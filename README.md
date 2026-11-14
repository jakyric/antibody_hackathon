# krishagen

A small, tested toolkit for the HER2 / trastuzumab CDRH3 problem: an entry-rule gate, an honest benchmark,
a pure-numpy binder model and a two-stage ranking. Plain Python with numpy and pytest. No GPU and no
framework, so anyone can run it.

This is one person's branch of the shared repository. It holds methods and tests only: no data, no model
outputs and no candidate lists.

The measured numbers behind the write-up are in `results/`, so any claim in it
can be recomputed rather than believed: `folding/calibration.csv` (does a score
order variants whose affinity is known), `folding/agreement.csv` (do the models
agree with each other), `folding/vs_trastuzumab.csv` (how each design compares
with the clinical baseline), `folding/scores.csv` (every score, every model) and
`epitope/epitope.csv` (contact overlap). Sequence columns are stripped by
`scripts/share_results.py` before anything is copied here, and a test asserts
that none survive.

**Start with [FINDINGS.md](FINDINGS.md).** It is the single record of what we measured and what it means,
including the things that did not work. The code here is what produced those numbers; the write-up is where
they are interpreted, with their intervals and their caveats. Where the two disagree, the write-up is the
one that has been checked.

## Layout

| Path | What it does |
| --- | --- |
| `src/abx/data.py` | Loads the FACS and SPR tables and makes the CDRH3 conventions explicit (the dataset's 10-mer excludes the closing Y) |
| `src/abx/gate.py` | The hard filters, run locally: format, length, unpaired cysteines, identity to trastuzumab, distance to the provided data. Whatever it cannot check is reported as `unknown`, never as a pass |
| `src/abx/benchmark.py` | Does a scorer rank known answers correctly? Reports both all rows and the legal region only, against random and edit-distance baselines |
| `src/abx/splits.py` | Splits of the FACS library that do not leak through near-neighbours |
| `src/abx/model.py` | A small binder model trained on FACS and tested on SPR, in numpy |
| `src/abx/generate.py`, `rank.py`, `pipeline.py` | Proposals, then gate, then filtering, then a two-stage ranking |
| `src/abx/metrics.py` | Spearman, AUROC, top-k precision and bootstrap intervals, without scipy |
| `src/abx/thermo.py` | Looks for places where a new VH-VL disulfide is geometrically possible, measured in the real structure. Builds the C-beta for glycine, which has none, so glycine positions are not silently skipped |
| `scripts/fold_panel.py`, `fold_queries.py` | Builds the evaluation panel: variants whose affinity was measured, spread across the range, plus the clinical baseline, and the identical folding input for every model |
| `scripts/fold_analyse.py` | Asks whether a folding model's confidence orders variants whose affinity is known, against the same baselines on the same variants |
| `scripts/epitope.py` | Contact epitope from a predicted complex, compared with the crystal structure and between designs. Geometry only, no confidence score |
| `scripts/disulfide_scan.py` | The disulfide scan, with the literature-backed pair marked |
| `tools/record_run.py` | Runs a command and records the command, commit, seed and input checksums in a run folder and a ledger |
| `tools/scan_repo.py` | Scans a repository for keys, credentials files and files over 50 MB before it is published |
| `tests/` | Tests for all of the above |

## Running it

```bash
make setup
make test
```

The tests and the benchmark read two public trastuzumab tables that are not included here: the FACS binding
classes (524,346 variants) and the SPR affinities (768 rows). They come from Chinery et al. (bioRxiv
10.1101/2024.03.26.586756). Put them in `data/` as `tz_her2_facs_3classes_524346.csv` and
`trastuzumab_her2_spr_all_768.csv`, then run `make verify`.

## Things worth knowing

- **Binding falls with distance from trastuzumab**, so a model can look good by learning edit distance alone.
  The benchmark therefore reports two numbers, all rows and the legal region, and says so when a model does not
  beat the edit-distance baseline there.
- **Splits must be checked for leakage.** The FACS library is a near-exhaustive mutagenesis, so a split on a few
  positions puts near-duplicates on both sides. See the notes at the top of `splits.py` for the numbers.
- **A class label is not an affinity.** Ranking therefore has two stages, each justified by a measured number
  (see `rank.py`).
- **State what you cannot check.** The gate returns `unknown` for anything it cannot verify locally.
- **Folding confidence is not a ranking function here.** Five models, the same panel, the same template and
  the same seed give five different answers to "how many of these beat the baseline", from none of them to
  most of them. Ask a model something with a known answer before ranking on it; `fold_analyse.py` and
  `epitope.py` are how we do that, and section 7 of the write-up has the numbers.
- **Check a filter against the thing it is supposed to catch.** Two of our own filters were wrong in the
  same way: one required an atom that glycine does not have, the other counted cysteines over the whole
  molecule when the rule applies per domain. Both returned clean answers to the wrong question.
