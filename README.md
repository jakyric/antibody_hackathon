# krishagen

A small, tested toolkit for the HER2 / trastuzumab CDRH3 problem: an entry-rule gate, an honest benchmark,
a pure-numpy binder model and a two-stage ranking. Plain Python with numpy and pytest. No GPU and no
framework, so anyone can run it.

This is one person's branch of the shared repository. It holds methods and tests only: no data, no model
outputs and no candidate lists.

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
