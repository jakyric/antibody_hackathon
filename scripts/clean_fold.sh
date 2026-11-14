#!/usr/bin/env bash
# Fold the liability-clean picks together with trastuzumab, through the same
# configuration as our existing panel run: same template, same chain layout,
# same seed. The parent is in the run so the epitope axis has its baseline,
# and so the run can be checked against the number we already have for it.
#
#   bash scripts/clean_fold.sh <model> [seed]
set -euo pipefail

MODEL="${1:?model}"
SEED="${2:-42}"
OUT="runs/clean/$MODEL-s$SEED"

mkdir -p runs/clean
cp -f runs/panel/template.cif runs/clean/template.cif

python scripts/clean_panel.py
python scripts/fold_queries.py \
  --panel candidates/clean_panel.csv \
  --reference runs/panel/reference_queries.json \
  --out runs/clean/queries.json

rm -rf "$OUT"
bash scripts/fold_run.sh "$MODEL" runs/clean/queries.json "$OUT" "$SEED"
