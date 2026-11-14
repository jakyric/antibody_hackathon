#!/usr/bin/env bash
# Fold the full panel with several fold models in sequence, on a single GPU.
# Each run goes through record_run.py, so it is listed in results/ledger.csv.
#
#   bash scripts/fold_sweep.sh <seed> <model> [model ...]
set -euo pipefail

# Run recorder: tools/ in the shared repo, the skills folder in our own working copy.
RECORD=tools/record_run.py
[ -f "$RECORD" ] || RECORD=.claude/skills/repro-run/scripts/record_run.py

SEED="${1:?seed}"; shift
QUERIES=runs/panel/queries_all.json

for MODEL in "$@"; do
  echo "=== $MODEL seed $SEED ==="
  python "$RECORD" \
    --name "fold-$MODEL-s$SEED" --seed "$SEED" \
    --input "$QUERIES" --input runs/panel/template.cif --input runs/panel/panel.csv \
    --note "fold model $MODEL, full panel (25 calibration + 20 design), 1 diffusion sample" \
    -- bash -c "bash scripts/fold_run.sh $MODEL $QUERIES \"\$(realpath --relative-to=\"\$PWD\" \"\$RUN_DIR\")/predictions\" \"\$RUN_SEED\"" \
    || echo "!! $MODEL failed, continuing"
done
echo "=== done ==="
