#!/usr/bin/env bash
# Fold the loops that are actually in the submission, to fill the epitope column.
#
# The panel is rebuilt from results/SUBMISSION.csv every time, so it can never
# describe an older generation of designs than the one we are submitting.
#
#   bash scripts/subpanel_fold.sh [seed ...]
set -euo pipefail

SEEDS=("$@")
[ "${#SEEDS[@]}" -gt 0 ] || SEEDS=(42)

python scripts/submission_panel.py --out runs/subpanel/panel.csv
cp runs/panel/template.cif runs/subpanel/template.cif
python scripts/fold_queries.py \
  --panel runs/subpanel/panel.csv \
  --reference runs/panel/reference_queries.json \
  --out runs/subpanel/queries.json

for SEED in "${SEEDS[@]}"; do
  echo "=== boltz2 seed $SEED ==="
  python .claude/skills/repro-run/scripts/record_run.py \
    --name "fold-subpanel-boltz2-s$SEED" --seed "$SEED" \
    --input runs/subpanel/queries.json --input runs/subpanel/panel.csv \
    --input results/SUBMISSION.csv \
    --note "Boltz-2 on the submitted loops, to fill the epitope column" \
    -- bash -c "bash scripts/fold_run.sh boltz2 runs/subpanel/queries.json \"\$(realpath --relative-to=\"\$PWD\" \"\$RUN_DIR\")/predictions\" \"\$RUN_SEED\"" \
    || echo "!! seed $SEED failed, continuing"
done
echo "=== done ==="
