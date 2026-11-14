#!/usr/bin/env bash
# Fold the +/- R66G panel with Boltz-2 at several seeds.
#
# Matched pairs, same model, same template, one seed at a time. The repeated
# seeds are not redundancy: two predictions of the SAME sequence give the noise
# floor against which a pair difference has to be judged.
#
#   bash scripts/r66g_fold.sh [seed ...]
set -euo pipefail

SEEDS=("${@:-42 77 123}")
[ "$#" -gt 0 ] || SEEDS=(42 77 123)

python scripts/r66g_panel.py --out runs/r66g/panel.csv
cp runs/panel/template.cif runs/r66g/template.cif
python scripts/fold_queries.py \
  --panel runs/r66g/panel.csv \
  --reference runs/panel/reference_queries.json \
  --out runs/r66g/queries.json

for SEED in "${SEEDS[@]}"; do
  echo "=== boltz2 seed $SEED ==="
  python .claude/skills/repro-run/scripts/record_run.py \
    --name "r66g-boltz2-s$SEED" --seed "$SEED" \
    --input runs/r66g/queries.json --input runs/r66g/template.cif \
    --input runs/r66g/panel.csv \
    --note "R66G matched pairs: trastuzumab + 3 designs, with and without VL R66G, boltz2" \
    -- bash -c "bash scripts/fold_run.sh boltz2 runs/r66g/queries.json \"\$(realpath --relative-to=\"\$PWD\" \"\$RUN_DIR\")/predictions\" \"\$RUN_SEED\"" \
    || echo "!! seed $SEED failed, continuing"
done
echo "=== done ==="

