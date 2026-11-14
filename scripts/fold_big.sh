#!/usr/bin/env bash
# The big calibration panel: 100 SPR variants with measured KD, plus
# trastuzumab, in Boltz-2. The purpose is a single number: how certain it is
# that fold confidence is related to measured affinity at all.
#
# On 25 variants we got rho +0.41 with an interval that barely clears zero.
# With 100 the interval becomes roughly half as wide, and then we will know
# whether the number is real or just random chance.
#
# Separate file because `$(...)` does not survive running in the background.
set -euo pipefail

# Run recorder: tools/ in the shared repo, the skills folder in our own working copy.
RECORD=tools/record_run.py
[ -f "$RECORD" ] || RECORD=.claude/skills/repro-run/scripts/record_run.py

SEED="${1:-42}"

exec python "$RECORD" \
  --name "fold-boltz2-big-s$SEED" --seed "$SEED" \
  --input runs/panel-big/queries_all.json --input runs/panel-big/template.cif \
  --input runs/panel-big/panel.csv \
  --note "Boltz-2 on 100 SPR variants with measured KD, for a narrower interval" \
  -- bash -c 'bash scripts/fold_run.sh boltz2 runs/panel-big/queries_all.json \
      "$(realpath --relative-to="$PWD" "$RUN_DIR")/predictions" "$RUN_SEED"'

