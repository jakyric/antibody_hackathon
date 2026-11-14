#!/usr/bin/env bash
# Same panel, same model, a different seed. The only thing this answers is how
# much two runs of the same model scatter — that is, whether a difference between
# two designs is larger than the noise at all.
#
# Must live in a file, not be passed in as arguments. If `$(...)` goes through a
# layer that quotes the arguments, it arrives as text: then you end up creating a
# directory that is literally named "$(realpath ...)", and the seed disappears
# without anything failing. That happened, and it is why this file exists.
#
#   bash scripts/fold_replicate.sh [seed]
set -euo pipefail

# Run recorder: tools/ in the shared repo, the skills folder in our own working copy.
RECORD=tools/record_run.py
[ -f "$RECORD" ] || RECORD=.claude/skills/repro-run/scripts/record_run.py

SEED="${1:-77}"

exec python "$RECORD" \
  --name "fold-boltz2-s$SEED" --seed "$SEED" \
  --input runs/panel/queries_all.json --input runs/panel/template.cif \
  --input runs/panel/panel.csv \
  --note "Boltz-2, same panel, seed $SEED: gives the noise level between runs" \
  -- bash -c 'bash scripts/fold_run.sh boltz2 runs/panel/queries_all.json \
      "$(realpath --relative-to="$PWD" "$RUN_DIR")/predictions" "$RUN_SEED"'
