#!/usr/bin/env bash
# Epitope overlap for the loops that are actually in the submission.
#
# `epitope.py` scans every fold run under runs/, which means re-reading every
# structure from every earlier panel. This points it at the submission runs
# only, so the number is about the designs we are submitting and nothing else.
#
#   bash scripts/submission_epitope.sh
set -euo pipefail

ONLY=runs/_submission_only
rm -rf "$ONLY"; mkdir -p "$ONLY"
n=0
for d in runs/*-fold-subpanel-*; do
  [ -d "$d/predictions" ] || continue
  ln -s "$(realpath "$d")" "$ONLY/$(basename "$d")"
  n=$((n + 1))
done
echo "$n submission fold run(s) linked into $ONLY"
[ "$n" -gt 0 ] || { echo "nothing to analyse"; exit 1; }

python scripts/epitope.py \
  --runs "$ONLY" \
  --panel runs/subpanel/panel.csv \
  --out results/epitope_submission

