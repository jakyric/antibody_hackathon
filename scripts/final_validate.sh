#!/usr/bin/env bash
# The organisers' own validator, run on the exact committed submission file.
#
# Two runs, because a gate that passes everything is indistinguishable from a
# gate that is switched off:
#   1. the ten submitted rows      -> must exit 0
#   2. the same ten plus trastuzumab -> must exit 1, naming trastuzumab
#
# Exit codes are captured from the official checker itself, not from a pipe. ${PIPESTATUS[0]}
# matters here: piping into `head` would report head's status and the log would
# claim a pass that never happened.
set -uo pipefail

GATE="${GATE:-.}"
ABG="${ABG:?set ABG to the official checker binary}"
cd "$GATE" || exit 2

echo "the official checker team-validate on the committed results/SUBMISSION.csv"
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "the official checker: $ABG"
echo "md5 of the validated file: $(md5sum SUBMISSION_committed.csv | cut -d' ' -f1)"
echo "reference: provided.duckdb, built from the provided HER2 data"
echo

echo "########## 1. THE TEN ALONE -- every row must pass ##########"
"$ABG" team-validate SUBMISSION_committed.csv --reference-db provided.duckdb 2>&1 | head -4
A=${PIPESTATUS[0]}
echo ">>> EXIT: $A   (0 = every row compliant)"
echo

echo "########## 2. PLUS TRASTUZUMAB -- the control MUST fail ##########"
"$ABG" team-validate SUB_with_tz.csv --reference-db provided.duckdb 2>&1 | head -4
B=${PIPESTATUS[0]}
echo ">>> EXIT: $B   (1 is CORRECT here: trastuzumab must be rejected)"
echo

if [ "$A" -eq 0 ] && [ "$B" -eq 1 ]; then
  echo "RESULT: gate behaves correctly -- ten pass, the known-bad control is caught."
  exit 0
fi
echo "RESULT: UNEXPECTED (ten=$A, control=$B). Do not submit on this result."
exit 1

