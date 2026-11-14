#!/usr/bin/env bash
# Can a judge actually run this repository? Tests the proposed `make setup` fix.
#
# The published Makefile does `python3 -m venv .venv` and then `pip install`.
# On a stock Debian/Ubuntu python that fails: the distribution splits ensurepip
# into a separate package, so the venv is created WITHOUT pip and every later
# step dies with "No module named pip". The judge sees a repository that does
# not build.
#
# This script runs the fixed sequence against the clean clone and reports
# whether the tests and the benchmark then pass.
#
#   bash scripts/check_clean_setup.sh [path-to-clean-clone]
set -uo pipefail

REPO="${1:-.}"
cd "$REPO" || { echo "no clone at $REPO"; exit 2; }
echo "repo: $REPO"
echo "commit: $(git log --oneline -1)"

rm -rf .venv
python3 -m venv .venv 2>/dev/null || python3 -m venv --without-pip .venv
if ! .venv/bin/python -m pip --version >/dev/null 2>&1; then
  echo "note: this python has no ensurepip (Debian/Ubuntu splits it out); bootstrapping pip"
  curl -sS https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py
  .venv/bin/python /tmp/get-pip.py -q
fi
.venv/bin/python -m pip -q install --upgrade pip
.venv/bin/python -m pip -q install -r requirements.txt
echo "pip: $(.venv/bin/python -m pip --version)"

echo
echo "=== make test ==="
make test 2>&1 | tail -12
echo ">>> TEST EXIT: ${PIPESTATUS[0]}"

echo
echo "=== make verify ==="
make verify 2>&1 | tail -50
echo ">>> VERIFY EXIT: ${PIPESTATUS[0]}"
