#!/usr/bin/env bash
# Pull + parse the VM results, then refresh leaderboard/Pareto/status, every ~90 s until 15:55 local or denovo/STOP.
cd "$(dirname "$0")/../.." || exit 1
PY="baseline/.venv/bin/python -W ignore"
while [ ! -e denovo/STOP ] && [ "$(date +%H%M)" -lt 1555 ]; do
  $PY denovo/vm_templated/collect.py 2>&1 | tail -n 1
  $PY denovo/rank.py 2>&1 | head -n 1
  $PY denovo/pareto.py 2>&1 | grep -E "^pareto" | cut -c1-200
  $PY denovo/status.py > /dev/null 2>&1
  sleep 90
done
echo "collect loop stopped $(date +%T)"
