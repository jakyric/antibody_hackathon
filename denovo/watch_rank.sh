#!/usr/bin/env bash
# Leaderboard + Pareto + status refresh loop. Stops at 16:00 local time or when denovo/STOP exists.
#   nohup bash denovo/watch_rank.sh > results/denovo/watch_rank.log 2>&1 &
cd "$(dirname "$0")/.." || exit 1
PY="denovo/.venv/bin/python -W ignore"
while [ ! -e denovo/STOP ] && [ "$(date +%H%M)" -lt "${PIPELINE_RANK_UNTIL_HHMM:-1600}" ]; do
  echo "=== $(date +%H:%M:%S)"
  $PY denovo/rank.py 2>&1 | tail -n 4
  $PY denovo/pareto.py 2>&1 | tail -n 4
  $PY denovo/status.py > /dev/null 2>&1
  sleep 300
done
echo "watch_rank stopped at $(date +%H:%M:%S)"
