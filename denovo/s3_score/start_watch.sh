#!/usr/bin/env bash
# Start the S3 continuous scorer in the background (idempotent: refuses if one is already running).
# Needs JAPANFOLD_API_KEY in the environment. Log: denovo/s3_score/watch.log, pid: denovo/s3_score/watch.pid.
# Stops by itself at 15:40 local time or when denovo/STOP exists; stop now with stop_watch.sh.
cd "$(dirname "$0")/../.." || exit 1
if [ -f denovo/s3_score/watch.pid ] && kill -0 "$(cat denovo/s3_score/watch.pid)" 2>/dev/null; then
  echo "already running: pid $(cat denovo/s3_score/watch.pid)"; exit 0
fi
setsid nohup denovo/.venv/bin/python denovo/s3_score/run_s3.py watch --tier2 "${TIER2:-dde,esm2,of3s1,af2ig}" --workers 8 --interval 90 \
  >> denovo/s3_score/watch.log 2>&1 < /dev/null &
echo $! > denovo/s3_score/watch.pid
echo "started pid $(cat denovo/s3_score/watch.pid), log denovo/s3_score/watch.log"
