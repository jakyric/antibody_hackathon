#!/usr/bin/env bash
# Stop the S3 watch now (in-flight JapanFold jobs keep running server-side; a restart adopts and parses them).
cd "$(dirname "$0")/../.." || exit 1
if [ -f denovo/s3_score/watch.pid ] && kill -0 "$(cat denovo/s3_score/watch.pid)" 2>/dev/null; then
  kill "$(cat denovo/s3_score/watch.pid)" && echo "stopped pid $(cat denovo/s3_score/watch.pid)"
else echo "not running"; fi
