#!/bin/bash
# compact campaign status: designs, failures, spend, S3-scored designs, time
cd "$(dirname "$0")/../.."
python3 - <<'PY'
import csv, json, os, time
n=sum(1 for _ in open("results/denovo/designs_rfd4jf.csv"))-1 if os.path.exists("results/denovo/designs_rfd4jf.csv") else 0
spent=sum(json.loads(x)["cost"] for x in open("denovo/s1_design/spend_ledger.jsonl")) if os.path.exists("denovo/s1_design/spend_ledger.jsonl") else 0
fails=sum(1 for _ in open("denovo/s1_design/failed_jobs.jsonl")) if os.path.exists("denovo/s1_design/failed_jobs.jsonl") else 0
try: sc=sum(1 for r in csv.DictReader(open("results/denovo/score_of3.csv")) if "rfd4jf" in r["design_id"])
except Exception: sc=0
print(f"{time.strftime('%H:%M:%S')} designs={n} spend=${spent:.3f} failed_jobs={fails} s3_scored={sc}")
PY
tail -n 2 denovo/s1_design/wave_auto.log 2>/dev/null; tail -n 1 denovo/s1_design/wave1.log
