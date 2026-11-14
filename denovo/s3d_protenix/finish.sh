#!/usr/bin/env bash
# Wait for the Nebius protenix-v2 run, fetch and parse it. Gives up (leaving the run alone) at 15:35 local time.
#   nohup bash denovo/s3d_protenix/finish.sh > denovo/s3d_protenix/finish.log 2>&1 &
# Resume by hand:  python3 ~/.claude/skills/nebius-scientific-ai/forge.py wait op-70cf0433a9fbb3f9
cd "$(dirname "$0")/../.." || exit 1
OP="${PTX_OP:?set PTX_OP to the Nebius run id printed by forge.py submit}"
F=~/.claude/skills/nebius-scientific-ai/forge.py
while [ "$(date +%H%M)" -lt "${PIPELINE_PTX_UNTIL_HHMM:-1535}" ] && [ ! -e denovo/STOP ]; do
  S=$(python3 "$F" status "$OP" 2>&1 | tail -n 1)
  echo "$(date +%H:%M:%S) $S"
  case "$S" in
    *"result ready: True"*) break ;;
    *failed*|*cancel*) echo "run ended without result"; exit 2 ;;
  esac
  sleep 60
done
case "$S" in *"result ready: True"*) ;; *) echo "not ready by deadline; resume with: python3 $F wait $OP"; exit 3 ;; esac
python3 "$F" fetch "$OP" --dest results/ptx-top && denovo/.venv/bin/python -W ignore denovo/s3d_protenix/parse_ptx.py
echo "done $(date +%H:%M:%S)"
