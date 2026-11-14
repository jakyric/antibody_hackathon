#!/usr/bin/env bash
# Start, inspect, stop and finish the whole de novo VHH pipeline (continuous, resumable loops).
#
#   bash denovo/run_pipeline.sh dry-run     # print the resolved plan and the exact commands; starts nothing, spends nothing
#   bash denovo/run_pipeline.sh start       # launch the loops in the background (logs and pids in denovo/.run/)
#   bash denovo/run_pipeline.sh status      # which loops are alive, funnel (status.py), last log lines, JapanFold balance
#   bash denovo/run_pipeline.sh stop        # graceful stop (writes denovo/STOP); in-flight server jobs finish on their own
#   bash denovo/run_pipeline.sh ptx         # optional: Protenix-v2 on Nebius for the current top designs (needs NEBIUS_CREDS)
#   bash denovo/run_pipeline.sh final       # stop, drain, rank + Pareto + status, copy small snapshots to results_denovo/, safety scan
#
# Settings (environment variables; the defaults are safe):
#   PIPELINE_UNTIL=HH:MM         everything is finished by then, local time (default: now + 2 h; needs GNU date)
#   PIPELINE_GEN_MARGIN=50       minutes before PIPELINE_UNTIL at which design generation stops launching new jobs
#   PIPELINE_BUDGET_DESIGN=12    USD cap for RFD4 design on JapanFold (about $0.003 per design, so this is generous)
#   PIPELINE_BUDGET_S3=20        USD cap for structure scoring on JapanFold
#   PIPELINE_TIER2=dde,esm2,of3s1,af2ig   second-tier scorers for the top third of designs (OpenDDE, ESMFold-2, 2nd OpenFold3 seed, AF2-IG)
#   PIPELINE_ALLOW_VM=0          1 = ALSO use the shared VM (monoyield on its GPU, VM abguard); needs PIPELINE_VM_HOST and
#                                PIPELINE_VM_KEY. Leave at 0 unless the VM's owner/teammates agreed: its GPU is shared.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=denovo/.venv/bin/python
RUN=denovo/.run; mkdir -p "$RUN"
CMD="${1:-}"

UNTIL="${PIPELINE_UNTIL:-$(date -d '+2 hours' +%H:%M)}"
hm() { date -d "@$(( $(date -d "$UNTIL" +%s) - $1 * 60 ))" +%H:%M; }   # UNTIL minus N minutes, e.g. hm 20
GEN_UNTIL="$(hm "${PIPELINE_GEN_MARGIN:-50}")"; S3_UNTIL="$(hm 20)"; S4_UNTIL="$(hm 15)"; S2_UNTIL="$(hm 15)"
BUDGET_DESIGN="${PIPELINE_BUDGET_DESIGN:-12}"; BUDGET_S3="${PIPELINE_BUDGET_S3:-20}"
TIER2="${PIPELINE_TIER2:-dde,esm2,of3s1,af2ig}"; ALLOW_VM="${PIPELINE_ALLOW_VM:-0}"
S4_ARGS=(--watch 240 --tools liab nb9 nanomelt tnp --where local --stop-at "$S4_UNTIL")
if [ "$ALLOW_VM" = 1 ]; then S4_ARGS=(--watch 240 --tools liab nb9 nanomelt tnp monoyield --where auto --allow-vm --dtype bf16 --stop-at "$S4_UNTIL"); fi

declare -A CMDS
CMDS[rules]="$PY denovo/s2_rules/rules.py --watch --interval 120 --until $S2_UNTIL"
CMDS[design]="$PY denovo/s1_design/run_design.py --auto --until $GEN_UNTIL --jobs-per-wave 24 --workers 40 --max-inflight 12 --budget $BUDGET_DESIGN"
CMDS[score]="PIPELINE_S3_UNTIL=$S3_UNTIL PIPELINE_BUDGET_S3=$BUDGET_S3 TIER2=$TIER2 bash denovo/s3_score/start_watch.sh   # background loop of its own"
CMDS[props]="$PY denovo/s4_props/props.py ${S4_ARGS[*]}"
CMDS[rank]="PIPELINE_RANK_UNTIL_HHMM=${UNTIL//:/} bash denovo/watch_rank.sh   # rank.py + pareto.py + status.py every 5 min"
ORDER=(rules design score props rank)

alive() { [ -f "$1" ] && kill -0 "$(cat "$1")" 2>/dev/null; }
need_env() {
  [ -x "$PY" ] || { echo "denovo/.venv missing: run  bash denovo/setup.sh" >&2; return 1; }
  [ -n "${JAPANFOLD_API_KEY:-}" ] || { echo "JAPANFOLD_API_KEY is not set (design and scoring need it)" >&2; return 1; }
  if [ "$ALLOW_VM" = 1 ] && { [ -z "${PIPELINE_VM_HOST:-}" ] || [ -z "${PIPELINE_VM_KEY:-}" ]; }; then
    echo "PIPELINE_ALLOW_VM=1 needs PIPELINE_VM_HOST and PIPELINE_VM_KEY" >&2; return 1; fi
}

case "$CMD" in
  dry-run)
    echo "plan (nothing is started, nothing is spent)"
    echo "  finish by $UNTIL | design stops launching at $GEN_UNTIL | scoring until $S3_UNTIL | props/rules until $S4_UNTIL"
    echo "  caps: design \$$BUDGET_DESIGN, scoring \$$BUDGET_S3 | tier 2: $TIER2 | VM: $([ "$ALLOW_VM" = 1 ] && echo ALLOWED || echo off)"
    echo "  keys: JAPANFOLD_API_KEY $([ -n "${JAPANFOLD_API_KEY:-}" ] && echo set || echo NOT set) | NEBIUS_CREDS $([ -n "${NEBIUS_CREDS:-}" ] && echo set || echo not set)"
    for k in "${ORDER[@]}"; do printf '  %-7s %s\n' "$k" "${CMDS[$k]}"; done
    if [ -x "$PY" ]; then echo; echo "campaign plan from run_design.py (dry run):"; "$PY" denovo/s1_design/run_design.py --auto --dry-run 2>&1 | head -12 || true
    else echo; echo "(denovo/.venv not created yet: bash denovo/setup.sh)"; fi
    ;;
  start)
    need_env || exit 1
    rm -f denovo/STOP
    for k in rules design props rank; do
      if alive "$RUN/$k.pid"; then echo "$k already running (pid $(cat "$RUN/$k.pid"))"; continue; fi
      # shellcheck disable=SC2086
      if [ "$k" = rank ]; then PIPELINE_RANK_UNTIL_HHMM="${UNTIL//:/}" setsid nohup bash denovo/watch_rank.sh >> "$RUN/$k.log" 2>&1 < /dev/null &
      else setsid nohup ${CMDS[$k]} >> "$RUN/$k.log" 2>&1 < /dev/null & fi
      echo $! > "$RUN/$k.pid"; echo "started $k (pid $!, log $RUN/$k.log)"
    done
    PIPELINE_S3_UNTIL="$S3_UNTIL" PIPELINE_BUDGET_S3="$BUDGET_S3" TIER2="$TIER2" bash denovo/s3_score/start_watch.sh
    echo "all loops started; follow with:  bash denovo/run_pipeline.sh status"
    ;;
  status)
    for k in rules design props rank; do alive "$RUN/$k.pid" && echo "$k: running (pid $(cat "$RUN/$k.pid"))" || echo "$k: not running"; done
    alive denovo/s3_score/watch.pid && echo "score: running (pid $(cat denovo/s3_score/watch.pid))" || echo "score: not running"
    [ -x "$PY" ] && "$PY" -W ignore denovo/status.py 2>/dev/null | sed -n 1,12p
    for l in "$RUN"/*.log denovo/s3_score/watch.log; do [ -f "$l" ] && { echo "--- $(basename "$l")"; tail -n 2 "$l" | cut -c1-200; }; done
    ;;
  stop)
    touch denovo/STOP
    echo "STOP written: loops finish their current cycle and exit; waiting up to 20 s"; sleep 20
    for p in "$RUN"/*.pid denovo/s3_score/watch.pid; do [ -f "$p" ] && alive "$p" && kill "$(cat "$p")" 2>/dev/null && echo "terminated $(basename "$p")"; done
    echo "stopped (jobs already submitted to JapanFold/Nebius keep running there; a restart adopts them, nothing is paid twice)"
    ;;
  ptx)
    [ -n "${NEBIUS_CREDS:-}" ] || { echo "NEBIUS_CREDS is not set" >&2; exit 1; }
    F="$HOME/.claude/skills/nebius-scientific-ai/forge.py"; [ -f "$F" ] || { echo "nebius-scientific-ai skill not installed" >&2; exit 1; }
    "$PY" denovo/s3d_protenix/build_input.py || exit 1
    OUT=$(python3 "$F" submit protenix-v2 --entry 'protenix-input|protenix-input-json/v1|nebius-runs/ptxtop-input.json|application/json' \
      --params '{"checkpoint":"protenix-v2","msa_mode":"none","sample_count":3,"model_seeds":[7]}' --label ptx-top --no-wait 2>&1); echo "$OUT" | tail -n 3
    OP=$(echo "$OUT" | grep -o 'op-[0-9a-f]\{16\}' | head -1); [ -n "$OP" ] || { echo "no run id found in the output" >&2; exit 1; }
    PTX_OP="$OP" setsid nohup bash denovo/s3d_protenix/finish.sh >> "$RUN/ptx.log" 2>&1 < /dev/null &
    echo "Protenix-v2 run $OP submitted; finish.sh fetches and parses it into results/denovo/score_ptx.csv (log $RUN/ptx.log)"
    ;;
  final)
    touch denovo/STOP
    echo "draining: waiting up to 3 min for the loops to exit"
    for i in $(seq 1 18); do
      n=0; for p in "$RUN"/*.pid denovo/s3_score/watch.pid; do [ -f "$p" ] && alive "$p" && n=$((n+1)); done
      [ "$n" = 0 ] && break; sleep 10; done
    "$PY" -W ignore denovo/rank.py | head -3
    "$PY" -W ignore denovo/pareto.py | head -14
    "$PY" -W ignore denovo/status.py > /dev/null
    "$PY" denovo/snapshot.py --src results/denovo --dst results_denovo
    python3 denovo/scan_public.py . || echo "SCAN FINDINGS above: fix them before committing or pushing"
    echo; echo "final list: results_denovo/submission_top20_pareto.csv  (full tables in results_denovo/, per-design reasons inside)"
    ;;
  *) sed -n 2,22p "$0"; exit 2;;
esac
