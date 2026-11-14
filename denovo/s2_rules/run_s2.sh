#!/usr/bin/env bash
# Start / inspect / stop the two S2 watchers.   bash denovo/s2_rules/run_s2.sh start|status|stop
cd "$(dirname "$0")/../.." || exit 1
PY=denovo/s2_rules/.venv/bin/python
LOG=denovo/s2_rules/watch.log
case "$1" in
  start)
    rm -f denovo/STOP
    nohup $PY denovo/s2_rules/rules.py --watch --interval 120 --until 15:45 >> $LOG 2>&1 &
    echo "rules watcher pid $!"
    # diversification is NOT RUN (framework stays intact, decision 14:15 CET); opt in with:  bash run_s2.sh start-fwdiv
    echo "(fwdiv watcher not started)"
    ;;
  start-fwdiv)   # documented fallback only
    nohup $PY denovo/s2_rules/diversify.py --watch --interval 240 --until 15:30 >> $LOG 2>&1 &
    echo "fwdiv watcher pid $!"
    ;;
  status)
    pgrep -af "s2_rules/(rules|diversify).py --watch" || echo "no S2 watcher running"
    tail -n 8 $LOG
    ;;
  stop)   # only the S2 watchers
    pkill -f "s2_rules/(rules|diversify).py --watch" 2>/dev/null; echo stopped
    ;;
  stop-all)   # the global STOP file: every stage that honours it exits (S2 watchers within ~2 s)
    touch denovo/STOP; sleep 3; pkill -f "s2_rules/(rules|diversify).py --watch" 2>/dev/null; echo "STOP file written"
    ;;
  *) echo "usage: $0 start|status|stop|stop-all"; exit 2;;
esac
