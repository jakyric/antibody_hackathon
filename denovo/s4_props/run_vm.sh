#!/usr/bin/env bash
# Run NanoMelt + TNP + liabilities + monoyield ON THE SHARED VM for every sequence in one CSV, and copy results back.
#
#   bash denovo/s4_props/run_vm.sh in.csv out_dir [tools] [dtype] [ncores]
#     in.csv  : design_id,vhh_seq[,cdr1,cdr2,cdr3]
#     out_dir : local folder that receives nanomelt.csv tnp.csv liab.csv my_monoyield.csv my_mdratio.csv status.txt + logs
#     tools   : comma list, default nanomelt,tnp,liab,nb9,monoyield
#     dtype   : auto|fp32|bf16 for monoyield (auto = fp32 only if >= 45 GB GPU memory is free, else bf16)
#     ncores  : TNP pool size on the VM (default 6)
# First time on a fresh VM:  bash denovo/s4_props/vm_install.sh push && bash denovo/s4_props/vm_install.sh install
# Needs only ssh/scp/rsync and the ssh key named by PIPELINE_VM_KEY (host in PIPELINE_VM_HOST).  Touches nothing outside ~/moritz/s4 on the VM.
set -euo pipefail
# User rule (14 Nov 14:33): nothing may start on the shared VM (GPU or otherwise) without the user's explicit permission.
[ "${S4_ALLOW_VM:-}" = 1 ] || { echo "run_vm.sh blocked: VM use needs the user's explicit permission (then: S4_ALLOW_VM=1 bash run_vm.sh ...)" >&2; exit 5; }
IN="$1"; OUT="$2"; TOOLS="${3:-nanomelt,tnp,liab,nb9,monoyield}"; DT="${4:-auto}"; NC="${5:-6}"
VM="${PIPELINE_VM_HOST:?set PIPELINE_VM_HOST=user@host to use the VM}"; KEY="${PIPELINE_VM_KEY:?set PIPELINE_VM_KEY=path of the ssh key}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CM="-o ControlMaster=auto -o ControlPath=$HOME/.ssh/cm-s4-%C -o ControlPersist=120"   # one TCP connection for all calls (sshd throttles bursts)
SSH=(ssh -i "$KEY" -o BatchMode=yes -o ServerAliveInterval=30 $CM "$VM")
retry() { local n=0 rc; until "$@"; do rc=$?; [ $rc -ne 255 ] && return $rc; n=$((n+1)); [ $n -ge 10 ] && return 255; sleep 6; done; }  # sshd on the shared VM sometimes resets new connections
mkdir -p "$OUT"
retry "${SSH[@]}" true || { echo "VM unreachable over ssh (connection resets); try again in a minute" >&2; exit 4; }
# User permission (14 Nov 14:40): monoyield on the VM GPU is allowed as long as it does not disturb other people's jobs.
# Guard: monoyield needs >= 28 GB of FREE GPU memory (bf16 ESMC-6B needs ~13 GB; leave generous headroom); otherwise skip
# this cycle (exit 6) and retry at the next interval. Never touches anyone else's processes.
if [[ ",$TOOLS," == *",monoyield,"* ]]; then
  FREE=$(retry "${SSH[@]}" "nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1" 2>/dev/null || echo 0)
  if [ "${FREE:-0}" -lt 28000 ]; then
    echo "monoyield skipped: only ${FREE:-0} MiB GPU memory free on the shared VM (others are using it); retry next cycle" >&2; exit 6
  fi
fi
if ! retry "${SSH[@]}" 'test -f ~/moritz/s4/INSTALL_OK'; then
  echo "VM env missing: run  bash $HERE/vm_install.sh push && bash $HERE/vm_install.sh install  first" >&2; exit 3
fi
retry bash "$HERE/vm_install.sh" push >/dev/null            # keep the VM copies of the scripts current (cheap rsync)
RUN="run-$(date -u +%Y%m%dT%H%M%S)-$$"
retry "${SSH[@]}" "mkdir -p ~/moritz/s4/runs/$RUN"
retry scp -q -i "$KEY" $CM "$IN" "$VM:moritz/s4/runs/$RUN/in.csv"
"${SSH[@]}" "bash ~/moritz/s4/denovo/s4_props/remote_s4.sh ~/moritz/s4/runs/$RUN '$TOOLS' '$DT' '$NC'" || true
for f in status.txt nanomelt.csv tnp.csv liab.csv nb9.csv my_monoyield.csv my_mdratio.csv liab.log tnp.log nanomelt.log my_monoyield.log my_mdratio.log; do
  retry scp -q -i "$KEY" $CM "$VM:moritz/s4/runs/$RUN/$f" "$OUT/$f" 2>/dev/null || true
done
cat "$OUT/status.txt"
if grep -q FAILED "$OUT/status.txt"; then echo "some tools FAILED, see $OUT/*.log" >&2; exit 1; fi
