#!/usr/bin/env bash
# Idempotent installer: NanoMelt + TNP (and their shared deps) in ONE micromamba env on the shared VM.
# Everything lives under ~/moritz/s4 (no sudo, nothing outside ~/moritz is touched).
# Usage (from the repo root, local machine):   bash denovo/s4_props/vm_install.sh push     # rsync sources
#                                              bash denovo/s4_props/vm_install.sh install  # remote, nohup, log in ~/moritz/s4/logs/install.log
#                                              bash denovo/s4_props/vm_install.sh status
set -euo pipefail
VM="${PIPELINE_VM_HOST:?set PIPELINE_VM_HOST=user@host to use the VM}"
KEY="${PIPELINE_VM_KEY:?set PIPELINE_VM_KEY=path of the ssh key}"
CM="-o ControlMaster=auto -o ControlPath=$HOME/.ssh/cm-s4-%C -o ControlPersist=120"
SSH="ssh -i $KEY -o BatchMode=yes $CM $VM"
retry() { local n=0 rc; until "$@"; do rc=$?; [ $rc -ne 255 ] && return $rc; n=$((n+1)); [ $n -ge 10 ] && return 255; sleep 6; done; }
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

case "${1:-}" in
push)
  retry $SSH 'mkdir -p ~/moritz/s4/{bin,logs,torch_cache,hf_cache,tools/TNP,tools/nanomelt,denovo/s4_props,results/denovo}'
  # same relative layout as the repo, so tnp_batch.py finds tools/TNP/src on both machines
  retry rsync -az -e "ssh -i $KEY -o BatchMode=yes $CM" "$ROOT/tools/TNP/src/" "$VM:moritz/s4/tools/TNP/src/"
  retry rsync -az -e "ssh -i $KEY -o BatchMode=yes $CM" --exclude notebooks --exclude '__pycache__' \
    "$ROOT/tools/nanomelt/nanomelt" "$ROOT/tools/nanomelt/setup.py" "$ROOT/tools/nanomelt/environment.yml" \
    "$ROOT/tools/nanomelt/LICENSE" "$VM:moritz/s4/tools/nanomelt/"
  retry rsync -az -e "ssh -i $KEY -o BatchMode=yes $CM" "$ROOT/denovo/s4_props/vm_install.sh" "$ROOT/denovo/s4_props/tnp_batch.py" \
    "$ROOT/denovo/s4_props/liab.py" "$ROOT/denovo/s4_props/nanobody_rules.py" "$ROOT/denovo/s4_props/remote_s4.sh" "$VM:moritz/s4/denovo/s4_props/" 2>/dev/null || true
  echo pushed
  ;;
install)
  $SSH 'cd ~/moritz/s4 && nohup bash denovo/s4_props/vm_install.sh _remote > logs/install.log 2>&1 & echo started'
  ;;
status)
  $SSH 'tail -n 15 ~/moritz/s4/logs/install.log; ls ~/moritz/s4/INSTALL_OK 2>&1'
  ;;
_remote)
  S="$HOME/moritz/s4"; cd "$S"
  export MAMBA_ROOT_PREFIX="$S/mamba"
  MM="$S/bin/micromamba"
  if [ ! -x "$MM" ]; then
    # official single static binary from the conda-forge/mamba project
    curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$S" bin/micromamba
  fi
  if [ ! -d "$MAMBA_ROOT_PREFIX/envs/s4" ]; then
    # versions mirror the local envs `tnp` and `nanomelt` (see tools/README.md); anarci/hmmer must come from bioconda
    "$MM" create -y -n s4 --override-channels -c conda-forge -c bioconda \
      python=3.12 pip anarci=2024.05.21 hmmer=3.4 biopython=1.88 openmm=8.6.1 pdbfixer=1.12 \
      freesasa=2.2.1 freesasa-python dssp 'numpy=2.5.*' 'scipy=1.18.*' pandas
  fi
  PIP="$MM run -n s4 python -m pip"
  $PIP install --no-input torch
  # --no-deps for ImmuneBuilder: its metadata would otherwise pull the unofficial PyPI `anarci` over bioconda's
  $PIP install --no-input --no-deps ImmuneBuilder==1.2
  $PIP install --no-input einops 'scikit-learn==1.9.1' 'xgboost==3.4.1' 'lightgbm==4.7.0' 'joblib==1.6.0' \
    'transformers==5.23.0' 'ablang==0.3.1' 'antiberty==1.0.0' 'fair-esm==2.0.0' 'kneed==0.9.0' 'kmedoids==0.5.5' matplotlib
  "$MM" run -n s4 python - <<'PY'
import torch, anarci, ImmuneBuilder, sklearn
print("torch", torch.__version__, "cuda", torch.cuda.is_available(), "sklearn", sklearn.__version__)
PY
  touch "$S/INSTALL_OK"
  echo INSTALL_DONE
  ;;
*)
  echo "usage: $0 push|install|status"; exit 2;;
esac
