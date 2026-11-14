#!/usr/bin/env bash
# ON THE VM: ESMFold-2 only for designs that have no JapanFold ESMFold-2 result, best first; claims chunks with a mkdir lock.
cd "$(dirname "$0")" || exit 1
Q=quay.io/apheris/oslo-hackathon
RUN=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)
for q in esmchunks/e_*.json; do
  [ -e STOP ] && break
  n=$(basename "$q" .json); n=${n#e_}
  mkdir "out/elock_$n" 2>/dev/null || continue
  while [ "$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)" -lt 12000 ]; do [ -e STOP ] && break 2; sleep 10; done
  echo "$(date +%T) esm2-only $n"
  "${RUN[@]}" $Q:apheris-esmfold2-v1.4.0-release python -m apheris_esmfold2.cli predict --input "$q" --output out/esm2_e$n --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_esm2_e$n.txt 2>&1
done
echo "$(date +%T) esm-only worker finished"
