#!/usr/bin/env bash
# Parallel worker (ON THE VM): claims the next unclaimed chunk (mkdir lock, atomic) in priority order and runs templated OpenFold3
# then untemplated ESMFold-2 on it. Several workers share the H100 (each OpenFold3 process needs about 10-15 GB); a worker waits
# whenever less than 20 GB of GPU memory is free, so it can never push other people's jobs out of memory. `touch STOP` stops all.
cd "$(dirname "$0")" || exit 1
Q=quay.io/apheris/oslo-hackathon
RUN=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)
for q in chunks/q_*.json; do
  [ -e STOP ] && break
  n=$(basename "$q" .json); n=${n#q_}
  mkdir "out/lock_$n" 2>/dev/null || continue
  while [ "$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)" -lt 20000 ]; do [ -e STOP ] && break 2; echo "$(date +%T) GPU memory tight, waiting"; sleep 15; done
  echo "$(date +%T) of3 $n"; [ -d out/of3_$n ] || "${RUN[@]}" $Q:apheris-openfold3-v1.3.0-of3ob python -m apheris_openfold3.cli predict --input "$q" --output out/of3_$n --template-paths template_her2.cif --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_of3_$n.txt 2>&1
  while [ "$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)" -lt 12000 ]; do [ -e STOP ] && break 2; sleep 10; done
  echo "$(date +%T) esm2 $n"; [ -d out/esm2_$n ] || "${RUN[@]}" $Q:apheris-esmfold2-v1.4.0-release python -m apheris_esmfold2.cli predict --input "$q" --output out/esm2_$n --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_esm2_$n.txt 2>&1
done
echo "$(date +%T) worker finished"
