#!/usr/bin/env bash
# Runs ON THE VM, one GPU container at a time, chunk by chunk (20 designs each): templated OpenFold3 (OpenBind weights, HER2 template)
# then untemplated ESMFold-2, so partial progress is usable at any moment. Stop with `touch STOP` in this folder.
cd "$(dirname "$0")" || exit 1
Q=quay.io/apheris/oslo-hackathon
RUN=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)
for q in chunks/q_*.json; do
  [ -e STOP ] && echo "STOP" && break
  n=$(basename "$q" .json); n=${n#q_}
  FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "${FREE:-0}" -lt 30000 ]; then echo "$(date +%T) GPU busy (free ${FREE} MiB), waiting"; sleep 20; fi
  [ -d out/of3_$n ] || { echo "$(date +%T) of3 $n"; "${RUN[@]}" $Q:apheris-openfold3-v1.3.0-of3ob python -m apheris_openfold3.cli predict --input "$q" --output out/of3_$n --template-paths template_her2.cif --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_of3_$n.txt 2>&1; }
  [ -d out/esm2_$n ] || { echo "$(date +%T) esm2 $n"; "${RUN[@]}" $Q:apheris-esmfold2-v1.4.0-release python -m apheris_esmfold2.cli predict --input "$q" --output out/esm2_$n --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_esm2_$n.txt 2>&1; }
done
echo "$(date +%T) finished"
