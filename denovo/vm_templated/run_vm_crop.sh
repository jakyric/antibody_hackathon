#!/usr/bin/env bash
# Parallel worker (ON THE VM) for the cropped-HER2 re-run: claims the next unclaimed chunk in crop/chunks (mkdir lock) and runs templated
# OpenFold3 (cropped template, no MSA) on it. Waits whenever less than 20 GB of GPU memory is free, so it never pushes other people's jobs out
# of memory. `touch crop/STOP` stops all crop workers.  Start N of them:  for i in 1 2 3 4; do nohup bash crop/run_vm_crop.sh > crop/w$i.log 2>&1 & done
cd "$(dirname "$0")/.." || exit 1
Q=quay.io/apheris/oslo-hackathon
RUN=(docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work)
for q in crop/chunks/q_*.json; do
  [ -e crop/STOP ] && break
  n=$(basename "$q" .json); n=${n#q_}
  mkdir "out/croplock_$n" 2>/dev/null || continue
  while [ "$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)" -lt 20000 ]; do [ -e crop/STOP ] && break 2; echo "$(date +%T) GPU memory tight, waiting"; sleep 15; done
  echo "$(date +%T) of3 crop $n"
  [ -d out/of3c_$n ] || "${RUN[@]}" $Q:apheris-openfold3-v1.3.0-of3ob python -m apheris_openfold3.cli predict --input "$q" --output out/of3c_$n --template-paths crop/template_her2_crop.cif --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip > logs_of3c_$n.txt 2>&1
done
echo "$(date +%T) worker finished"
