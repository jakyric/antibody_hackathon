#!/usr/bin/env bash
# Fold one panel with one folding model, in its own Docker image. Everyone gets the same
# panel, the same seed and the same template, so the numbers can be compared between them.
#
#   bash scripts/fold_run.sh <model> <queries.json> <output-dir> [seed]
#
# model: boltz2 | openfold3-of3ob | openfold3-of3p2 | protenix-v2 | protenix-v1 |
#        boltz1 | esmfold2 | mol1-base
# Models without a template (boltz1, esmfold2) are run without --template-paths.
set -euo pipefail

MODEL="${1:?model}"
QUERIES="${2:?queries.json}"
OUT="${3:?output-dir}"
SEED="${4:-42}"
TEMPLATE="${TEMPLATE:-$(dirname "$QUERIES")/template.cif}"

Q=quay.io/apheris/oslo-hackathon
case "$MODEL" in
  boltz2)           IMG=$Q:apheris-boltz-v1.3.0-boltz2;              TPL=1 ;;
  boltz1)           IMG=$Q:apheris-boltz-v1.3.0-boltz1;              TPL=0 ;;
  openfold3-of3ob)  IMG=$Q:apheris-openfold3-v1.3.0-of3ob;           TPL=1 ;;
  openfold3-of3p2)  IMG=$Q:apheris-openfold3-v1.3.0-of3p2;           TPL=1 ;;
  protenix-v2)      IMG=$Q:apheris-protenix-v1.3.0-v2;               TPL=1 ;;
  protenix-v1)      IMG=$Q:apheris-protenix-v1.3.0-v1-default;       TPL=1 ;;
  mol1-base)        IMG=$Q:apheris-mol1-v1.3.0-base;                 TPL=1 ;;
  esmfold2)         IMG=$Q:apheris-esmfold2-v1.4.0-release;          TPL=0 ;;
  *) echo "unknown model: $MODEL" >&2; exit 2 ;;
esac

mkdir -p "$OUT"
CMD=$(docker inspect -f '{{join .Config.Cmd " "}}' "$IMG")

args=(predict --input "$QUERIES" --output "$OUT"
      --write-detailed-scores true --num-diffusion-samples 1
      --seeds "$SEED" --query-failures skip)
if [ "$TPL" = 1 ]; then
  args+=(--template-paths "$TEMPLATE")
fi

echo "model=$MODEL image=$IMG seed=$SEED template=$TPL"
docker image inspect -f '{{index .RepoDigests 0}}' "$IMG" 2>/dev/null || true

# shellcheck disable=SC2086
exec docker run --rm --gpus all --shm-size 16g --user "$(id -u):$(id -g)" \
  -e HOME=/tmp -v /etc/passwd:/etc/passwd:ro -v "$PWD:/work" -w /work \
  "$IMG" $CMD "${args[@]}"

