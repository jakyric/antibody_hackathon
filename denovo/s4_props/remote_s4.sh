#!/usr/bin/env bash
# Runs ON THE VM.  Usage: remote_s4.sh <workdir> <tools> [dtype] [ncores]
#   workdir : contains in.csv (design_id,vhh_seq[,cdr1,cdr2,cdr3]); results are written next to it
#   tools   : comma list of nanomelt,tnp,liab,nb9,monoyield   (default all)
#   dtype   : auto (fp32 if >= 45 GB GPU memory is free, else bf16) | fp32 | bf16     (monoyield only)
#   ncores  : TNP downstream pool size (default 6)
# One GPU process at a time; everything under ~/moritz.
set -uo pipefail
W="$(cd "$1" && pwd)"; TOOLS="${2:-nanomelt,tnp,liab,nb9,monoyield}"; DT="${3:-auto}"; NC="${4:-6}"
S="$HOME/moritz/s4"; MM="$S/bin/micromamba"; export MAMBA_ROOT_PREFIX="$S/mamba"
PY="$MM run -n s4"
export TORCH_HOME="$S/torch_cache" HF_HOME="$S/hf_cache" OMP_NUM_THREADS=4
cd "$W"
echo "start $(date -u +%FT%TZ) tools=$TOOLS" > status.txt
has() { [[ ",$TOOLS," == *",$1,"* ]]; }

if has liab; then
  $PY python "$S/denovo/s4_props/liab.py" --in in.csv --out liab.csv > liab.log 2>&1 && echo "liab ok" >> status.txt || echo "liab FAILED" >> status.txt
fi
if has nb9; then
  $PY python "$S/denovo/s4_props/nanobody_rules.py" --in in.csv --out nb9.csv > nb9.log 2>&1 && echo "nb9 ok" >> status.txt || echo "nb9 FAILED" >> status.txt
fi
if has tnp; then
  CUDA_VISIBLE_DEVICES= $PY python "$S/denovo/s4_props/tnp_batch.py" --in in.csv --out tnp.csv --pdb-dir "$W/tnp_pdbs" --ncores "$NC" > tnp.log 2>&1 \
    && echo "tnp ok" >> status.txt || echo "tnp FAILED" >> status.txt
fi
if has nanomelt; then
  $PY python - <<'PYEOF' > nanomelt.fa
import csv
for r in csv.DictReader(open("in.csv")):
    print(f">{r['design_id']}\n{r['vhh_seq']}")
PYEOF
  split -l 100 -d nanomelt.fa nmchunk_      # 50 sequences per chunk (2 fasta lines each)
  rm -f nanomelt.csv; ok=1
  for c in nmchunk_??; do
    PYTHONPATH="$S/tools/nanomelt" $PY python -m nanomelt predict -i "$c" -align -o "$c.csv" >> nanomelt.log 2>&1 || ok=0
    if [ -f "$c.csv" ]; then if [ -f nanomelt.csv ]; then tail -n +2 "$c.csv" >> nanomelt.csv; else cp "$c.csv" nanomelt.csv; fi; fi
  done
  [ $ok = 1 ] && echo "nanomelt ok" >> status.txt || echo "nanomelt FAILED (partial results kept)" >> status.txt
fi
if has monoyield; then
  FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
  if [ "$DT" = auto ]; then [ "${FREE:-0}" -ge 45000 ] && DT=fp32 || DT=bf16; fi
  [ "$DT" = fp32 ] && DTARG=float32 || DTARG=bfloat16
  echo "monoyield dtype=$DT gpu_free_MiB=$FREE" >> status.txt
  $PY python - <<'PYEOF'
import csv
rows = list(csv.DictReader(open("in.csv")))
with open("my_in.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["name", "seq"])
    for r in rows: w.writerow([r["design_id"], r["vhh_seq"]])
PYEOF
  (
    source "$HOME/moritz/monoyield-predictor/env.sh"
    cd "$HOME/moritz/monoyield-predictor"
    REV=45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a
    for M in monoyield mdratio; do
      python scripts_esm34/predict_esmc_residue.py --external_csv "$W/my_in.csv" --seq_col seq \
        --model_dir "weights/esmc_nb_$M" --model_name biohub/ESMC-6B --model_revision "$REV" \
        --dtype "$DTARG" --batch_size 4 --out_csv "$W/my_$M.csv" > "$W/my_$M.log" 2>&1 \
        && echo "monoyield-$M ok ($DT)" >> "$W/status.txt" || echo "monoyield-$M FAILED" >> "$W/status.txt"
    done
  )
fi
echo "done $(date -u +%FT%TZ)" >> status.txt
