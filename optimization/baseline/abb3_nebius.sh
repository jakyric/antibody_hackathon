#!/usr/bin/env bash
# Model every nebius-runs/paper/inputs/p*.fasta (H + L chains) with ABodyBuilder3 on Nebius, 4 at a time, then fetch.
# Run from anywhere in the repo:  bash optimization/baseline/abb3_nebius.sh   (the script changes into optimization/ itself)
set -u
F=~/.claude/skills/nebius-scientific-ai/forge.py
cd "$(dirname "$0")/.."
mkdir -p results/paper/abb3
ls nebius-runs/paper/inputs/p*.fasta | xargs -n1 basename | sed 's/\.fasta$//' | \
  xargs -P 4 -I@@ sh -c "python3 $F submit abodybuilder3 --entry 'sequences|antibody-fasta/v1|nebius-runs/paper/inputs/@@.fasta' --params '{}' --label abb3-@@ > results/paper/abb3/@@.log 2>&1"
echo "all submitted runs finished; fetching"
for id in $(ls nebius-runs/paper/inputs/p*.fasta | xargs -n1 basename | sed 's/\.fasta$//'); do
  op=$(grep -o 'op-[0-9a-f]\{16\}' results/paper/abb3/$id.log | head -1)
  [ -n "$op" ] && python3 $F fetch $op --dest results/paper/abb3/$id > results/paper/abb3/$id.fetch.log 2>&1 && echo "fetched $id ($op)"
done
