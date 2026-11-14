#!/usr/bin/env bash
# Fold every complex in nebius-runs/refold${SUF}/ with Boltz-2 on the Nebius platform, 4 at a time, then fetch the results.
# Run from the project root:  bash baseline/refold_nebius.sh
# Needs NEBIUS_CREDS (set by .claude/settings.local.json) and the nebius-scientific-ai skill.
set -u
F=~/.claude/skills/nebius-scientific-ai/forge.py
cd "$(dirname "$0")/.."
SUF=${VARIANT:+_$VARIANT}   # VARIANT=contacts reruns into refold_contacts/
mkdir -p results/refold${SUF}/nebius
ls nebius-runs/refold${SUF}/*.json | xargs -n1 basename | sed 's/\.json$//' | \
  xargs -P 4 -I{} sh -c "python3 $F run boltz2 --args-file nebius-runs/refold${SUF}/{}.json --label refold${SUF}-{} > results/refold${SUF}/nebius/{}.log 2>&1"
echo "all submitted runs finished; fetching"
for id in $(ls nebius-runs/refold${SUF}/*.json | xargs -n1 basename | sed 's/\.json$//'); do
  op=$(grep -o 'op-[0-9a-f]\{16\}' results/refold${SUF}/nebius/$id.log | head -1)
  [ -n "$op" ] && python3 $F fetch $op --dest results/refold${SUF}/nebius/$id > results/refold${SUF}/nebius/$id.fetch.log 2>&1 && echo "fetched $id ($op)"
done
