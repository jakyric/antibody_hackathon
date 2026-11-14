#!/usr/bin/env bash
# One-command setup for the de novo VHH pipeline. Idempotent; never prints a secret and never spends anything.
#
#   bash denovo/setup.sh            # create denovo/.venv, check keys and targets, run the no-cost smoke test
#   bash denovo/setup.sh --check    # read-only: report what is present/missing, change nothing (exit 1 only on hard failures)
#   bash denovo/setup.sh --tools    # additionally create the micromamba envs for NanoMelt and TNP (CPU; ~5 GB download)
#   bash denovo/setup.sh --panel    # additionally build the known-VHH panel for the novelty gate (public downloads, ~1 GB)
#
# Environment variables (presence is checked, values are never printed):
#   JAPANFOLD_API_KEY   required for design (RFD4) and scoring (OpenFold3, OpenDDE, ESMFold-2)
#   NEBIUS_CREDS        optional: path of the Nebius credentials file (Protenix-v2); needs the nebius-scientific-ai skill
#   PIPELINE_VM_HOST / PIPELINE_VM_KEY   optional: user@host and ssh key of a shared VM (abguard, monoyield); VM use stays off
#                       unless you ask for it explicitly (see run_pipeline.sh)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
CHECK=0; TOOLS=0; PANEL=0
for a in "$@"; do case "$a" in --check) CHECK=1;; --tools) TOOLS=1;; --panel) PANEL=1;; -h|--help) sed -n 2,16p "$0"; exit 0;; *) echo "unknown option $a"; exit 2;; esac; done
FAIL=0
ok()   { printf '  [ ok ] %s\n' "$*"; }
warn() { printf '  [warn] %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; FAIL=1; }
PY=denovo/.venv/bin/python

echo "== python"
if command -v python3 >/dev/null && python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
  ok "python3 $(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
else bad "python3 >= 3.10 not found"; fi

echo "== virtual environment (denovo/.venv)"
if [ -x "$PY" ]; then ok "exists"
elif [ "$CHECK" = 1 ]; then warn "not created yet: run  bash denovo/setup.sh"
else
  if command -v uv >/dev/null; then uv venv denovo/.venv --python 3.12 >/dev/null 2>&1 || uv venv denovo/.venv >/dev/null
  else python3 -m venv denovo/.venv; fi
  if command -v uv >/dev/null; then uv pip install --python "$PY" -r denovo/requirements.txt >/dev/null
  else "$PY" -m pip install -q -r denovo/requirements.txt; fi
  [ -x "$PY" ] && ok "created and requirements installed" || bad "could not create denovo/.venv"
fi
if [ -x "$PY" ]; then
  "$PY" - <<'EOF' && ok "imports: numpy pandas scipy pyarrow rapidfuzz Bio requests immunum" || bad "a required package is missing: rerun  bash denovo/setup.sh"
import numpy, pandas, scipy, pyarrow, rapidfuzz, Bio, requests, immunum  # noqa: F401
EOF
fi

echo "== keys (presence only)"
[ -n "${JAPANFOLD_API_KEY:-}" ] && ok "JAPANFOLD_API_KEY is set" || warn "JAPANFOLD_API_KEY is not set: design and scoring cannot run (get one at https://japanfold.aiand.com/account/token)"
if [ -n "${NEBIUS_CREDS:-}" ] && [ -f "${NEBIUS_CREDS}" ]; then ok "NEBIUS_CREDS points to a file"; else warn "NEBIUS_CREDS not set: the optional Protenix-v2 stage is unavailable"; fi
[ -f "$HOME/.claude/skills/nebius-scientific-ai/forge.py" ] && ok "nebius-scientific-ai skill found" || warn "nebius-scientific-ai skill not installed (optional, Protenix-v2 only)"
[ -n "${PIPELINE_VM_HOST:-}" ] && ok "PIPELINE_VM_HOST is set (VM paths still need an explicit opt-in)" || ok "no VM configured (default: VM off)"

echo "== targets (HER2 domain IV, derived from public PDB 1N8Z)"
miss=0; for f in READY her2_domIV.cif her2_domIV.a3m her2_domIV.fasta epitope.json numbering_map.json ref_1n8z.cif; do [ -s "denovo/targets/$f" ] || { miss=1; warn "missing denovo/targets/$f"; }; done
if [ "$miss" = 0 ]; then ok "targets present"
elif [ "$CHECK" = 1 ]; then bad "targets incomplete"
elif [ -f optimization/data/her2/her2_fv_1n8z.cif ] && [ -x "$PY" ]; then "$PY" denovo/targets/make_targets.py >/dev/null && ok "targets rebuilt" || bad "make_targets.py failed"
else bad "targets incomplete and optimization/data/her2/her2_fv_1n8z.cif is not available to rebuild them (see optimization/README.md, Data)"; fi

echo "== frameworks and design inputs"
for f in s1_design/frameworks/bcii10.json s1_design/frameworks/capla.json s1_design/inputs/bcii10_on_her2.cif s1_design/inputs/capla_caa_on_her2.cif; do
  [ -s "denovo/$f" ] && ok "denovo/$f" || bad "missing denovo/$f"; done

echo "== source files compile"
if python3 -m compileall -q denovo -x '(\.venv|/cache/)' >/dev/null 2>&1; then ok "all python files compile"; else bad "syntax error in a pipeline file"; fi
find denovo -name __pycache__ -type d -prune -exec rm -rf {} + 2>/dev/null

echo "== optional: NanoMelt / TNP environments"
if command -v micromamba >/dev/null; then
  for e in nanomelt tnp; do
    if micromamba env list 2>/dev/null | grep -qE "^\s*$e\s"; then ok "micromamba env $e exists"
    elif [ "$TOOLS" = 1 ] && [ "$CHECK" = 0 ]; then
      micromamba create -y -n "$e" -f "denovo/env/$e.yml" >/dev/null 2>&1 && ok "created env $e" || bad "could not create env $e"
      case $e in tnp) repo=https://github.com/oxpig/TNP; dir=tools/TNP; sha=ea619d1;; nanomelt) repo=https://gitlab.doc.ic.ac.uk/sormanni-lab/nanomelt; dir=tools/nanomelt; sha=dc6a141;; esac
      [ -d "$dir/.git" ] || git clone -q "$repo" "$dir"; git -C "$dir" checkout -q "$sha" && micromamba run -n "$e" pip install -q -e "$dir" && ok "installed $dir at $sha"
    else warn "micromamba env $e missing: bash denovo/setup.sh --tools (Tm and TNP are skipped without it, liabilities and nine rules still run)"; fi
  done
else warn "micromamba not installed: NanoMelt/TNP unavailable (optional)"; fi

echo "== optional: known-VHH panel for the novelty gate"
if [ -s denovo/s2_rules/cache/panel.csv ]; then ok "panel.csv present"
elif [ "$PANEL" = 1 ] && [ "$CHECK" = 0 ] && [ -x "$PY" ]; then
  mkdir -p denovo/s2_rules/cache
  curl -fsSL -o denovo/s2_rules/cache/pdb_seqres.txt.gz https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz && \
    "$PY" denovo/s2_rules/build_panel.py && "$PY" denovo/s2_rules/build_panel2.py && ok "panel built" || bad "panel build failed (see denovo/s2_rules/README.md)"
else warn "no panel (denovo/s2_rules/cache/panel.csv): the nearest-known-VHH identity column stays empty; bash denovo/setup.sh --panel builds it"; fi
[ -s denovo/s2_rules/cache/abguard_cdrh3.parquet ] && ok "abguard reference CDRH3s present" || \
  warn "abguard reference CDRH3s absent (organisers' abguard wheel, reference_data/): the prior-art CDRH3 check of the local replica is skipped"

echo "== no-cost smoke test"
if [ -x "$PY" ] && [ "$CHECK" = 0 ]; then
  "$PY" denovo/s1_design/run_design.py --auto --dry-run >/dev/null 2>&1 && ok "run_design.py --auto --dry-run (plans the campaign, spends nothing)" || warn "run_design.py dry-run did not complete (needs JAPANFOLD_API_KEY?)"
  "$PY" -W ignore denovo/rank.py >/dev/null 2>&1; ok "rank.py loads (no designs yet is expected)"
elif [ "$CHECK" = 1 ]; then ok "skipped in --check mode"; fi

echo
if [ "$FAIL" = 0 ]; then echo "setup OK. Next:  bash denovo/run_pipeline.sh dry-run"; else echo "setup finished with FAILURES above"; fi
exit "$FAIL"
