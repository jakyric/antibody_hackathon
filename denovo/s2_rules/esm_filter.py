"""ESM-2 650M naturalness check of the diversified variants (hosted on Nebius, one job for the whole batch).

For every variant we send the parent sequence + the substitution set; ESM-2 masked-marginals returns the summed
log-odds (mutant vs parent residue) in the parent context. A variant whose substitutions the language model dislikes
(per-substitution mean strongly negative) is a likelier folding/expression loser: it is flagged `esm_ok = 0`.

  .venv/bin/python denovo/s2_rules/esm_filter.py submit  --designs results/denovo/designs_fwdiv.csv --parents results/denovo/designs_*.csv --tag fwdiv1
  .venv/bin/python denovo/s2_rules/esm_filter.py collect --tag fwdiv1 --out denovo/s2_rules/esm_scores_fwdiv1.csv

Recorded runs live in nebius-runs/ (forge.py convention); re-running `submit` follows the existing run, never a second job.
"""
from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_lib as L  # noqa: E402

dn = L.dn
FORGE = Path.home() / ".claude/skills/nebius-scientific-ai/forge.py"
ROOT = dn.ROOT
RUNS = ROOT / "nebius-runs"


def subs(parent: str, var: str) -> str:
    return ":".join(f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(parent, var)) if a != b)


def build_tar(designs_csv: str, parent_csvs: list[str], tag: str) -> Path:
    var = dn.read_csv(designs_csv)
    parents = {}
    for pat in parent_csvs:
        for f in glob.glob(pat):
            for r in dn.read_csv(f):
                parents[r["design_id"]] = r["vhh_seq"]
    work = HERE / "cache" / f"esm_{tag}"
    work.mkdir(parents=True, exist_ok=True)
    need = sorted({v["parent_id"] for v in var})
    (work / "sequences.fasta").write_text("".join(f">{p}\n{parents[p]}\n" for p in need))
    lines = ["id,sequence_id,mutations"]
    for v in var:
        m = subs(parents[v["parent_id"]], v["vhh_seq"])
        if m:
            lines.append(f'{v["design_id"]},{v["parent_id"]},{m}')
    (work / "variants.csv").write_text("\n".join(lines) + "\n")
    tar = RUNS / f"s2-esm-{tag}.tar"
    RUNS.mkdir(exist_ok=True)
    with tarfile.open(tar, "w") as t:
        for n in ("sequences.fasta", "variants.csv"):
            t.add(work / n, arcname=n)
    return tar


def submit(a):
    tar = build_tar(a.designs, a.parents, a.tag)
    cmd = [sys.executable, str(FORGE), "submit", "esm2-650m", "--operation", "score-variants",
           "--entry", f"workload|protein-sequence-scoring-tar/v1|{tar}", "--params", json.dumps({"strategy": "masked-marginals"}),
           "--label", f"s2-esm-{a.tag}", "--no-wait"]
    print(subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).stdout[-800:])


def collect(a):
    cmd = [sys.executable, str(FORGE), "runs"]
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).stdout
    op = next((l.split()[3] for l in out.splitlines() if f"s2-esm-{a.tag}" in l), None)
    if not op:
        sys.exit("run not found")
    dest = ROOT / "results" / "denovo_s2" / f"esm-{a.tag}"
    r = subprocess.run([sys.executable, str(FORGE), "fetch", op, "--dest", str(dest)], cwd=ROOT, capture_output=True, text=True)
    print(r.stdout[-300:], r.stderr[-300:])
    print("fetched to", dest)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["submit", "collect"])
    ap.add_argument("--designs")
    ap.add_argument("--parents", nargs="*", default=[str(dn.OUT / "designs_*.csv")])
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    {"submit": submit, "collect": collect}[a.cmd](a)
