"""Stage 3: organisers' novelty / liability gate (`abguard team-validate`).

abguard compares the whole CDR-H3 (edit distance, so insertions and deletions count) with every known CDR-H3 in its
bundled reference (610k sequences) and requires an even number of heavy-chain Cys. Rows that fail are dropped.
Modes: ssh (run on the machine that holds the build), local (binary on this machine), skip (no gate; for tests).
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from .common import out, save_summary


def run(cfg: dict) -> None:
    v = cfg["validate"]
    src, dst = out(cfg, "panel_in.csv"), out(cfg, "panel_valid.csv")
    n_in = len(pd.read_csv(src))
    if v["mode"] == "skip":
        shutil.copy(src, dst)
    elif v["mode"] == "local":
        log = subprocess.run([v["abguard_bin"], "team-validate", str(src), "--drop", "-o", str(dst)],
                             capture_output=True, text=True).stderr
        out(cfg, "abguard.log").write_text(log)
    elif v["mode"] == "ssh":
        key = str(Path(v["ssh_key"]).expanduser())
        ssh = ["ssh", "-i", key, "-o", "BatchMode=yes", v["ssh_host"]]
        wd = v["workdir"]
        subprocess.run(ssh + [f"mkdir -p {wd}"], check=True)
        subprocess.run(["scp", "-q", "-i", key, str(src), f"{v['ssh_host']}:{wd}/panel_in.csv"], check=True)
        cmd = f"cd {wd} && {v['abguard_bin']} team-validate panel_in.csv --drop -o panel_valid.csv > abguard.log 2>&1; cat abguard.log"
        log = subprocess.run(ssh + [cmd], capture_output=True, text=True).stdout
        out(cfg, "abguard.log").write_text(log)
        subprocess.run(["scp", "-q", "-i", key, f"{v['ssh_host']}:{wd}/panel_valid.csv", str(dst)], check=True)
    else:
        raise ValueError(f"validate.mode {v['mode']!r}")
    valid = pd.read_csv(dst)
    failed = set()
    log_f = out(cfg, "abguard.log")
    if log_f.exists():
        failed = set(re.findall(r"\[novelty\] ([A-Z]+)_", log_f.read_text()))
    print(f"  abguard: {len(valid)}/{n_in} sequences valid" + (f"; novelty failures on H3: {sorted(failed)}" if failed else ""))
    # a CDR-H3 counts as valid when its all-WT-Vernier genotype passes (abguard's CDR-3 string includes the H94 anchor,
    # so a Vernier swap at H94 can turn a failing loop into a passing sequence)
    wt_rows = valid[valid.sequence_id.str.endswith("_WT")]
    save_summary(cfg, n_valid=len(valid), n_valid_h3=int(len(wt_rows)))
