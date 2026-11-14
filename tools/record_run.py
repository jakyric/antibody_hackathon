#!/usr/bin/env python3
"""Run a command and record everything needed to repeat it.

Usage:
    record_run.py --name <slug> [--seed N] [--input FILE ...] [--note TEXT] -- <command ...>

Creates runs/<UTC-time>-<slug>/ with run.json, stdout.log and stderr.log, and appends a line
to results/ledger.csv. The command gets $RUN_DIR and $RUN_SEED. Environment variables are not
saved, since they may contain keys.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

LEDGER_FIELDS = ["run_id", "name", "started", "seconds", "exit_code", "git_sha", "dirty", "seed", "note", "command"]


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state(cwd: Path) -> tuple[str | None, bool | None]:
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip())
        return sha, dirty
    except (OSError, subprocess.CalledProcessError):
        return None, None


def slugify(name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "-", name).strip("-")
    if not slug:
        raise ValueError("--name must contain at least one character that is a letter or a digit")
    return slug[:60]


def pump(source, sinks) -> None:
    for line in iter(source.readline, b""):
        for sink in sinks:
            sink.write(line)
            sink.flush()


def record(name: str, command: list[str], root: Path, seed: int | None = None, inputs: list[str] | None = None, note: str = "") -> tuple[Path, int]:
    started = datetime.now(timezone.utc)
    run_id = f"{started.strftime('%Y%m%dT%H%M%SZ')}-{slugify(name)}"
    run_dir = root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    input_hashes = {}
    for item in inputs or []:
        path = (root / item) if not Path(item).is_absolute() else Path(item)
        input_hashes[item] = sha256_of(path) if path.is_file() else None

    sha, dirty = git_state(root)
    manifest = {
        "run_id": run_id,
        "name": name,
        "command": command,
        "cwd": str(root),
        "git_sha": sha,
        "dirty": dirty,
        "seed": seed,
        "inputs": input_hashes,
        "note": note,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "started": started.isoformat(),
    }
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")

    env = dict(os.environ, RUN_DIR=str(run_dir), RUN_SEED="" if seed is None else str(seed))
    begin = time.monotonic()
    with (run_dir / "stdout.log").open("wb") as out_log, (run_dir / "stderr.log").open("wb") as err_log:
        process = subprocess.Popen(command, cwd=root, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        threads = [
            threading.Thread(target=pump, args=(process.stdout, [out_log, sys.stdout.buffer])),
            threading.Thread(target=pump, args=(process.stderr, [err_log, sys.stderr.buffer])),
        ]
        for thread in threads:
            thread.start()
        code = process.wait()
        for thread in threads:
            thread.join()
    seconds = round(time.monotonic() - begin, 1)

    manifest.update({"ended": datetime.now(timezone.utc).isoformat(), "seconds": seconds, "exit_code": code})
    (run_dir / "run.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")

    ledger = root / "results" / "ledger.csv"
    ledger.parent.mkdir(parents=True, exist_ok=True)
    fresh = not ledger.exists()
    with ledger.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS)
        if fresh:
            writer.writeheader()
        writer.writerow({
            "run_id": run_id, "name": name, "started": manifest["started"], "seconds": seconds, "exit_code": code,
            "git_sha": sha or "", "dirty": "" if dirty is None else int(dirty), "seed": "" if seed is None else seed,
            "note": note, "command": " ".join(command),
        })
    return run_dir, code


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--" not in argv:
        print("Missing «--» before the command. Usage: record_run.py --name <slug> -- <command ...>", file=sys.stderr)
        return 2
    split = argv.index("--")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", required=True)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--note", default="")
    parser.add_argument("--root", default=".", help="the project directory (default: current directory)")
    args = parser.parse_args(argv[:split])
    command = argv[split + 1:]
    if not command:
        print("No command after «--».", file=sys.stderr)
        return 2
    run_dir, code = record(args.name, command, Path(args.root).resolve(), args.seed, args.input, args.note)
    print(f"\nThe run was saved in {run_dir.relative_to(Path(args.root).resolve())} (exit code {code})", file=sys.stderr)
    return code


if __name__ == "__main__":
    sys.exit(main())
