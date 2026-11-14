#!/usr/bin/env python3
"""Looks for secrets and oversized files before anything is published to a public repo.

Checks what git tracks and what would be committed (not ignored files):
  - files that look like key or login files
  - private keys and known key formats in text files
  - files over 50 MB

Never prints the value it found, only the file and line. Exit code 1 if anything was found.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

MAX_BYTES = 50 * 1024 * 1024
SECRET_NAMES = re.compile(r"(^|/)(\.env(\..*)?|id_(rsa|ed25519|ecdsa|dsa)(\.pub)?|.*\.pem|.*\.p12|.*\.pfx|.*\.private\.json|.*credentials.*\.json|.*secrets?\..*)$", re.I)
PATTERNS = [
    ("private key", re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----")),
    ("AWS key", re.compile(r"\bAKIA[0-9A-Z]{16}\b")),
    ("API key (sk-)", re.compile(r"\bsk-[A-Za-z0-9_-]{20,}")),
    ("GitHub token", re.compile(r"\bgh[pousr]_[A-Za-z0-9]{30,}")),
    ("key in field", re.compile(r"""["']?(api[_-]?key|secret|token|password|passwd)["']?\s*[:=]\s*["'][^"'\s]{12,}["']""", re.I)),
]
SKIP_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".zip", ".gz", ".zst", ".tar", ".parquet", ".pdb", ".cif", ".npy", ".npz", ".pkl", ".ico"}


def git_files(root: Path) -> list[str]:
    tracked = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, text=True, check=True).stdout
    others = subprocess.run(["git", "ls-files", "-z", "--others", "--exclude-standard"], cwd=root, capture_output=True, text=True, check=True).stdout
    return sorted({f for f in (tracked + others).split("\0") if f})


def scan(root: Path) -> list[str]:
    findings = []
    for name in git_files(root):
        path = root / name
        if not path.is_file():
            continue
        if SECRET_NAMES.search(name):
            findings.append(f"{name}: looks like a key or login file")
        size = path.stat().st_size
        if size > MAX_BYTES:
            findings.append(f"{name}: {size // (1024 * 1024)} MB, over 50 MB")
            continue
        if path.suffix.lower() in SKIP_SUFFIXES or size > 5 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(errors="ignore").splitlines()
        except OSError:
            continue
        for number, line in enumerate(lines, 1):
            for label, pattern in PATTERNS:
                if pattern.search(line):
                    findings.append(f"{name}:{number}: {label}")
    return findings


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    try:
        findings = scan(root)
    except subprocess.CalledProcessError:
        print("This is not a git repo.", file=sys.stderr)
        return 2
    for finding in findings:
        print(finding)
    if findings:
        print(f"\n{len(findings)} findings. Nothing should be published until they are fixed.", file=sys.stderr)
        return 1
    print("No findings.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
