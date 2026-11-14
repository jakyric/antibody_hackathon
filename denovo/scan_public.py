#!/usr/bin/env python3
"""Pre-commit safety scan for this PUBLIC repository (run by sync/final; run it yourself before any push).

    python3 denovo/scan_public.py [repo-root]        # exit 1 if anything suspicious is found

Checks every text file under denovo/, .claude/, results_denovo/ and the root README/.gitignore for: key-like strings,
bearer tokens, the value of $JAPANFOLD_API_KEY, VM addresses and key paths, absolute home paths, credential-file names,
personal e-mail addresses, and files larger than 2 MB. Prints file:line and the rule name, never the matching text.
"""
import os
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
RULES = {
    "vm-address": r"\b\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}\b",
    "ssh-key-path": r"\.ssh/(?!cm-)[A-Za-z0-9_.-]+",  # cm-* are ssh ControlPath sockets, not keys
    "absolute-home-path": r"/(home|Users)/[A-Za-z0-9_.-]+/",
    "credential-file": r"credentials\.private|participant-0\d|settings\.local\.json",
    "bearer-token": r"Bearer\s+[A-Za-z0-9._\-]{12,}",
    "secret-assignment": r"(?i)(api[_-]?key|secret|token|passw(or)?d)\s*[:=]\s*[\"'][^\"'\s$]{8,}[\"']",
    "known-key-format": r"AKIA[0-9A-Z]{16}|ghp_[A-Za-z0-9]{30,}|sk-[A-Za-z0-9]{20,}|-----BEGIN [A-Z ]*PRIVATE KEY-----",
    "email": r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+\.[A-Za-z]{2,}",
}
ALLOW_EMAIL = {"noreply@anthropic.com", "git@github.com"}
SKIP_PARTS = {".git", ".venv", "__pycache__", "cache", "runs", "raw"}
TEXT = {".py", ".sh", ".md", ".txt", ".json", ".csv", ".yml", ".yaml", ".html", ".css", ".toml", ""}
key = os.environ.get("JAPANFOLD_API_KEY", "")
files = []
for top in ("denovo", ".claude", "results_denovo", "README.md", ".gitignore"):
    p = root / top
    if p.is_file():
        files.append(p)
    elif p.is_dir():
        files += [q for q in p.rglob("*") if q.is_file() and not (set(q.relative_to(root).parts) & SKIP_PARTS)]
bad = 0
for f in sorted(files):
    rel = f.relative_to(root)
    if rel.as_posix() == "denovo/scan_public.py":
        continue
    if f.stat().st_size > 2_000_000:
        print(f"{rel}: LARGE FILE ({f.stat().st_size / 1e6:.1f} MB)")
        bad += 1
    if f.suffix not in TEXT:
        continue
    try:
        lines = f.read_text(errors="replace").splitlines()
    except Exception:
        continue
    for i, line in enumerate(lines, 1):
        if len(key) >= 8 and key in line:
            print(f"{rel}:{i}: contains the JapanFold key value")
            bad += 1
        for name, pat in RULES.items():
            for m in re.finditer(pat, line):
                s = m.group(0)
                if name == "email" and s in ALLOW_EMAIL:
                    continue
                if name == "vm-address" and f.suffix in (".yml", ".yaml") and re.search(r"=\s*" + re.escape(s), line):
                    continue  # conda spec versions such as libcufile=1.14.1.1 look like IPv4 addresses
                if name == "credential-file" and rel.name == ".gitignore":
                    continue  # an ignore RULE naming the file is the safeguard, not a leak
                if name == "vm-address" and (s.startswith(("0.", "127.", "10.")) or re.fullmatch(r"\d+\.\d+\.\d+\.\d+", s) and
                                             re.search(r"[A-Za-z]\s*" + re.escape(s) + r"|" + re.escape(s) + r"\s*[A-Za-z]", line)):
                    continue  # version strings / numbers such as 1.2.3.4 inside prose are rare; keep the filter conservative
                print(f"{rel}:{i}: {name}")
                bad += 1
print(f"scan: {len(files)} files, {bad} finding(s)")
sys.exit(1 if bad else 0)
