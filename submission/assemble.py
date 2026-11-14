"""Assemble the team submission CSV (one file, both tracks) and validate it against the organisers' format rules.

    python submission/assemble.py        # writes ../fabsolutely_binding.csv (and docs/) from submission/parts/*.csv

Inputs: parts/optimization.csv (the 15 optimization rows as first submitted, ids d01..d15) and the de novo part files
parts/{moritz_final,krisha_final,guled_final}.csv (5 + 5 + 5, with the per-design evidence and caveat columns). Output columns: sequence_id,kind,heavy_sequence,light_sequence,design_description,track.
IDs: opt-J01.., dn-M01..M05, dn-K01..K05, dn-G01..G05 (a one-letter source code that only the team can decode).
"""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=str(ROOT.parent / "fabsolutely_binding.csv"))
ap.add_argument("--opt", default=str(ROOT / "parts" / "optimization.csv"))  # the 15 rows as first submitted (ids d01..d15)
a = ap.parse_args()
AA = set("ACDEFGHIKLMNPQRSTVWY")
rows = []
for i, r in enumerate(csv.DictReader(open(a.opt)), 1):
    rows.append({"sequence_id": f"opt-J{i:02d}", "kind": r["kind"], "heavy_sequence": r["heavy_sequence"], "light_sequence": r["light_sequence"],
                 "design_description": r["design_description"], "track": "optimization"})
for code, name in (("M", "moritz_final"), ("K", "krisha_final"), ("G", "guled_final")):
    p = ROOT / "parts" / f"{name}.csv"
    if not p.exists():
        print(f"MISSING part file {p}", file=sys.stderr)
        continue
    for j, r in enumerate(sorted(csv.DictReader(open(p)), key=lambda x: int(x["src_rank"])), 1):
        rows.append({"sequence_id": f"dn-{code}{j:02d}", "kind": r["kind"], "heavy_sequence": r["heavy_sequence"], "light_sequence": r["light_sequence"],
                     "design_description": r["design_description"], "track": "denovo"})
bad = []
seen = {}
for r in rows:
    for col in ("heavy_sequence", "light_sequence"):
        s = r[col]
        if s and (set(s) - AA):
            bad.append((r["sequence_id"], col, "non-standard residues " + "".join(sorted(set(s) - AA))))
        if len(s) > 250:
            bad.append((r["sequence_id"], col, f"length {len(s)} > 250"))
        if s and s.count("C") % 2:
            bad.append((r["sequence_id"], col, f"odd cysteine count {s.count('C')}"))
    if r["kind"] not in ("vhh", "scfv", "paired"):
        bad.append((r["sequence_id"], "kind", r["kind"]))
    if r["kind"] != "paired" and r["light_sequence"]:
        bad.append((r["sequence_id"], "light_sequence", "must be empty unless paired"))
    if r["kind"] == "paired" and not r["light_sequence"]:
        bad.append((r["sequence_id"], "light_sequence", "paired needs a light chain"))
    key = (r["heavy_sequence"], r["light_sequence"])
    if key in seen:
        bad.append((r["sequence_id"], "duplicate", f"identical to {seen[key]}"))
    seen[key] = r["sequence_id"]
    if re.search(r"[\n\r]", r["design_description"]):
        bad.append((r["sequence_id"], "design_description", "newline"))
ids = [r["sequence_id"] for r in rows]
assert len(ids) == len(set(ids)), "duplicate sequence_id"
for out_path in (a.out, str(Path(a.out).parent / "docs" / "fabsolutely_binding.csv")):  # docs/ copy: GitHub Pages cannot read outside docs/
  if not Path(out_path).parent.exists():
    continue
  with open(out_path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description", "track"])
    w.writeheader()
    w.writerows(rows)
n = {t: sum(r["track"] == t for r in rows) for t in ("optimization", "denovo")}
print(f"wrote {a.out}: {len(rows)} rows {n}; kinds {sorted({r['kind'] for r in rows})}")
print("format checks:", "OK" if not bad else "")
for b in bad:
    print("  FAIL", *b)
