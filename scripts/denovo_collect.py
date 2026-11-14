#!/usr/bin/env python3
"""Collect the de novo designs from several runs, with provenance, and filter them.

Each run (weight + seed) lives in its own runs/<id>/ folder with run.json.
Here they are merged into a single table where every design can be traced
back to model, weight, seed and run id, and where the novelty numbers sit
next to it.

THE FILTERING, in order (a design must pass all of them):
  1. the whole construct is legal: no unpaired cysteines, scFv under 250
  2. the CDR identity to trastuzumab is under 70 %
  3. CDR-H3 lies at least 3 point mutations from EVERYTHING in the provided data
     (the de novo requirement is 2; we use 3, the same margin as the
     optimisation track, because the lab runs a larger reference set than
     the official checker includes)
  4. identical constructs are counted once

Item 3 is "filter out known binders": the cautious settings produce
loops that are already in the FACS screen, and they are discarded here,
before anything is folded.

    python scripts/denovo_collect.py --runs-dir runs --out results/denovo
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from abx.data import TZ_CDRH3_10, provided_cdrh3s  # noqa: E402
from abx.gate import ProvidedIndex, count_unpaired_cysteines  # noqa: E402
from abx.sequences import TRASTUZUMAB_VH  # noqa: E402

from denovo_gate import evaluate  # noqa: E402

MIN_PROVIDED_DIST = 3      # de novo requires 2; we require 3 for margin
MAX_CDR_IDENTITY = 70.0    # strictly below
MAX_SCFV_LEN = 250


def run_metadata(run_dir: Path) -> dict:
    """Weight, seed and run id, read from run.json — not guessed from the name."""
    meta = {"run_id": run_dir.name, "weight": "", "seed": "", "commit": ""}
    rj = run_dir / "run.json"
    if rj.exists():
        info = json.loads(rj.read_text())
        meta["seed"] = str(info.get("seed", ""))
        meta["commit"] = str(info.get("git_sha") or "")[:12]
        name = str(info.get("name", run_dir.name))
        # the name is denovo-<weight>-s<seed>
        if name.startswith("denovo-") and "-s" in name:
            meta["weight"] = name[len("denovo-"):name.rindex("-s")]
    return meta


AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def sequences_from_cif(path: Path) -> str:
    """"A:...;B:..." read straight out of a generated structure.

    The metrics table is only written once the whole run is finished. If a
    run is interrupted, the designs are still in inference_samples/. The
    sequences are then read from there instead, so the work is not lost.
    """
    cols: list[str] = []
    chains: dict[str, list[str]] = {}
    seen: set[tuple[str, str]] = set()
    for line in path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1])
            continue
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        if len(f) < len(cols):
            continue
        rec = dict(zip(cols, f))
        chain = rec.get("auth_asym_id") or rec.get("label_asym_id")
        num = rec.get("auth_seq_id") or rec.get("label_seq_id")
        comp = rec.get("label_comp_id") or rec.get("auth_comp_id", "")
        key = (chain, num)
        if key in seen:
            continue
        seen.add(key)
        chains.setdefault(chain, []).append(AA3.get(comp.upper(), "X"))
    return ";".join(f"{c}:{''.join(s)}" for c, s in sorted(chains.items()))


def designs_in(run_dir: Path):
    metrics = run_dir / "design" / "generation" / "inference_metrics" / "per_example_metrics.csv"
    if metrics.exists():
        with metrics.open(newline="") as fh:
            yield from csv.DictReader(fh)
        return
    samples = sorted((run_dir / "design" / "generation" / "inference_samples").glob("*/*_generated.cif"))
    if samples:
        print(f"  {run_dir.name}: no metrics table, reading {len(samples)} "
              f"structures directly (interrupted run)")
    for cif in samples:
        yield {"example_id": cif.parent.name, "sequences": sequences_from_cif(cif)}


def parse_chains(field: str) -> dict[str, str]:
    out = {}
    for part in field.split(";"):
        if ":" in part:
            cid, seq = part.split(":", 1)
            out[cid.strip()] = seq.strip()
    return out


def evaluate_vhh(h: str) -> dict:
    """The numbers for a VHH. One chain, its own framework, its own novelty assessment.

    A VHH does not sit on trastuzumab's framework, and that is the whole
    point: the overall identity to trastuzumab is what the official checker
    uses to decide the track, and it should be well below 70 %. CDR3 is
    taken between the same conserved anchors as in any VH domain.
    """
    from abx.gate import identity, levenshtein

    i = h.find("YYC")
    j = h.find("WGQG", i + 3) if i >= 0 else -1
    cdr3 = h[i + 3:j] if i >= 0 and j > i else ""
    return {
        "vh_len": len(h),
        "vl_len": 0,
        "scfv_len": len(h),
        "unpaired_cys": count_unpaired_cysteines(h),
        "identity_whole_pct": round(100 * identity(h, TRASTUZUMAB_VH), 1),
        "identity_cdr_pct": "",          # not defined against a different framework
        "cdrh3": cdr3,
        "cdrh3_min_dist_provided": levenshtein(cdr3, TZ_CDRH3_10) if cdr3 else "",
        "cdrh3_identity_tz_pct": round(100 * identity(cdr3, TZ_CDRH3_10), 1) if cdr3 else "",
    }


def verdict_vhh(row: dict) -> tuple[bool, str]:
    """The filtering for the VHH track.

    The requirements are the de novo requirements: under 70 % identity to a
    known antibody (here measured against trastuzumab; the official checker
    has the reference database), at least 2 CDR3 mutations from
    trastuzumab — we require 3 for margin — no unpaired cysteines, and
    under 250 residues.
    """
    if row["unpaired_cys"]:
        return False, "unpaired cysteine"
    if int(row["scfv_len"]) > MAX_SCFV_LENGTH_GUARD:
        return False, f"longer than {MAX_SCFV_LENGTH_GUARD}"
    if float(row["identity_whole_pct"]) >= MAX_CDR_IDENTITY:
        return False, "sequence identity to trastuzumab is not under 70 %"
    if not row["cdrh3"]:
        return False, "no CDR3 found between the framework anchors"
    if int(row["cdrh3_min_dist_provided"]) < MIN_PROVIDED_DIST:
        return False, "CDR3 too close to trastuzumab's"
    return True, ""


MAX_SCFV_LENGTH_GUARD = MAX_SCFV_LEN


def verdict(row: dict) -> tuple[bool, str]:
    """Does the design pass the filtering? Also return the reason it failed."""
    if row["unpaired_cys"]:
        return False, "unpaired cysteine"
    if row["scfv_len"] == "" or int(row["scfv_len"]) > MAX_SCFV_LEN:
        return False, f"scFv longer than {MAX_SCFV_LEN}"
    if row["identity_cdr_pct"] == "" or float(row["identity_cdr_pct"]) >= MAX_CDR_IDENTITY:
        return False, "CDR identity to trastuzumab not under 70 %"
    if not row["cdrh3"]:
        return False, "no CDRH3 found between the framework anchors"
    d = row["cdrh3_min_dist_provided"]
    if d == "" or int(d) < MIN_PROVIDED_DIST:
        return False, f"known binder: CDRH3 only {d} mutation(s) from provided data"
    return True, ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", type=Path, default=Path("runs"))
    ap.add_argument("--pattern", default="denovo-",
                    help="runs with this in the name; several patterns separated by commas")
    ap.add_argument("--out", type=Path, default=Path("results/denovo"))
    args = ap.parse_args()

    patterns = [p.strip() for p in args.pattern.split(",") if p.strip()]
    run_dirs = sorted(
        d for d in args.runs_dir.iterdir()
        if d.is_dir() and any(p in d.name for p in patterns)
        and ("smoke" not in d.name or any("smoke" in p for p in patterns))
    )
    if not run_dirs:
        raise SystemExit(f"no runs with '{args.pattern}' in {args.runs_dir}")

    novelty = ProvidedIndex(provided_cdrh3s())
    rows: list[dict] = []
    for run_dir in run_dirs:
        meta = run_metadata(run_dir)
        for rec in designs_in(run_dir):
            chains = parse_chains(rec.get("sequences", ""))
            vl, vh = chains.get("A"), chains.get("B")
            h = chains.get("H")
            is_vhh = bool(h) and not (vl and vh)
            if not is_vhh and (not vl or not vh):
                continue
            row = {
                "design_id": f"{meta['run_id'][:15]}-{rec.get('example_id', '')}",
                **meta,
                "kind": "vhh" if is_vhh else "scfv",
                "vh_sequence": h if is_vhh else vh,
                "vl_sequence": "" if is_vhh else vl,
                **(evaluate_vhh(h) if is_vhh else evaluate(vl, vh, novelty)),
                "hotspot_obedience": rec.get("hotspot_obedience.hotspot_obedience", ""),
                "pre_oracle_pass": rec.get("pre_oracle_pass", ""),
                "sequence_adherence": rec.get("sequence_adherence.sequence_adherence", ""),
            }
            ok, why = verdict_vhh(row) if is_vhh else verdict(row)
            row["passes_local_filter"] = int(ok)
            row["filter_reason"] = why
            rows.append(row)

    if not rows:
        raise SystemExit("no designs found in the runs")

    # identical constructs are counted once; the first is kept
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["vh_sequence"], row["vl_sequence"])
        if key in seen and row["passes_local_filter"]:
            row["passes_local_filter"] = 0
            row["filter_reason"] = "duplicate of an earlier design"
        seen.add(key)

    args.out.mkdir(parents=True, exist_ok=True)
    allfile = args.out / "candidates_all.csv"
    passfile = args.out / "candidates_pass.csv"
    fields = list(rows[0])
    with allfile.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    kept = [r for r in rows if r["passes_local_filter"]]
    with passfile.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(kept)

    print(f"{len(rows)} designs from {len(run_dirs)} runs -> {allfile}")
    print(f"{len(kept)} pass the filtering -> {passfile}")
    reasons: dict[str, int] = {}
    for r in rows:
        if not r["passes_local_filter"]:
            reasons[r["filter_reason"]] = reasons.get(r["filter_reason"], 0) + 1
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {why}")
    for run_dir in run_dirs:
        n = sum(1 for r in rows if r["run_id"] == run_dir.name)
        k = sum(1 for r in kept if r["run_id"] == run_dir.name)
        print(f"  {run_dir.name}: {k}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
