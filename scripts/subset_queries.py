#!/usr/bin/env python3
"""Keep only the queries that belong to designs we actually intend to fold.

refold-prepare picks the top N by the design model's own scores and knows nothing
about the gate. We only fold what has passed. This filtering is the reason
no GPU time is spent on designs that could not be delivered anyway.

    python scripts/subset_queries.py --queries prepared/holo/queries.json \
        --keep results/denovo/candidates_pass.csv --out prepared/holo/queries_gated.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queries", type=Path, required=True)
    ap.add_argument("--keep", type=Path, required=True,
                    help="CSV with a design_id column")
    ap.add_argument("--id-column", default="design_id")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    with args.keep.open(newline="") as fh:
        wanted = {r[args.id_column] for r in csv.DictReader(fh)}
    # our design_ids are "<run-id>-<example_id>"; queries use example_id
    suffixes = {w.split("-", 1)[-1] for w in wanted} | wanted

    data = json.loads(args.queries.read_text())
    keep = [q for q in data["queries"]
            if q["query_id"] in wanted
            or any(q["query_id"] == s or q["query_id"].startswith(s) for s in suffixes)]
    if not keep:
        raise SystemExit("none of the queries matched the list; check the ids before folding")

    data["queries"] = keep
    args.out.write_text(json.dumps(data))
    print(f"{len(keep)} of {len(json.loads(args.queries.read_text())['queries'])} "
          f"queries kept -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
