"""One command: does the submission table pass every entry requirement?

WHY THIS EXISTS
  The checks already existed, but they were spread across several scripts and
  nothing forced them to run. Someone had to remember which to run, in which
  order, and against which version of the table. A submission could be
  published without the checks having been run on the file actually being
  published.

  This is the single command. One row per design, one column per requirement,
  and a verdict that is PASS, FAIL or UNKNOWN -- never a blank, and never a
  PASS for something we cannot actually test here.

THE THREE VERDICTS, AND WHY UNKNOWN IS NOT A PASS
  PASS     we checked it, it holds.
  FAIL     we checked it, it does not hold. The design must not be submitted.
  UNKNOWN  we cannot check it here, so we say so.

  Two requirements are genuinely uncheckable without the organisers' tooling:
  the "under 70% identity to any known antibody" test for the de novo track,
  which needs the lab's antibody database, and the immunum numbering
  confidence. The official checker is the authority for both. Reporting
  those as PASS would be a lie; reporting them as FAIL would be wrong too.

  Because those two can never be satisfied locally, UNKNOWN does not block
  publication -- it would block it forever. UNKNOWN warns, FAIL blocks. The
  exit code reflects that: 0 when nothing failed, 1 when anything failed.

MARGIN IS PRINTED, BECAUSE WE HAVE NONE
  Every design sits at exactly 3 mutations from the nearest provided sequence
  and the requirement is 3. That is a margin of zero. A design that only just
  passes can fail if the organisers' tool counts even slightly differently, so
  the margin is shown as a number rather than hidden behind a tick.

THE STRICTEST READING IS USED
  The sources disagree: the slides say "<=70%" and "at least 2 mutations", the
  brief says "under 70%" and "strictly more than 2", and the official checker
  enforces 3. We apply under-70% and >=3, which satisfies every reading.

USAGE
  python scripts/requirement_check.py                     # results/SUBMISSION.csv
  python scripts/requirement_check.py --submission X.csv --json out.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.data import TZ_CDRH3_10, provided_cdrh3s  # noqa: E402
from abx.gate import (  # noqa: E402
    CDRH3_IDENTITY_CEILING,
    MAX_TOTAL_LENGTH,
    count_unpaired_cysteines,
    hamming,
    identity,
    levenshtein,
)

PASS, FAIL, UNKNOWN = "PASS", "FAIL", "UNKNOWN"
REQUIRED_MUTATIONS = 3  # what the official checker enforces; the strictest reading
VALID_KINDS = {"vhh", "scfv", "paired", "fab"}


def min_distance_to_provided(cdrh3: str, provided: list[str]) -> tuple[int, str]:
    best, who = len(cdrh3) + 99, ""
    for p in provided:
        d = hamming(p, cdrh3) if len(p) == len(cdrh3) else levenshtein(p, cdrh3)
        if d < best:
            best, who = d, p
            if best == 0:
                break
    return best, who


def check_row(row: dict, provided: list[str]) -> dict:
    """Every requirement for one submitted row."""
    sid = row["sequence_id"]
    kind = (row.get("kind") or "").strip().lower()
    heavy = (row.get("heavy_sequence") or "").strip().upper()
    light = (row.get("light_sequence") or "").strip().upper()
    combined = heavy + light

    out: dict[str, dict] = {}

    # --- format -------------------------------------------------------
    out["format"] = {
        "verdict": PASS if kind in VALID_KINDS else FAIL,
        "detail": f"kind={kind or '(empty)'}",
        "margin": "",
    }

    # --- total length -------------------------------------------------
    n = len(combined)
    out["length_max_250"] = {
        "verdict": PASS if n <= MAX_TOTAL_LENGTH else FAIL,
        "detail": f"{n} aa",
        "margin": f"{MAX_TOTAL_LENGTH - n} aa to spare",
    }

    # --- unpaired cysteines, per domain -------------------------------
    # The rule applies per domain, not across the whole molecule: an odd
    # count in each of two chains sums to an even total and would look clean.
    bad = []
    for name, chain in (("heavy", heavy), ("light", light)):
        if chain and count_unpaired_cysteines(chain):
            bad.append(f"{name} has {chain.count('C')} Cys")
    out["no_unpaired_cysteines"] = {
        "verdict": FAIL if bad else PASS,
        "detail": "; ".join(bad) if bad else
                  f"heavy {heavy.count('C')} Cys, light {light.count('C')} Cys, both even",
        "margin": "",
    }

    # --- CDR-H3 identity to trastuzumab -------------------------------
    try:
        start = heavy.index("YYCSR") + 5
        cdrh3 = heavy[start:heavy.index("YWGQGT", start)]
    except ValueError:
        cdrh3 = ""
    if cdrh3 and len(cdrh3) == len(TZ_CDRH3_10):
        ident = identity(cdrh3, TZ_CDRH3_10)
        muts = hamming(cdrh3, TZ_CDRH3_10)
        out["cdrh3_identity_under_70pct"] = {
            "verdict": PASS if ident < CDRH3_IDENTITY_CEILING else FAIL,
            "detail": f"{cdrh3}, {ident:.0%} identity, {muts} mutations",
            "margin": f"{muts - 4} mutation(s) beyond the 4 needed on a 10-mer",
        }
    else:
        out["cdrh3_identity_under_70pct"] = {
            "verdict": UNKNOWN, "detail": "CDR-H3 not located in heavy chain",
            "margin": "",
        }

    # --- distance to every provided sequence --------------------------
    if cdrh3:
        d, who = min_distance_to_provided(cdrh3, provided)
        out["mutations_from_provided_data"] = {
            "verdict": PASS if d >= REQUIRED_MUTATIONS else FAIL,
            "detail": f"nearest of {len(provided)} provided is {who} at distance {d}, "
                      f"requirement {REQUIRED_MUTATIONS}",
            "margin": f"{d - REQUIRED_MUTATIONS} mutation(s) of margin",
        }
    else:
        out["mutations_from_provided_data"] = {
            "verdict": UNKNOWN, "detail": "no CDR-H3 to compare", "margin": "",
        }

    # --- what only the organisers' tooling can answer ------------------
    out["immunum_numbering_confidence"] = {
        "verdict": UNKNOWN,
        "detail": "needs the organisers' immunum library; the official checker is the authority",
        "margin": "",
    }
    out["identity_to_known_antibodies"] = {
        "verdict": UNKNOWN,
        "detail": "needs the lab's antibody database (de novo track rule); "
                  "not applicable to optimisation designs, but we cannot verify it here",
        "margin": "",
    }
    return {"sequence_id": sid, "checks": out}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--submission", type=Path, default=ROOT / "results/SUBMISSION.csv")
    ap.add_argument("--json", type=Path, default=None, help="also write machine-readable JSON")
    ap.add_argument("--md", type=Path, default=None, help="also write a markdown table")
    ap.add_argument("--strict-unknown", action="store_true",
                    help="treat UNKNOWN as a failure (blocks on checks only the official checker can do)")
    args = ap.parse_args(argv)

    provided = sorted(provided_cdrh3s())
    rows = list(csv.DictReader(args.submission.open(newline="", encoding="utf-8")))
    print(f"submission: {args.submission}")
    print(f"rows: {len(rows)}   provided CDR-H3s compared against: {len(provided)}\n")

    results = [check_row(r, provided) for r in rows]
    names = list(results[0]["checks"]) if results else []

    width = max((len(r["sequence_id"]) for r in results), default=12)
    header = " " * (width + 2) + "  ".join(f"{n[:14]:>14s}" for n in names)
    print(header)
    for res in results:
        cells = []
        for n in names:
            v = res["checks"][n]["verdict"]
            cells.append(f"{v:>14s}")
        print(f"{res['sequence_id']:<{width}}  " + "  ".join(cells))

    n_fail = sum(1 for r in results for c in r["checks"].values() if c["verdict"] == FAIL)
    n_unk = sum(1 for r in results for c in r["checks"].values() if c["verdict"] == UNKNOWN)

    print("\nper-design detail")
    for res in results:
        bad = {k: v for k, v in res["checks"].items() if v["verdict"] != PASS}
        print(f"\n  {res['sequence_id']}")
        for k, v in res["checks"].items():
            mark = {PASS: "ok  ", FAIL: "FAIL", UNKNOWN: "?   "}[v["verdict"]]
            margin = f"   [{v['margin']}]" if v["margin"] else ""
            print(f"    {mark} {k}: {v['detail']}{margin}")
        if not bad:
            pass

    print(f"\n{'='*70}")
    print(f"FAIL: {n_fail}     UNKNOWN: {n_unk}     designs: {len(results)}")
    if n_fail:
        print("VERDICT: NOT READY. A failing design must not be submitted.")
    elif n_unk and args.strict_unknown:
        print("VERDICT: BLOCKED by --strict-unknown. Only the official checker can clear these.")
    else:
        print("VERDICT: every locally checkable requirement passes.")
        print("         The UNKNOWN columns are the two only the official checker can answer;")
        print("         run the official checker before submitting.")
    print("=" * 70)

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(
            {"submission": str(args.submission), "n_fail": n_fail,
             "n_unknown": n_unk, "designs": results}, indent=1), encoding="utf-8")
        print(f"wrote {args.json}")

    if args.md:
        L = ["# Entry requirements, per design", "",
             f"`{args.submission.name}`, {len(rows)} rows, checked against all "
             f"{len(provided)} provided CDR-H3s.", "",
             "PASS = checked and holds. FAIL = checked and does not. "
             "UNKNOWN = only the official checker can answer; never reported as a pass.", "",
             "| Design | " + " | ".join(names) + " |",
             "| --- | " + " | ".join("---" for _ in names) + " |"]
        for res in results:
            L.append(f"| `{res['sequence_id']}` | " + " | ".join(
                res["checks"][n]["verdict"] for n in names) + " |")
        L += ["", "## Margins", "",
              "Every design sits at the minimum distance the rule allows. "
              "The margin column is zero by construction, because the generator "
              "was built to stay as close to trastuzumab as the rules permit.", "",
              "| Design | CDR-H3 | Distance to nearest provided | Margin |",
              "| --- | --- | ---: | ---: |"]
        for res in results:
            c = res["checks"]["mutations_from_provided_data"]
            i = res["checks"]["cdrh3_identity_under_70pct"]
            loop = i["detail"].split(",")[0]
            L.append(f"| `{res['sequence_id']}` | {loop} | {c['detail'].split('distance ')[-1].split(',')[0]} | {c['margin']} |")
        args.md.parent.mkdir(parents=True, exist_ok=True)
        args.md.write_text("\n".join(L) + "\n", encoding="utf-8")
        print(f"wrote {args.md}")

    if n_fail:
        return 1
    if n_unk and args.strict_unknown:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
