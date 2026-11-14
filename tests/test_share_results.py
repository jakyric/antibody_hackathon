"""Tests for what gets released into the public repository.

This is the one test where a failure doesn't just give a wrong number, but puts
our candidate sequences out for all the other teams. It checks the rule, not
that someone remembered to remove a column.
"""

import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sr = load("share_results", "scripts/share_results.py")


def test_sequence_columns_are_dropped():
    rows = [{"id": "d1", "cdrh3": "APGDGFQIMD", "iptm": "0.93"}]
    out, dropped = sr.strip_sequences(rows, ["id", "cdrh3", "iptm"])
    assert out == [{"id": "d1", "iptm": "0.93"}]
    assert dropped == ["cdrh3"]


def test_whole_chains_are_dropped_too():
    fields = ["design_id", "vh_sequence", "vl_sequence", "scfv_sequence",
              "heavy_sequence", "light_sequence", "score"]
    out, dropped = sr.strip_sequences([{f: "x" for f in fields}], fields)
    assert list(out[0]) == ["design_id", "score"]
    assert len(dropped) == 5


def test_uppercase_does_not_help():
    out, dropped = sr.strip_sequences([{"CDRH3": "A", "n": "1"}], ["CDRH3", "n"])
    assert dropped == ["CDRH3"] and list(out[0]) == ["n"]


def test_id_is_kept():
    """A design ID without a sequence is verifiable by the team and useless to others."""
    out, _ = sr.strip_sequences([{"id": "tz-opt-05", "iptm": "0.9"}], ["id", "iptm"])
    assert out[0]["id"] == "tz-opt-05"


def test_no_shared_table_contains_a_sequence():
    """The check on real data: every file we actually share, column by column.

    An amino acid sequence is long and consists only of amino acid letters. If
    something like that shows up in a shared table, a column has slipped through.
    """
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    for src_rel, _ in sr.TABLES:
        src = ROOT / "results" / src_rel
        if not src.exists():
            continue
        with src.open(newline="") as fh:
            reader = csv.DictReader(fh)
            rows, fields = list(reader), list(reader.fieldnames or [])
        kept, _ = sr.strip_sequences(rows, fields)
        for row in kept:
            for key, value in row.items():
                v = str(value).strip().upper()
                assert not (len(v) >= 8 and set(v) <= aa), \
                    f"{src_rel}: column {key!r} looks like a sequence: {v[:20]}"
