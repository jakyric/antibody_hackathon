"""S2 rule harness: hard rules + novelty/identity gates for every design.

    .venv/bin/python denovo/s2_rules/rules.py                       # all results/denovo/designs_*.csv -> rules.csv
    .venv/bin/python denovo/s2_rules/rules.py --in a.csv b.csv --out x.csv [--no-vm]

Columns (one row per design_id), grouped:
  local hard rules : len_ok, invalid_aa, cys_count, noncanon_cys, unpaired_cys (abguard definition: odd count),
                     canon_cys_ok (Cys at IMGT 23 and 104), numbering_ok, imm_chain, imm_conf, cdrh3, cdrh3_len,
                     hallmarks (IMGT 42/49/50/52), cdrh3_dist_trast (Levenshtein, min over IMGT 13-mer and 10-mer core)
  abguard replica  : identity_to_trastuzumab (abguard formula), novelty_class_local, cdrh3_min_dist_ref (+ nearest
                     source/id) against the 610k CDRH3s bundled in abguard, cdrh3_novel_local
  known panel      : identity_nearest_known (whole Fv, PDB heavy domains + human IGHV germlines), nearest_known_id,
                     nearest_known_source, identity_nearest_vhh, ident_lt70, ident_lt68
  real abguard (VM): abguard_pass (is_valid with its bundled reference), abguard_errors, abguard_class,
                     abguard_identity_trast, abguard_cdrh3_min_dist, abguard_nearest, abguard_n_high_liab, abguard_liab_tags
  verdicts         : basic_pass (local hard rules), denovo_vm_pass (abguard_pass AND class de_novo),
                     denovo_strict_pass (denovo_vm_pass AND ident_lt70), gap_vm_ok_ident_ge70 (the gap to watch)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_lib as L  # noqa: E402

dn = L.dn
AAS = set("ACDEFGHIKLMNPQRSTVWY")


def local_rules(seq: str) -> dict:
    a = L.analyse(seq)
    num = a.get("numbering", {})
    cys_idx = [lab for lab, aa in num.items() if aa == "C"]
    cys_total = seq.count("C")
    noncanon = [lab for lab in cys_idx if lab not in ("23", "104")]
    # second-disulfide tolerance: an even number of non-canonical Cys (both partners present) is tolerated, but flagged
    unpaired = cys_total % 2 == 1
    hall = "".join(num.get(p, "-") for p in L.HALLMARK_POS)
    d_imgt = L.Levenshtein.distance(a.get("cdr3", ""), L.TRAS_CDR3)
    d_core = L.Levenshtein.distance(a.get("cdr3", ""), L.TRAS_CDRH3_CORE)
    return dict(
        vhh_len=len(seq), len_ok=int(len(seq) <= 250), invalid_aa=int(bool(set(seq) - AAS)),
        cys_count=cys_total, noncanon_cys=len(noncanon), unpaired_cys=int(unpaired),
        canon_cys_ok=int("23" in cys_idx and "104" in cys_idx),
        numbering_ok=int(a["chain"] == "H" and not a["error"] and a["conf"] >= 0.1), imm_chain=a["chain"],
        imm_conf=round(a["conf"], 3), cdrh3=a.get("cdr3", ""), cdrh3_len=len(a.get("cdr3", "")), hallmarks=hall,
        cdrh3_dist_trast=min(d_imgt, d_core), cdrh3_dist_trast_imgt=d_imgt,
        identity_to_trastuzumab=round(L.lev_identity(a["fv"], L.TRAS_FV), 4) if a.get("fv") else "",
        fv=a.get("fv", ""),
    )


def run(designs: list[dict], use_vm: bool = True, panel: L.Panel | None = None) -> list[dict]:
    panel = panel or L.Panel()
    rows = []
    for d in designs:
        r = dict(design_id=d["design_id"], parent_id=d.get("parent_id", ""), generator=d.get("generator", ""),
                 is_control=int(d.get("generator") == "control" or d["design_id"].startswith("ctl-")),
                 vhh_seq=d["vhh_seq"])
        r.update(local_rules(d["vhh_seq"]))
        rows.append(r)
    ok = [r for r in rows if r["fv"]]
    # known-antibody panel (whole Fv). A design never counts as "known" because of its own parent: we do NOT exclude it,
    # a diversified child must be far from its parent too.
    near = panel.nearest([r["fv"] for r in ok])
    for r, nn in zip(ok, near):
        r.update({k: (round(v, 4) if isinstance(v, float) else v) for k, v in nn.items()})
        r["ident_lt70"] = int(nn["identity_nearest_known"] < 0.70)
        r["ident_lt68"] = int(nn["identity_nearest_known"] < 0.68)
        r["novelty_class_local"] = "novel" if float(r["identity_to_trastuzumab"]) >= 0.70 else "de_novo"
    # abguard's CDRH3 prior-art gate, local replica (610k CDRH3s from FACS/SPR tables, PDB, patents, plabdab, therasabdab...)
    md = L.cdrh3_min_dist([r["cdrh3"] for r in ok])
    for r, (dist, ref, src, rid) in zip(ok, md):
        r.update(cdrh3_min_dist_ref=dist, cdrh3_nearest_ref=ref, cdrh3_nearest_ref_src=f"{src}:{rid}" if rid not in ("", "None", "nan") else src)
        need = 3 if r["novelty_class_local"] == "novel" else 2
        r["cdrh3_novel_local"] = int(dist >= need and r["cdrh3_dist_trast_imgt"] >= need)
    # real abguard on the VM
    if use_vm and ok:
        try:
            res = L.abguard_vm({r["design_id"]: r["vhh_seq"] for r in rows})
            for r in rows:
                g = res[r["design_id"]]
                r.update(abguard_pass=int(bool(g.get("is_valid"))), abguard_errors=";".join(g.get("errors") or []),
                         abguard_class=g.get("novelty_class") or "", abguard_identity_trast=round(g.get("identity_to_trastuzumab") or 0, 4),
                         abguard_cdrh3_min_dist=g.get("cdrh3_min_distance"), abguard_nearest=f'{g.get("cdrh3_nearest_source")}:{g.get("cdrh3_nearest_id")}',
                         abguard_n_high_liab=g.get("n_high_severity"), abguard_liab_tags=";".join(g.get("liability_tags") or []))
        except Exception as e:  # keep going: local columns are still useful
            print("WARNING: abguard on the VM failed:", str(e)[:300], file=sys.stderr)
    for r in rows:
        basic = (r["len_ok"] and not r["invalid_aa"] and not r["unpaired_cys"] and r.get("canon_cys_ok") and r["numbering_ok"]
                 and r.get("cdrh3_novel_local", 0))
        r["basic_pass"] = int(bool(basic))
        if "abguard_pass" in r:
            r["denovo_vm_pass"] = int(r["abguard_pass"] == 1 and r["abguard_class"] == "de_novo")
            r["denovo_strict_pass"] = int(r["denovo_vm_pass"] == 1 and r.get("ident_lt70", 0) == 1)
            r["gap_vm_ok_ident_ge70"] = int(r["denovo_vm_pass"] == 1 and r.get("ident_lt70", 0) == 0)
    return rows


COLS = ["design_id", "parent_id", "generator", "is_control", "vhh_len", "len_ok", "invalid_aa", "cys_count", "noncanon_cys",
        "unpaired_cys", "canon_cys_ok", "numbering_ok", "imm_chain", "imm_conf", "cdrh3", "cdrh3_len", "hallmarks",
        "cdrh3_dist_trast", "cdrh3_dist_trast_imgt", "identity_to_trastuzumab", "novelty_class_local", "cdrh3_min_dist_ref",
        "cdrh3_nearest_ref", "cdrh3_nearest_ref_src", "cdrh3_novel_local", "identity_nearest_known", "nearest_known_id",
        "nearest_known_source", "identity_nearest_vhh", "nearest_vhh_id", "ident_lt70", "ident_lt68", "abguard_pass",
        "abguard_errors", "abguard_class", "abguard_identity_trast", "abguard_cdrh3_min_dist", "abguard_nearest",
        "abguard_n_high_liab", "abguard_liab_tags", "basic_pass", "denovo_vm_pass", "denovo_strict_pass",
        "gap_vm_ok_ident_ge70", "vhh_seq"]


BOOL_COLS = ["is_control", "len_ok", "invalid_aa", "unpaired_cys", "canon_cys_ok", "numbering_ok", "cdrh3_novel_local", "ident_lt70",
             "ident_lt68", "abguard_pass", "basic_pass", "denovo_vm_pass", "denovo_strict_pass", "gap_vm_ok_ident_ge70"]
NUM_COLS = ["identity_nearest_known", "identity_nearest_vhh", "identity_to_trastuzumab", "abguard_identity_trast", "imm_conf"]


def to_csv_rows(rows: list[dict]) -> list[dict]:
    """booleans as True/False (pandas reads them as bool), identity columns always numeric (1.0 = unnumberable, worst case)."""
    out = []
    for r in rows:
        r = dict(r)
        for k in BOOL_COLS:
            v = r.get(k, "")
            r[k] = "" if v in ("", None) else ("True" if str(v) in ("1", "True", "true") else "False")
        if r.get("identity_nearest_known", "") in ("", None):
            r["identity_nearest_known"] = 1.0
        out.append(r)
    return out


def from_csv_rows(rows: list[dict]) -> list[dict]:
    for r in rows:
        for k in BOOL_COLS:
            if r.get(k) in ("True", "False"):
                r[k] = int(r[k] == "True")
    return rows


def cycle(out_path: str, use_vm: bool, panel, log=print) -> int:
    """One incremental pass: compute rules only for new/changed designs (or ones whose abguard result is still missing)."""
    path = Path(out_path)
    have = {r["design_id"]: r for r in from_csv_rows(dn.read_csv(path))} if path.exists() else {}
    designs = [d for d in dn.load_designs() if d.get("vhh_seq")]
    todo = [d for d in designs if d["design_id"] not in have or have[d["design_id"]]["vhh_seq"] != d["vhh_seq"]
            or (use_vm and have[d["design_id"]].get("abguard_pass") in ("", None))]
    current = {d["design_id"] for d in designs}
    stale = [k for k in have if k not in current]  # design tables can lose rows (e.g. parked files): keep rules.csv in step
    if not todo and not stale:
        return 0
    new = run(todo, use_vm=use_vm, panel=panel) if todo else []
    for r in new:
        have[r["design_id"]] = r
    rows = [have[d["design_id"]] for d in designs if d["design_id"] in have]
    dn.write_csv(path, to_csv_rows(rows), COLS)
    log(f"{dn.utcnow()} rules: +{len(new)} designs, -{len(stale)} stale, {len(rows)} total -> {path}")
    return len(new)


def watch(a):
    import time
    panel = L.Panel()
    log = lambda m: print('[rules] ' + str(m), flush=True)
    log(f"rules watch start pid={__import__('os').getpid()} interval={a.interval}s until={a.until} out={a.out}")
    while True:
        if L.stop_requested():
            log("STOP file found, exiting"); break
        if L.past_deadline(a.until):
            log(f"deadline {a.until} reached, exiting"); break
        try:
            cycle(a.out, not a.no_vm, panel, log)
        except Exception as e:  # never die in a loop: log and retry next cycle
            log(f"{dn.utcnow()} rules cycle error: {type(e).__name__}: {str(e)[:300]}")
        for _ in range(a.interval):
            if L.stop_requested() or L.past_deadline(a.until):
                break
            time.sleep(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--watch", action="store_true", help="loop: process only new/changed designs every --interval s")
    ap.add_argument("--interval", type=int, default=120)
    ap.add_argument("--until", default="15:45", help="local (CET) HH:MM at which the watcher exits")
    ap.add_argument("--in", dest="inp", nargs="*", default=None, help="designs csv(s); default all results/denovo/designs_*.csv")
    ap.add_argument("--out", default=str(dn.OUT / "rules.csv"))
    ap.add_argument("--no-vm", action="store_true", help="skip the real abguard (ssh to the VM)")
    a = ap.parse_args()
    if not __import__('os').environ.get('PIPELINE_VM_HOST'):
        a.no_vm = True  # no VM configured: local replica only
    if a.watch:
        return watch(a)
    designs = []
    if a.inp:
        for f in a.inp:
            designs += dn.read_csv(f)
    else:
        designs = dn.load_designs()
    designs = [d for d in designs if d.get("vhh_seq")]
    rows = run(designs, use_vm=not a.no_vm)
    dn.write_csv(a.out, to_csv_rows(rows), COLS)
    n = len(rows)
    print(f"{n} designs -> {a.out}")
    for k in ("basic_pass", "denovo_vm_pass", "denovo_strict_pass", "gap_vm_ok_ident_ge70", "ident_lt70", "abguard_pass"):
        vals = [r.get(k) for r in rows if r.get(k) not in (None, "")]
        if vals:
            print(f"  {k}: {sum(int(v) for v in vals)}/{len(vals)}")


if __name__ == "__main__":
    main()
