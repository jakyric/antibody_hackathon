#!/usr/bin/env python3
"""S4 orchestrator: yield / Tm / developability profile for every VHH design, incremental and resumable.

    python3 denovo/s4_props/props.py                      # one pass over results/denovo/designs_*.csv + baselines
    python3 denovo/s4_props/props.py --watch 120          # keep polling for new designs every 120 s
    python3 denovo/s4_props/props.py --where vm           # everything on the shared VM (one remote call per tool)
    python3 denovo/s4_props/props.py --in-csv x.csv --no-designs --out /tmp/p.csv    # ad-hoc sequences (design_id,vhh_seq)

Tools and columns (all keyed by design_id; blank = not computed):
  nanomelt  nm_tm                                    apparent Tm (C), NanoMelt ensemble, local env `nanomelt` / VM env `s4`
  tnp       tnp_cdr_len tnp_cdr3_len tnp_cdr3_compact tnp_psh tnp_ppc tnp_pnc tnp_tetrad tnp_flag_* tnp_n_red tnp_n_amber
  liab      liab_* (motifs, pI, charge, CDR hydrophobicity) and liab_f_* flags, liab_n_flags
  nb9       nb9_score nb9_applicable nb6_score nb6_applicable nb9_flags nb_has_cys_met   (monoyield-paper CDR2/CDR3 rules; see nanobody_rules.py)
  monoyield my_yield my_yield_sd (log10 mg/L monomer, E. coli) my_mdratio my_mdratio_sd (log10 monomer/dimer), my_dtype
--where   (VM paths below are DISABLED unless --allow-vm is passed; the user must grant it)
          auto (default): liab, nanomelt, tnp local with VM fallback; monoyield VM (needs ESMC-6B, 25 GB)
          local: no VM at all (monoyield skipped)      vm: all four tools on the VM
          both : as auto, plus the VM also computes nanomelt/tnp/liab into *_vm columns (cross-machine regression check)
--max-tnp N  only run TNP on N pending designs (baselines first, then highest nm_tm); TNP is cheap here (about 3 s per
          sequence on a GPU) but the flag bounds the worst case.
Cost tiers (--tier on, default): liab + nanomelt run on EVERY design; tnp + monoyield only on the "candidate set" =
controls + baselines + designs that pass S2's gate (`denovo_vm_pass` in rules.csv) AND are in the top --top-frac (0.35) of the
gate-passers by `of3_consensus_rank_score` (results/denovo/score_of3.csv, higher = better).  No scores yet -> controls only.
generator=fwdiv rows are ignored (framework diversification was cancelled).
--watch N loops every min(N, 300) s until --stop-at (default 15:45 local) or the file denovo/STOP exists; monoyield is
called at most once per --my-interval s (600) as ONE batched VM call.  Writes are atomic (tmp + rename).
Re-running never recomputes a finished (design_id, vhh_seq) pair; a changed sequence under the same id is recomputed.
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
import dn_common as dn  # noqa: E402

MM = shutil.which("micromamba") or str(Path.home() / ".local/bin/micromamba")
TOOLS = ["liab", "nb9", "nanomelt", "tnp", "monoyield"]
DONE_COL = {"liab": "liab_len", "nb9": "nb9_applicable", "nanomelt": "nm_tm", "tnp": "tnp_cdr_len", "monoyield": "my_yield"}
ERR_COL = {"liab": None, "nb9": None, "nanomelt": "nm_error", "tnp": "tnp_error", "monoyield": None}
TOOL_COLS = {
    "liab": None,  # filled from the tool output (liab_*)
    "nb9": None,   # nb9_score nb9_applicable nb6_score nb6_applicable nb9_flags nb_has_cys_met nb_cdr2_core nb_cdr3_loop
    "nanomelt": ["nm_tm", "nm_error"],
    "tnp": ["tnp_cdr_len", "tnp_cdr3_len", "tnp_cdr3_compact", "tnp_psh", "tnp_ppc", "tnp_pnc", "tnp_tetrad",
            "tnp_flag_l", "tnp_flag_l3", "tnp_flag_c", "tnp_flag_psh", "tnp_flag_ppc", "tnp_flag_pnc",
            "tnp_n_red", "tnp_n_amber", "tnp_error"],
    "monoyield": ["my_yield", "my_yield_sd", "my_mdratio", "my_mdratio_sd", "my_dtype"],
}
ALLOW_VM = [False]
TNP_ENV: list = []  # ["env", "CUDA_VISIBLE_DEVICES="] with --tnp-cpu
IGNORE_GENERATORS = {"fwdiv", "boltzgen", "boltzgen-cys"}
LOG = lambda m: print(f"[props {time.strftime('%H:%M:%S')}] {m}", flush=True)  # noqa: E731


# ------------------------------------------------------------------ candidates and state
def candidates(a) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    if not a.no_baselines:
        for r in dn.read_csv(HERE / "baselines.csv"):
            rows[r["design_id"]] = {"design_id": r["design_id"], "vhh_seq": r["vhh_seq"], "generator": "baseline",
                                    "framework": r["name"], "is_control": 1}
    if not a.no_designs:
        for r in dn.load_designs():
            if not r.get("vhh_seq"):  # e.g. designs_boltzgen_metrics.csv also matches designs_*.csv
                continue
            rows.setdefault(r["design_id"], {
                "design_id": r["design_id"], "vhh_seq": r["vhh_seq"], "generator": r.get("generator", ""),
                "framework": r.get("framework", ""), "parent_id": r.get("parent_id", ""),
                "is_control": int(r.get("generator") == "control" or r["design_id"].startswith("ctl-")),
                "cdr1": r.get("cdr1", ""), "cdr2": r.get("cdr2", ""), "cdr3": r.get("cdr3", "")})
    for f in a.in_csv or []:
        for r in dn.read_csv(f):
            rows.setdefault(r["design_id"], {"design_id": r["design_id"], "vhh_seq": r["vhh_seq"],
                                             "generator": r.get("generator", "adhoc"), "framework": r.get("framework", ""),
                                             "is_control": int(r["design_id"].startswith("ctl-"))})
    ids = set(a.ids) if a.ids else None
    out = {}
    seen_seq: dict[str, str] = {}
    for k, r in rows.items():
        if ids and k not in ids:
            continue
        if r.get("generator") in IGNORE_GENERATORS:  # cancelled generators: fwdiv, BoltzGen (user stopped it); old rows stay, ignored
            continue
        r["vhh_seq"] = r["vhh_seq"].strip().upper()
        if not r["is_control"]:  # identical sequences under different ids (re-used seeds): compute once, first id wins
            if r["vhh_seq"] in seen_seq:
                continue
            seen_seq[r["vhh_seq"]] = k
        out[k] = r
    return out


def load_state(path: Path) -> dict[str, dict]:
    return {r["design_id"]: r for r in dn.read_csv(path)} if path.exists() else {}


IDENTITY_COLS = {"design_id", "generator", "framework", "is_control", "stage_ts", "tools_done", "vhh_seq"}


def reconcile(state: dict, cands: dict) -> int:
    """Rows whose design_id left the designs table (S1a de-duplicated re-used seeds) are removed from props.csv; any finished tool
    value is first moved to the id that now carries the same sequence, so no TNP/monoyield result is lost.  Rows of ignored
    generators (BoltzGen, fwdiv) are left alone."""
    by_seq = {c["vhh_seq"]: k for k, c in cands.items() if not c["is_control"]}
    moved = dropped = 0
    for did in [d for d, st in state.items() if d not in cands and st.get("generator") not in IGNORE_GENERATORS
                and st.get("generator") != "baseline"]:
        st = state.pop(did)
        keep = by_seq.get(st.get("vhh_seq", ""))
        if keep:
            dst = state.setdefault(keep, {"design_id": keep, "vhh_seq": st["vhh_seq"], "generator": cands[keep]["generator"],
                                          "framework": cands[keep]["framework"], "is_control": cands[keep]["is_control"],
                                          "stage_ts": dn.utcnow()})
            for k, v in st.items():
                if k not in IDENTITY_COLS and v != "" and dst.get(k, "") == "":
                    dst[k] = v
                    moved += 1
        dropped += 1
    return dropped


def pending(tool: str, cands: dict, state: dict, retry_failed: bool) -> list[dict]:
    done_col, err_col = DONE_COL[tool], ERR_COL[tool]
    out = []
    for did, c in cands.items():
        s = state.get(did)
        if s and s.get("vhh_seq") == c["vhh_seq"]:
            if s.get(done_col, "") != "":
                continue
            if err_col and s.get(err_col, "") != "" and not retry_failed:
                continue
        if dn.check_seq(c["vhh_seq"]) and not c["is_control"]:
            continue  # nonsense sequence: S2 rules reject it, do not waste compute
        out.append(c)
    return out


def write_in_csv(rows: list[dict], path: Path) -> None:
    dn.write_csv(path, [{"design_id": r["design_id"], "vhh_seq": r["vhh_seq"], "cdr1": r.get("cdr1", ""),
                         "cdr2": r.get("cdr2", ""), "cdr3": r.get("cdr3", "")} for r in rows],
                 ["design_id", "vhh_seq", "cdr1", "cdr2", "cdr3"])


def expensive_ids(cands: dict, a) -> set:
    """design_ids that earn the slow tools (tnp, monoyield): controls + baselines, plus designs that pass S2's gate
    (`denovo_vm_pass` in rules.csv, when that file exists) AND sit in the top --top-frac of the gate-passers by S3's
    `of3_consensus_rank_score` (score_of3.csv, higher = better).  No scores yet -> controls only."""
    keep = {k for k, c in cands.items() if c["is_control"]}
    f = dn.OUT / "score_of3.csv"
    if not f.exists():
        return keep
    gate = None
    rf = dn.OUT / "rules.csv"
    if rf.exists():
        rr = dn.read_csv(rf)
        if rr and a.gate_col in rr[0]:
            gate = {r["design_id"] for r in rr if str(r[a.gate_col]).strip().lower() in ("true", "1", "1.0")}
    sc = {}
    for r in dn.read_csv(f):
        try:
            sc[r["design_id"]] = float(r[a.score_col])
        except (KeyError, ValueError, TypeError):
            continue
    real = {k: v for k, v in sc.items() if k in cands and not cands[k]["is_control"] and v == v
            and (gate is None or k in gate)}
    if not real:
        return keep
    vals = sorted(real.values(), reverse=a.score_higher_better)
    n_top = max(1, int(round(a.top_frac * len(vals))))
    cut = vals[n_top - 1]
    ok = (lambda v: v >= cut) if a.score_higher_better else (lambda v: v <= cut)
    return keep | {k for k, v in real.items() if ok(v)}


# ------------------------------------------------------------------ runners; each returns {design_id: {col: value}}
def mm_run(env: str, *cmd: str, log: Path | None = None) -> None:
    full = [MM, "run", "-n", env, *cmd]
    with open(log or os.devnull, "w") as fh:
        r = subprocess.run(full, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError(f"{' '.join(cmd[:3])} exit {r.returncode}; see {log}")


def run_local(tool: str, rows: list[dict], work: Path, suffix: str = "") -> dict:
    work.mkdir(parents=True, exist_ok=True)
    inp = work / f"{tool}_in.csv"
    write_in_csv(rows, inp)
    if tool == "liab":
        mm_run("tnp", "python", str(HERE / "liab.py"), "--in", str(inp), "--out", str(work / "liab.csv"), log=work / "liab.log")
        return parse(tool, work / "liab.csv", suffix)
    if tool == "nb9":
        mm_run("tnp", "python", str(HERE / "nanobody_rules.py"), "--in", str(inp), "--out", str(work / "nb9.csv"), log=work / "nb9.log")
        return parse(tool, work / "nb9.csv", suffix)
    if tool == "tnp":
        pdbs = work.parent / "tnp_pdbs"
        mm_run("tnp", *TNP_ENV, "python", str(HERE / "tnp_batch.py"), "--in", str(inp), "--out", str(work / "tnp.csv"),
               "--pdb-dir", str(pdbs), "--ncores", "8", log=work / "tnp.log")
        return parse(tool, work / "tnp.csv", suffix)
    if tool == "nanomelt":
        return run_nanomelt_local(rows, work, suffix)
    raise ValueError(tool)


def run_nanomelt_local(rows: list[dict], work: Path, suffix: str = "", chunk: int = 32) -> dict:
    """NanoMelt on the CPU (user decision: never on the 8 GB laptop GPU, where 248 sequences at once ran out of memory), in chunks of
    ~32 sequences, 8 threads (about 18 s for 4 sequences incl. model load).  A failing chunk is bisected down to the single
    offending sequence, which gets nm_error instead of poisoning the batch."""
    res: dict = {}
    counter = [0]

    def attempt(sub: list[dict]) -> None:
        counter[0] += 1
        fa, outc = work / f"nanomelt_{counter[0]}.fa", work / f"nanomelt_{counter[0]}.csv"
        fa.write_text("".join(f">{r['design_id']}\n{r['vhh_seq']}\n" for r in sub))
        try:
            mm_run("nanomelt", "env", "CUDA_VISIBLE_DEVICES=", "OMP_NUM_THREADS=8", "nanomelt", "predict", "-i", str(fa), "-align", "-o", str(outc),
                   log=work / f"nanomelt_{counter[0]}.log")
            res.update(parse("nanomelt", outc, suffix))
        except Exception as e:  # noqa: BLE001
            if len(sub) == 1:
                LOG(f"nanomelt: {sub[0]['design_id']} cannot be scored ({e})")
                res[sub[0]["design_id"]] = {"nm_tm" + suffix: "", "nm_error" + suffix: "nanomelt failed on this sequence"}
            else:
                h = len(sub) // 2
                attempt(sub[:h])
                attempt(sub[h:])

    for i in range(0, len(rows), chunk):
        attempt(rows[i:i + chunk])
    return res


def run_vm(tool: str, rows: list[dict], work: Path, suffix: str = "", dtype: str = "auto") -> dict:
    work.mkdir(parents=True, exist_ok=True)
    inp = work / f"{tool}_in.csv"
    write_in_csv(rows, inp)
    outd = work / f"vm_{tool}"
    if not ALLOW_VM[0]:
        raise RuntimeError("VM use is not allowed without --allow-vm (user permission)")
    r = subprocess.run(["bash", str(HERE / "run_vm.sh"), str(inp), str(outd), tool, dtype, "6"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                       env={**os.environ, "S4_ALLOW_VM": "1"})
    (work / f"vm_{tool}.log").write_text(r.stdout)
    if r.returncode:
        raise RuntimeError(f"VM {tool} failed (exit {r.returncode}); see {work / f'vm_{tool}.log'}")
    if tool == "monoyield":
        return parse(tool, outd, suffix)
    return parse(tool, outd / {"liab": "liab.csv", "nb9": "nb9.csv", "tnp": "tnp.csv", "nanomelt": "nanomelt.csv"}[tool], suffix)


def _f(x):
    try:
        return f"{float(x):.4g}" if x not in ("", None) else ""
    except ValueError:
        return x


def parse(tool: str, path: Path, suffix: str = "") -> dict:
    out: dict[str, dict] = {}
    sfx = (lambda c: c + suffix)
    if tool == "monoyield":
        y, m = path / "my_monoyield.csv", path / "my_mdratio.csv"
        st = (path / "status.txt").read_text() if (path / "status.txt").exists() else ""
        dt = "fp32" if "ok (fp32)" in st else ("bf16" if "ok (bf16)" in st else "")
        for f, col in ((y, "my_yield"), (m, "my_mdratio")):
            if f.exists():
                for r in dn.read_csv(f):
                    d = out.setdefault(r["name"], {})
                    d[sfx(col)] = _f(r["pred_mean"])
                    d[sfx(col + "_sd")] = _f(r["pred_std"])
                    d[sfx("my_dtype")] = dt
        return out
    rows = dn.read_csv(path)
    for r in rows:
        if tool == "nanomelt":
            out[r["ID"]] = {sfx("nm_tm"): _f(r["NanoMelt Tm (C)"])}
        else:
            did = r.pop("design_id")
            out[did] = {sfx(k): (_f(v) if k.startswith(("tnp_p", "tnp_cdr3_c")) else v) for k, v in r.items()}
    return out


# ------------------------------------------------------------------ main
def merge(state: dict, cands: dict, res: dict, tool: str) -> None:
    for did, cols in res.items():
        c = cands.get(did)
        if not c:
            continue
        s = state.setdefault(did, {})
        if s.get("vhh_seq") != c["vhh_seq"]:
            s.clear()  # sequence changed under the same id: drop every stale value
        s.update({"design_id": did, "vhh_seq": c["vhh_seq"], "generator": c["generator"], "framework": c["framework"],
                  "is_control": c["is_control"], "stage_ts": dn.utcnow()})
        s.update(cols)


def save(state: dict, path: Path) -> None:
    cols = ["design_id", "is_control", "generator", "framework", "vhh_seq", "stage_ts"]
    seen = set(cols)
    for s in state.values():
        for k in s:
            if k not in seen:
                cols.append(k)
                seen.add(k)
    for r in state.values():  # honest completion marker: a tool counts as done only if its value exists
        r["tools_done"] = ";".join(t for t in TOOLS if r.get(DONE_COL[t], "") != "")
    if "tools_done" not in seen:
        cols.append("tools_done")
    dn.write_csv(path, sorted(state.values(), key=lambda s: (int(s.get("is_control", 0)) == 0, s["design_id"])), cols)


def one_pass(a) -> int:
    out = Path(a.out)
    work = Path(a.work) if a.work else out.parent / "props_work"
    state = load_state(out)
    if a.redo:  # forget finished work of these tools (e.g. after a definition change); applies to the first pass only
        for t in a.redo:
            for r in state.values():
                r[DONE_COL[t]] = ""
        save(state, out)
        a.redo = None
    cands = candidates(a)
    nstale = reconcile(state, cands)
    if nstale:
        LOG(f"reconcile: dropped {nstale} rows whose design_id left the designs table (values carried to the id holding the same sequence)")
        save(state, out)
    ndone = 0
    where = a.where
    exp = expensive_ids(cands, a) if a.tier == "on" else set(cands)
    LOG(f"{len(cands)} designs; candidate set for tnp/monoyield: {len(exp)}")
    last_my = work / "monoyield_last.txt"
    for tool in a.tools:
        todo = pending(tool, cands, state, a.retry_failed)
        if tool in ("tnp", "monoyield"):
            todo = [r for r in todo if r["design_id"] in exp]
        if tool == "monoyield" and todo and a.watch and not a.final and last_my.exists():
            wait = a.my_interval - (time.time() - float(last_my.read_text() or 0))
            if wait > 0:
                LOG(f"monoyield: {len(todo)} pending, next VM call allowed in {wait:.0f}s")
                continue
        if tool == "tnp" and a.max_tnp is not None and len(todo) > a.max_tnp:
            key = lambda r: (not r["is_control"], -float(state.get(r["design_id"], {}).get("nm_tm") or -1))  # noqa: E731
            todo = sorted(todo, key=key)[: a.max_tnp]
        if not todo:
            continue
        mode = ("vm" if tool == "monoyield" or where == "vm" else "local")
        if tool == "monoyield" and where == "local":
            LOG("monoyield is VM-only (ESMC-6B does not fit the local 8 GB GPU); skipped under --where local")
            continue
        LOG(f"{tool}: {len(todo)} pending, running on {mode}")
        res, t0 = None, time.time()
        try:
            res = run_vm(tool, todo, work / tool, dtype=a.dtype) if mode == "vm" else run_local(tool, todo, work / tool)
        except Exception as e:  # noqa: BLE001
            LOG(f"{tool} on {mode} FAILED: {e}")
            if mode == "local" and where in ("auto", "both"):
                LOG(f"{tool}: falling back to the VM")
                try:
                    res = run_vm(tool, todo, work / tool, dtype=a.dtype)
                except Exception as e2:  # noqa: BLE001
                    LOG(f"{tool} VM fallback FAILED too: {e2}")
        if tool == "monoyield":
            work.mkdir(parents=True, exist_ok=True)
            last_my.write_text(str(time.time()))
        if res:
            merge(state, cands, res, tool)
            save(state, out)
            ndone += len(res)
            LOG(f"{tool}: {len(res)} rows in {time.time() - t0:.0f}s -> {out}")
        if where == "both" and tool in ("liab", "nanomelt", "tnp"):
            try:
                resv = run_vm(tool, todo, work / (tool + "_vmcheck"), suffix="_vm")
                merge(state, cands, resv, tool)
                save(state, out)
            except Exception as e:  # noqa: BLE001
                LOG(f"{tool} VM cross-check FAILED: {e}")
    return ndone


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(dn.OUT / "props.csv"))
    ap.add_argument("--work", default=None, help="scratch dir (default <out dir>/props_work)")
    ap.add_argument("--in-csv", action="append", help="extra sequences: design_id,vhh_seq (repeatable)")
    ap.add_argument("--no-designs", action="store_true", help="ignore results/denovo/designs_*.csv")
    ap.add_argument("--no-baselines", action="store_true", help="skip denovo/s4_props/baselines.csv")
    ap.add_argument("--ids", nargs="*", help="only these design_ids")
    ap.add_argument("--tools", nargs="+", default=TOOLS, choices=TOOLS)
    ap.add_argument("--where", choices=["auto", "local", "vm", "both"], default="auto")
    ap.add_argument("--dtype", choices=["auto", "fp32", "bf16"], default="bf16", help="monoyield precision on the VM")
    ap.add_argument("--tier", choices=["on", "off"], default="on", help="off = run every tool on every design")
    ap.add_argument("--score-col", default="of3_consensus_rank_score")
    ap.add_argument("--gate-col", default="denovo_vm_pass", help="boolean column of rules.csv that must be true for tnp/monoyield")
    ap.add_argument("--score-lower-better", dest="score_higher_better", action="store_false")
    ap.add_argument("--top-frac", type=float, default=0.35)
    ap.add_argument("--my-interval", type=int, default=600, help="min seconds between monoyield VM calls in --watch mode")
    ap.add_argument("--stop-at", default="15:45", help="local HH:MM at which --watch ends")
    ap.add_argument("--max-tnp", type=int, default=None)
    ap.add_argument("--retry-failed", action="store_true")
    ap.add_argument("--allow-vm", action="store_true",
                    help="USER PERMISSION REQUIRED: without it nothing runs on the shared VM (no monoyield, no VM NanoMelt/TNP)")
    ap.add_argument("--tnp-cpu", action="store_true", help="model with NanoBodyBuilder2 on the CPU (default: local GPU, peak 0.6 GiB)")
    ap.add_argument("--redo", nargs="+", choices=TOOLS, help="recompute these tools for every design (first pass)")
    ap.add_argument("--watch", type=int, default=0, help="loop forever, polling every N seconds")
    ap.add_argument("--final", action="store_true", help="with --watch: ignore the monoyield interval")
    a = ap.parse_args()
    if a.tnp_cpu:
        TNP_ENV[:] = ["env", "CUDA_VISIBLE_DEVICES="]
    if not a.allow_vm:  # user rule 14:33: nothing on the VM GPU without explicit permission
        if a.where != "local":
            LOG(f"--where {a.where} needs --allow-vm (user permission); forcing --where local")
            a.where = "local"
        if "monoyield" in a.tools:
            LOG("monoyield runs only on the VM GPU: paused (needs --allow-vm); keeping rows already in props.csv")
            a.tools = [t for t in a.tools if t != "monoyield"]
    ALLOW_VM[0] = a.allow_vm
    if a.watch:
        a.watch = min(a.watch, 300)
        hh, mm = map(int, a.stop_at.split(":"))
        stop_file = dn.DENOVO / "STOP"
        Path(HERE / "watch.pid").write_text(str(os.getpid()))
        LOG(f"watching every {a.watch}s until {a.stop_at} local or {stop_file}; pid {os.getpid()}")
        while True:
            now = datetime.now()
            if stop_file.exists() or (now.hour, now.minute) >= (hh, mm):
                LOG("stop condition reached, exiting")
                break
            try:
                n = one_pass(a)
                LOG(f"pass done ({n} rows)")
            except Exception as e:  # noqa: BLE001  never die inside the watch loop
                LOG(f"pass crashed: {type(e).__name__}: {e}")
            for _ in range(a.watch):
                if stop_file.exists():
                    break
                time.sleep(1)
        return
    one_pass(a)


if __name__ == "__main__":
    main()
