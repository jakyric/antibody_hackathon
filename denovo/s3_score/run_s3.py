"""S3 runner: refold VHH + HER2 domain IV complexes on JapanFold with OpenFold3 (primary) and orthogonal models,
parse every structure/PAE into interface metrics, and keep score_of3.csv / score_orth.csv up to date.

  denovo/.venv/bin/python denovo/s3_score/run_s3.py score  --models of3,dde,esm2   # everything unscored now
  denovo/.venv/bin/python denovo/s3_score/run_s3.py watch  --models of3,dde,esm2   # loop: poll designs_*.csv
  denovo/.venv/bin/python denovo/s3_score/run_s3.py calib  --models of3,dde,esm2   # calibration complexes only
  denovo/.venv/bin/python denovo/s3_score/run_s3.py assemble                         # rebuild the CSVs from raw metrics

Re-running is free: jobs are cached by content hash (dn_common.JapanFold.submit) and parsed metrics are append-only in
results/denovo/s3_raw/metrics.jsonl. The key is read from JAPANFOLD_API_KEY and never printed.
"""
from __future__ import annotations

import argparse
import os
from datetime import datetime
import concurrent.futures as cf
import hashlib
import json
import math
import sys
import threading
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "denovo"))
sys.path.insert(0, str(HERE))
import dn_common as dn  # noqa: E402
import geom as g  # noqa: E402

RAW = dn.OUT / "s3_raw"
RAW.mkdir(exist_ok=True)
METRICS = RAW / "metrics.jsonl"
CALIB = HERE / "calib_set.csv"
GATES = HERE / "gates.json"
LOCK = threading.Lock()

HER2 = (ROOT / "optimization" / "baseline" / "refold_inputs" / "her2_domIV.txt").read_text().strip()   # UniProt P04626 489-607
HEWL = ("KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVS"
        "DGNGMNAWVAWRNRCKGTDVQAWIRGCRL")                                                 # hen egg-white lysozyme, decoy antigen
assert len(HER2) == 119

_P = dict(write_pae=True, seed=0, output_format="cif", diffusion_samples=3)
MODELS = {   # key -> (JapanFold model id, column prefix, params)
    "of3":      ("openfold3",    "of3",   dict(_P, use_msa_server=True)),
    "dde":      ("opendde-abag", "dde",   dict(_P, use_msa_server=True)),
    "ddeg":     ("opendde",      "ddeg",  dict(_P, use_msa_server=True)),
    "ddenm":    ("opendde-abag", "ddenm", dict(_P, use_msa_server=False)),
    "esm2":     ("esmfold2",     "esm2",  dict(_P, use_msa_server=False)),
    "of3s1":    ("openfold3",    "of3s1", dict(_P, use_msa_server=True, seed=1)),   # second seed, same recipe: seed consistency
    "bz2":      ("boltz2",       "bz2",   dict(_P, use_msa_server=True)),
    "rf3":      ("rf3",          "rf3",   dict(_P, use_msa_server=True)),
    # AF2-IG re-scores a given complex: needs coordinates, takes no seed / samples / MSA
    "af2ig":    ("af2ig",        "af2ig", dict(write_pae=True, output_format="cif")),
}
BATCH = 40
BATCH_FOR = {"af2ig": 24}    # AF2-IG inputs carry ~150 kB of coordinates per target
PARSE_V = 4   # bump when parse_target changes: `reparse` redoes older entries from the stored structures


# ------------------------------------------------------------------ work items
IGNORED = ("fwdiv", "boltzgen")     # framework diversification cancelled; BoltzGen stopped by the user (prefix match)


def load_designs():
    """dn.load_designs() minus rows that are not design rows (designs_boltzgen_metrics.csv also matched the glob) and
    minus the generators that were cancelled."""
    return [d for d in dn.load_designs() if d.get("generator") and d.get("vhh_seq") and not d["generator"].startswith(IGNORED)]



def _yaml(chains):
    return "sequences:\n" + "".join(f"  - protein: {{id: {c}, sequence: {s}}}\n" for c, s in chains)


def item(tid, kind, vhh, vl, ant, ant_name, cdrs, design_cif="", group="design", label=None):
    chains = [("A", vhh)] + ([("B", vl)] if kind == "fv" else []) + [("C" if kind == "fv" else "B", ant)]
    return {"id": tid, "kind": kind, "vhh": vhh, "vl": vl, "ant": ant, "ant_name": ant_name, "cdrs": cdrs,
            "design_cif": design_cif, "yaml": _yaml(chains), "group": group, "label": label}


def design_items(decoy=True):
    out = []
    for d in load_designs():
        pj = json.loads(d.get("params_json") or "{}")
        cd = g.find_cdrs(d["vhh_seq"], d.get("cdr1", ""), d.get("cdr2", ""), d.get("cdr3", ""))
        if pj.get("kind") == "fv":
            out.append(item(f"{d['design_id']}__her2", "fv", d["vhh_seq"], pj["vl_seq"], HER2, "her2", {}, group="design"))
            continue
        if dn.check_seq(d["vhh_seq"]):
            continue
        out.append(item(f"{d['design_id']}__her2", "vhh", d["vhh_seq"], "", HER2, "her2", cd, d.get("design_cif", "")))
        if decoy:
            out.append(item(f"{d['design_id']}__hewl", "vhh", d["vhh_seq"], "", HEWL, "hewl", cd))
    return out


def af2ig_items(items):
    """AF2-IG inputs (target + binder coordinates in one structure, no MSA): the generator's own design when `design_cif`
    exists, else the OpenFold3 top model of that complex. Re-scores that pose; deterministic."""
    st = load_state()
    out = []
    for it in items:
        if it["kind"] != "vhh":
            continue
        src = ROOT / it["design_cif"] if it["design_cif"] and (ROOT / it["design_cif"]).exists() else None
        tag = "design"
        if src is None:
            r = st.get((it["id"], "of3"))
            if not r or r["status"] != "ok":
                continue
            src, tag = ROOT / r["metrics"]["cif"], "of3"
        try:
            her, vhh = g.split_complex(src, it["ant"])
            if len(vhh.ca().res) != len(it["vhh"]):
                continue
        except Exception:
            continue
        content = ("target:\n  structure: |\n" + "".join("    " + l + "\n" for l in g.to_pdb([("A", her), ("B", vhh)]).splitlines()) +
                   "  chain: A\nbinder:\n  sequence: " + it["vhh"] + "\n  chain: B\n")
        j = dict(it)
        j.update(yaml=content, af2ig=True, src_path=str(src), af2ig_src=tag)
        out.append(j)
    return out


def calib_items():
    out = []
    if not CALIB.exists():
        return out
    for r in dn.read_csv(CALIB):
        ant = r["antigen_seq"] or HER2
        cd = g.find_cdrs(r["vhh_seq"]) if r["kind"] == "vhh" else {}
        out.append(item(f"{r['calib_id']}__{r['antigen_name']}", r["kind"], r["vhh_seq"], r["vl_seq"], ant,
                        r["antigen_name"], cd, group=r["group"], label=int(r["label"])))
    return out


# ------------------------------------------------------------------ state
def load_state():
    st = {}
    if METRICS.exists():
        for line in METRICS.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                st[(r["key"], r["model"])] = r
    return st


def append(rows):
    with LOCK, open(METRICS, "a") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


# ------------------------------------------------------------------ parsing one target
def parse_af2ig(it, row, dest: Path):
    """AF2-IG output: chain A = target (HER2 residues present in the input, pooled), chain B = binder."""
    tid = it["id"]
    cifs = sorted(dest.glob(f"{tid}*.cif"), key=lambda p: (p.name != f"{tid}.cif", p.name))
    a = g.load_cif(cifs[0])
    her, vhh = a.of_chain("A"), a.of_chain("B")
    nh, nv = len(her.ca().res), len(vhh.ca().res)
    m = {"cif": str(cifs[0].relative_to(ROOT)), "src": it.get("af2ig_src"), "n_her_res": nh,
         "iptm": fnum(row.get("iptm")), "ptm": fnum(row.get("ptm")), "complex_plddt": fnum(row.get("plddt", row.get("complex_plddt"))),
         "interface_pae": fnum(row.get("interface_pae")), "confidence": fnum(row.get("confidence_score"))}
    ca = vhh.ca()
    m["binder_plddt"] = float(ca.b.mean() * (100.0 if ca.b.max() <= 1.0 else 1.0))
    pairs = g.contact_pairs(vhh, her)
    cdr_idx = [i for x, y in it["cdrs"].values() for i in range(x, y)]
    npz = dest / f"{tid}_pae.npz"
    if npz.exists():
        pae = np.load(npz)["pae"]
        idx = np.r_[nh:nh + nv, 0:nh]
        m.update(g.pae_metrics(pae[np.ix_(idx, idx)], nv, pairs, cdr_idx, [], n_her=nh))
    m["n_contact_pairs"] = len(pairs)
    # does AF2 keep the generator's pose? binder CA RMSD after superposing on the target (the RFD4 oracle criterion is <= 2 A)
    try:
        her_in, vhh_in = g.split_complex(it["src_path"], it["ant"])
        hi, vi = her_in.ca(), vhh_in.ca()
        if len(hi.xyz) == nh and len(vi.xyz) == nv:
            R, t = g.kabsch(her.ca().xyz, hi.xyz)
            m["pose_rmsd_vhh"] = g.rmsd(vhh.ca().xyz @ R.T + t, vi.xyz)
            if "cdr3" in it["cdrs"]:
                x, y = it["cdrs"]["cdr3"]
                m["pose_rmsd_cdr3"] = g.rmsd((vhh.ca().xyz @ R.T + t)[x:y], vi.xyz[x:y])
    except Exception as e:
        m["pose_error"] = str(e)[:80]
    return m


def parse_target(it, row, dest: Path):
    if it.get("af2ig"):
        return parse_af2ig(it, row, dest)
    tid = it["id"]
    cifs = sorted(dest.glob(f"{tid}*.cif"), key=lambda p: (p.name != f"{tid}.cif", p.name))
    if not cifs:
        raise FileNotFoundError(f"no structure for {tid} in {dest}")
    samples = [g.load_cif(p) for p in cifs]
    if it.get("swap"):   # AF2-IG writes the target (HER2) first, binder second: relabel to the S3 convention A=binder, B=HER2
        for s in samples:
            s.chain = np.where(s.chain == "A", "B", np.where(s.chain == "B", "A", s.chain))
    a0 = samples[0]
    fv = it["kind"] == "fv"
    her_c, extra = ("C", ("B",)) if fv else ("B", ())
    n_bind = len(it["vhh"]) + (len(it["vl"]) if fv else 0)
    is_her2 = it["ant_name"] == "her2"
    m = {"n_samples": len(samples), "cif": str(cifs[0].relative_to(ROOT))}
    runs = row.get("all_runs") or []
    if row.get("interface_pae") is not None:
        m["interface_pae"] = fnum(row.get("interface_pae"))
    m.update(iptm=fnum(row.get("iptm")), ptm=fnum(row.get("ptm")), complex_plddt=fnum(row.get("complex_plddt")),
             confidence=fnum(row.get("confidence_score")),
             iptm_sd=float(np.std([fnum(r.get("iptm")) for r in runs])) if len(runs) > 1 else float("nan"),
             msa_depth=row.get("msa_depth"))
    ca_len = {c: len(a0.of_chain(c).ca().res) for c in a0.chains()}
    if ca_len.get("A") != len(it["vhh"]) or ca_len.get(her_c) != len(it["ant"]):
        m["warn"] = f"chain lengths {ca_len}"
    vhh_atoms = a0.of_chain("A") if not fv else a0.sel(np.isin(a0.chain, ["A", "B"]))
    pairs = g.contact_pairs(a0.of_chain("A"), a0.of_chain(her_c))
    if fv:   # pooled VH + VL contacts, PAE rows of VL come after the VH rows
        pairs = pairs | {(i + len(it["vhh"]), j) for i, j in g.contact_pairs(a0.of_chain("B"), a0.of_chain(her_c))}
    cdr_idx = [i for a, b in it["cdrs"].values() for i in range(a, b)]
    epi = g.reference()["epitope"] if is_her2 else []
    if is_her2:
        m.update(g.complex_metrics(a0, it["cdrs"], vhh_chain="A", her_chain=her_c, extra_ab_chains=extra))
    else:
        m["n_contact_pairs"] = len(pairs)
        ca = a0.of_chain("A").ca()
        m["binder_plddt"] = float(ca.b.mean() * (100.0 if ca.b.max() <= 1.0 else 1.0))
    npz = dest / f"{tid}_pae.npz"
    if npz.exists():
        pae = np.load(npz)["pae"]
        if it.get("swap"):
            nh = len(it["ant"])
            idx = np.r_[nh:pae.shape[0], 0:nh]
            pae = pae[np.ix_(idx, idx)]
        m.update(g.pae_metrics(pae, n_bind, pairs, [] if fv else cdr_idx, epi, n_her=len(it["ant"])))
    # across samples: pose stability (VHH, HER2 frame)
    if is_her2 and not fv and len(samples) > 1:
        rm, cd3, ep = [], [], [m.get("epi_overlap", float("nan"))]
        for s in samples[1:]:
            pv = g.pose_vs(a0, s, it["cdrs"].get("cdr3"))
            rm.append(pv["vhh_rmsd_her2frame"])
            cd3.append(pv.get("cdr3_centroid_dist", float("nan")))
            ep.append(g.complex_metrics(s, it["cdrs"])["epi_overlap"])
        m.update(pose_rmsd_max=float(np.max(rm)), cdr3_cent_max=float(np.nanmax(cd3)), epi_overlap_min=float(np.min(ep)))
    if it["design_cif"] and is_her2 and not fv:
        p = ROOT / it["design_cif"]
        if p.exists():
            try:
                m.update(g.own_design_rmsd(p, a0, HER2))
            except Exception as e:  # keep the row, note the failure
                m["own_error"] = str(e)[:80]
    return m


# ------------------------------------------------------------------ jobs
_jf = None


def jf():
    global _jf
    if _jf is None:
        _jf = dn.JapanFold()
    return _jf


def run_job(mk, items, jname=None):
    model, _, params = MODELS[mk]
    names = sorted(i["id"] for i in items)
    jname = jname or f"s3-{mk}-" + hashlib.sha1((json.dumps(params, sort_keys=True) + ",".join(names)).encode()).hexdigest()[:10]
    body = {"model": model, "name": jname, "targets": [{"name": i["id"], "content": i["yaml"]} for i in items],
            "params": params}
    try:
        jid = jf().submit(jname, "/v1/predictions", body)
    except RuntimeError as e:      # e.g. HTTP 429 "Max queue length": back off, the next cycle plans it again
        print(f"[{mk}] submit refused ({str(e)[:120]}), will retry later", flush=True)
        return
    print(f"[{mk}] {jname}: job {jid}, {len(items)} targets", flush=True)
    try:
        j = jf().wait(jid, max_s=3300, log=lambda s: print(f"[{mk}] {s}", flush=True))
    except TimeoutError as e:
        print(f"[{mk}] {e}", flush=True)
        return
    if j["status"] != "succeeded":
        err = json.dumps(j)[:300]
        print(f"[{mk}] job {jid} ended {j['status']}: {err}", flush=True)
        append([{"key": i["id"], "model": mk, "status": "job_failed", "job_id": jid, "error": err, "ts": dn.utcnow()} for i in items])
        return
    res = jf().results(jid)
    dest = RAW / mk / jid
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "s3_results.json").write_text(json.dumps(res))
    if not any(dest.glob("*.cif")) if dest.exists() else True:
        jf().download(jid, dest)
    rows_by_id = {r.get("id"): r for r in res.get("rows", [])}
    out = []
    for it in items:
        row = rows_by_id.get(it["id"])
        if not row or row.get("status") != "ok":
            out.append({"key": it["id"], "model": mk, "status": "target_failed", "job_id": jid,
                        "error": json.dumps(row)[:300] if row else "missing row", "ts": dn.utcnow()})
            continue
        try:
            mt = parse_target(it, row, dest)
            if it.get("af2ig_src"):
                mt["src"] = it["af2ig_src"]
            out.append({"key": it["id"], "model": mk, "status": "ok", "job_id": jid, "v": PARSE_V, "metrics": mt,
                        "runtime_s": row.get("runtime_s"), "ts": dn.utcnow()})
        except Exception as e:
            out.append({"key": it["id"], "model": mk, "status": "parse_failed", "job_id": jid, "error": repr(e)[:300], "ts": dn.utcnow()})
    append(out)
    n_ok = sum(r["status"] == "ok" for r in out)
    print(f"[{mk}] job {jid}: {n_ok}/{len(items)} parsed", flush=True)
    if res.get("partial"):
        print(f"[{mk}] PARTIAL {json.dumps(res['partial'])}", flush=True)


def run_job_safe(mk, items, jname=None):
    """run_job with retries on transient errors (HTTP 429 'max queue length', DNS blips); the job name is cached, so a retry resumes the same job."""
    for attempt in range(5):
        try:
            return run_job(mk, items, jname)
        except Exception as e:
            print(f"[{mk}] transient error ({repr(e)[:140]}), retry {attempt + 1}/5 in 30 s", flush=True)
            time.sleep(30)


def score(models, items, workers=6, max_fail=2):
    st = load_state()
    todo = {}
    for mk in models:
        fresh = []
        for it in (af2ig_items(items) if mk == "af2ig" else items):
            r = st.get((it["id"], mk))
            fails = sum(1 for (k, m), x in st.items() if k == it["id"] and m == mk and x["status"] != "ok")
            if r and r["status"] == "ok":
                continue
            if r and fails >= max_fail:
                continue
            fresh.append(it)
        todo[mk] = fresh
    jobs = [(mk, v[i:i + BATCH_FOR.get(mk, BATCH)]) for mk, v in todo.items() for i in range(0, len(v), BATCH_FOR.get(mk, BATCH))]
    if not jobs:
        print("nothing to score")
        return 0
    c0 = jf().credits()
    print(f"{sum(len(v) for v in todo.values())} model-targets in {len(jobs)} jobs; credits ${c0:.3f}", flush=True)
    with cf.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_job, mk, its) for mk, its in jobs]
        for f in futs:
            try:
                f.result()
            except Exception as e:
                print("job error:", repr(e)[:300], flush=True)
    c1 = jf().credits()
    with open(HERE / "spend.md", "a") as fh:
        fh.write(f"- {dn.utcnow()} models={','.join(models)} targets={sum(len(v) for v in todo.values())} jobs={len(jobs)} "
                 f"credits ${c0:.3f} -> ${c1:.3f} (spent ${c0 - c1:.3f})\n")
    print(f"credits ${c1:.3f} (spent ${c0 - c1:.3f})")
    return len(jobs)


# ------------------------------------------------------------------ scoring tables
STOP_FILE = ROOT / "denovo" / "STOP"
STOP_HHMM = tuple(int(x) for x in __import__("os").environ.get("PIPELINE_S3_UNTIL", "15:40").split(":"))                      # stop launching new jobs at 15:40 local (CET)
BUDGET = float(__import__("os").environ.get("PIPELINE_BUDGET_S3", "25"))                            # USD of JapanFold credit this stage may spend
DEFAULT_METRIC = {"of3": ("ipae_tm", "high")}   # used for percentiles until calibrate.py writes gates.json
OF3_CONSENSUS = (("ipsae_chn15", "high"), ("ipae_tm", "high"), ("min_ipae", "low"), ("iptm", "high"))   # a-priori triple; see decisions.md
PREF = {mk: MODELS[mk][1] for mk in MODELS}


def gates():
    return json.loads(GATES.read_text()) if GATES.exists() else {}


def rules_pass():
    """design_id -> bool from S2's results/denovo/rules.csv (denovo_vm_pass AND canon_cys_ok); missing id = no verdict yet."""
    f = dn.OUT / "rules.csv"
    out = {}
    if f.exists():
        try:
            for x in dn.read_csv(f):
                t = lambda k: str(x.get(k, "")).strip().lower() in ("1", "true")
                out[x["design_id"]] = t("denovo_vm_pass") and t("canon_cys_ok")   # de novo gate AND canonical disulfide present
        except Exception:
            pass
    return out


def pct_of(v, pop, direction="high"):
    """Share of the population strictly worse than v (ties count half): 1 = best, in [0, 1]."""
    pop = [x for x in pop if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not pop or v is None or (isinstance(v, float) and math.isnan(v)):
        return float("nan")
    if direction == "low":
        worse, tie = sum(x > v for x in pop), sum(x == v for x in pop)
    else:
        worse, tie = sum(x < v for x in pop), sum(x == v for x in pop)
    return (worse + 0.5 * tie) / len(pop)


# ------------------------------------------------------------------ de-duplication by vhh_seq
SEQ_CACHE = RAW / "id_seq.json"


def seq_cache(designs):
    """design_id -> vhh_seq for every id ever scored (ids that left the table are recovered from their stored structure)."""
    cache = json.loads(SEQ_CACHE.read_text()) if SEQ_CACHE.exists() else {}
    cache.update({d: row["vhh_seq"] for d, row in designs.items()})
    return cache


def fill_cache(st, designs, cache):
    changed = False
    for (key, mk), r in st.items():
        did = key.rsplit("__", 1)[0]
        if did in cache or r["status"] != "ok" or mk == "af2ig" or "cif" not in r.get("metrics", {}) or key.startswith(("cal-", "ctl-tz")):
            continue
        try:
            cache[did] = g.load_cif(ROOT / r["metrics"]["cif"]).seq("A")
            changed = True
        except Exception:
            pass
    if changed or not SEQ_CACHE.exists():
        SEQ_CACHE.write_text(json.dumps(cache))
    return cache


def canon_ids(designs):
    """One id per vhh_seq among the table's designs (earliest created_utc, then id): the only one that is ever scored or ranked."""
    by = {}
    for d in sorted(designs, key=lambda d: (designs[d].get("created_utc", ""), d)):
        by.setdefault(designs[d]["vhh_seq"], []).append(d)
    return {v[0]: v[1:] for v in by.values()}, by


def canonical_view(st, designs):
    """State re-keyed onto the current table: a score belongs to the sequence, so it is copied to every id in the table that carries
    the same vhh_seq, and scores of ids that left the table still count. Returns (view, stats)."""
    cache = fill_cache(st, designs, seq_cache(designs))
    _, by = canon_ids(designs)
    view = {}
    for (key, mk), r in st.items():          # direct hits first
        did = key.rsplit("__", 1)[0]
        if did in designs:
            view[(key, mk)] = r
    for (key, mk), r in st.items():
        did, ant = key.rsplit("__", 1) if "__" in key else (key, "")
        seq = cache.get(did)
        if seq is None or r["status"] != "ok":
            continue
        for t in by.get(seq, []):
            view.setdefault((f"{t}__{ant}", mk), r)
    scored = [k.rsplit("__", 1)[0] for (k, mk), r in st.items() if mk == "of3" and r["status"] == "ok" and k.endswith("__her2")
              and not k.startswith(("cal-", "ctl-"))]
    seqs = [cache.get(i) for i in scored if cache.get(i)]
    stats = {"of3_rows": len(scored), "unique_seqs": len(set(seqs)), "dup_rows": len(seqs) - len(set(seqs)),
             "ids_left_table": sum(1 for i in set(scored) if i not in designs)}
    return view, stats


def build(st=None):
    """Per-design wide dict from the raw metric store. Returns (per, designs)."""
    st = st or load_state()
    gt = gates().get("models", {})
    designs = {d["design_id"]: d for d in load_designs()}
    st, _stats = canonical_view(st, designs)      # duplicates share one score; ids that left the table drop out
    dups, _by = canon_ids(designs)
    dup_of = {x: c for c, xs in dups.items() for x in xs}
    rp = rules_pass()
    per = {}
    for (key, mk), r in st.items():
        if r["status"] != "ok" or "__" not in key:
            continue
        did, ant = key.rsplit("__", 1)
        if did not in designs:
            continue
        p, m = PREF[mk], r["metrics"]
        d = per.setdefault(did, {})
        if ant == "her2":
            for k, v in m.items():
                if k in ("warn", "own_error", "msa_depth", "n_samples") or isinstance(v, str) and k not in ("cif", "src"):
                    continue
                d[f"{p}_{k}"] = v
        else:
            for k in ("iptm", "ipsae15", "min_ipae", "n_contact_pairs"):
                if k in m:
                    d[f"{p}_{ant}_{k}"] = m[k]
    for did, d in per.items():
        if "of3_iptm" in d and "of3_hewl_iptm" in d:
            d["of3_spec_iptm_delta"] = d["of3_iptm"] - d["of3_hewl_iptm"]
    # cross-model pose agreement on the top sample of each model, anchored on OF3
    for did, d in per.items():
        if '"fv"' in designs[did].get("params_json", ""):
            continue
        cifs = {mk: st[(f"{did}__her2", mk)]["metrics"]["cif"] for mk in MODELS
                if (f"{did}__her2", mk) in st and st[(f"{did}__her2", mk)]["status"] == "ok"}
        if len(cifs) > 1:
            cd = g.find_cdrs(designs[did]["vhh_seq"], designs[did].get("cdr1", ""), designs[did].get("cdr2", ""), designs[did].get("cdr3", ""))
            ref_mk = "of3" if "of3" in cifs else sorted(cifs)[0]
            try:
                ra = g.load_cif(ROOT / cifs[ref_mk])
            except Exception:
                continue
            rm, c3 = [], []
            for mk, c in cifs.items():
                if mk == ref_mk or mk in ("af2ig", "of3s1"):
                    continue
                try:
                    pv = g.pose_vs(ra, g.load_cif(ROOT / c), cd.get("cdr3"))
                except Exception:
                    continue
                d[f"x_{PREF[mk]}_vs_{PREF[ref_mk]}_vhh_rmsd"] = pv["vhh_rmsd_her2frame"]
                d[f"x_{PREF[mk]}_vs_{PREF[ref_mk]}_cdr3_dist"] = pv.get("cdr3_centroid_dist", float("nan"))
                rm.append(pv["vhh_rmsd_her2frame"])
                c3.append(pv.get("cdr3_centroid_dist", float("nan")))
            if rm:
                d["xmodel_vhh_rmsd_mean"] = float(np.mean(rm))
                d["xmodel_cdr3_dist_max"] = float(np.nanmax(c3))
                d["xmodel_n_agree"] = int(sum(x < 8.0 for x in c3))   # other models whose CDR3 sits within 8 A of OF3's
    # seed consistency: the same OF3 recipe with a second seed
    for did, d in per.items():
        c0, c1 = (st.get((f"{did}__her2", mk)) for mk in ("of3", "of3s1"))
        if c0 and c1 and c0["status"] == "ok" and c1["status"] == "ok" and '"fv"' not in designs[did].get("params_json", ""):
            cd = g.find_cdrs(designs[did]["vhh_seq"], designs[did].get("cdr1", ""), designs[did].get("cdr2", ""), designs[did].get("cdr3", ""))
            try:
                pv = g.pose_vs(g.load_cif(ROOT / c0["metrics"]["cif"]), g.load_cif(ROOT / c1["metrics"]["cif"]), cd.get("cdr3"))
                d["of3_seed_vhh_rmsd"] = pv["vhh_rmsd_her2frame"]
                d["of3_seed_cdr3_dist"] = pv.get("cdr3_centroid_dist", float("nan"))
                d["of3_seed_iptm_min"] = min(c0["metrics"]["iptm"], c1["metrics"]["iptm"])
                d["of3_seed_consistent"] = int(d["of3_seed_cdr3_dist"] < 6.0 and min(c0["metrics"].get("epi_overlap", 0), c1["metrics"].get("epi_overlap", 0)) >= 0.5)
            except Exception:
                pass
    # percentile populations = non-control designs that passed the VM abguard de novo gate (S2 rules.csv)
    # percentile populations = unique, rule-passing, non-control designs (duplicates would double-count)
    cand = [d for d in per if designs[d]["generator"] != "control" and rp.get(d) and d not in dup_of]
    # OF3 consensus: mean percentile of three confidence metrics against the non-control designs scored so far
    pops = {met: [per[c].get(f"of3_{met}") for c in cand] for met, _ in OF3_CONSENSUS}
    for did, d in per.items():
        ps = [pct_of(d.get(f"of3_{met}"), pops[met], dirn) for met, dirn in OF3_CONSENSUS if f"of3_{met}" in d]
        ps = [x for x in ps if not math.isnan(x)]
        if ps:
            d["of3_consensus_rank_score"] = float(np.mean(ps))
    # per-model primary-metric percentile and gate pass flags
    models_seen = {mk for (k, mk), r in st.items() if r["status"] == "ok" and mk not in ("of3s1",)}
    for mk in models_seen:
        p = PREF[mk]
        gate = gt.get(mk)
        met, dirn = (gate["metric"], gate["direction"]) if gate else DEFAULT_METRIC.get(mk, ("ipae_tm", "high"))
        pop = [per[c].get(f"{p}_{met}") for c in cand]
        for did, d in per.items():
            v = d.get(f"{p}_{met}")
            if v is None:
                continue
            d[f"{p}_pct"] = pct_of(v, pop, dirn)
            if gate and not (isinstance(v, float) and math.isnan(v)):
                ok = (v >= gate["thr"]) if dirn == "high" else (v <= gate["thr"])
                if gate.get("epi_min") is not None and f"{p}_epi_overlap" in d:
                    ok = ok and d[f"{p}_epi_overlap"] >= gate["epi_min"]
                d[f"{p}_pass"] = int(ok)
    for did, d in per.items():
        pcts = [d[f"{PREF[mk]}_pct"] for mk in models_seen if f"{PREF[mk]}_pct" in d and not math.isnan(d[f"{PREF[mk]}_pct"])]
        d["consensus_rank_score"] = float(np.mean(pcts)) if pcts else float("nan")
        d["n_models_scored"] = len(pcts)
        passes = [d[f"{PREF[mk]}_pass"] for mk in models_seen if f"{PREF[mk]}_pass" in d]
        d["n_models_pass"] = int(sum(passes)) if passes else float("nan")
        if any(k.startswith("of3_") for k in d):
            d["of3_recipe"] = "openfold3 JapanFold; MSA server on; no template; 3 samples; seed 0; top-ranked sample"
        d["tier"] = 2 if any(k.startswith(("dde_", "esm2_", "af2ig_", "of3s1_")) for k in d) else 1
        d["denovo_vm_pass"] = int(rp[did]) if did in rp else float("nan")
        d["dup_of"] = dup_of.get(did, "")
    bgf = dn.OUT / "designs_boltzgen_metrics.csv"   # BoltzGen's own (Boltz-2 based) metrics: informational, never a gate
    if bgf.exists():
        for br in dn.read_csv(bgf):
            if br.get("design_id") in per:
                per[br["design_id"]].update({f"bg_{k}": v for k, v in br.items() if k != "design_id"})
    return per, designs


COMMON = ["design_id", "dup_of", "generator", "framework", "hotspot_set", "parent_id", "is_control", "denovo_vm_pass", "tier", "consensus_rank_score",
          "n_models_pass", "n_models_scored", "of3_consensus_rank_score", "stage_ts"]


def assemble():
    per, designs = build()
    ts = dn.utcnow()
    rows_of3, rows_orth = [], []
    for did, d in sorted(per.items()):
        ds = designs[did]
        base = {"design_id": did, "dup_of": d.get("dup_of", ""), "generator": ds["generator"], "framework": ds.get("framework", ""),
                "hotspot_set": ds.get("hotspot_set", ""), "parent_id": ds.get("parent_id", ""),
                "is_control": int(ds["generator"] == "control"), "stage_ts": ts,
                **{k: d.get(k, float("nan")) for k in ("denovo_vm_pass", "tier", "consensus_rank_score", "n_models_pass", "n_models_scored",
                                                       "of3_consensus_rank_score")}}
        rows_of3.append({**base, **{k: v for k, v in d.items() if k.startswith("of3_")}})
        rows_orth.append({**base, **{k: v for k, v in d.items() if not k.startswith("of3_")},
                          **{k: v for k, v in d.items() if k in ("of3_pass", "of3_pct")}})
    for path, rows in ((dn.OUT / "score_of3.csv", rows_of3), (dn.OUT / "score_orth.csv", rows_orth)):
        if rows:
            cols = COMMON + sorted({k for r in rows for k in r} - set(COMMON))
            dn.write_csv(path, rows, cols)
    return len(rows_of3)


# ------------------------------------------------------------------ reparse (new metrics on stored structures, no API cost)
def reparse():
    st = load_state()
    base = {i["id"]: i for i in design_items() + calib_items()}
    af = {i["id"]: i for i in af2ig_items(list(base.values()))}
    out, n = [], 0
    by_job = {}
    for (key, mk), r in st.items():
        if r["status"] == "ok" and r.get("v", 1) < PARSE_V and key in (af if mk == "af2ig" else base):
            by_job.setdefault((mk, r["job_id"]), []).append(key)
    for (mk, jid), keys in by_job.items():
        items = af if mk == "af2ig" else base
        dest = RAW / mk / jid
        rf = dest / "s3_results.json"
        res = json.loads(rf.read_text()) if rf.exists() else jf().results(jid)
        if not rf.exists():
            rf.write_text(json.dumps(res))
        rows = {r.get("id"): r for r in res.get("rows", [])}
        for k in keys:
            try:
                mt = parse_target(items[k], rows[k], dest)
                if items[k].get("af2ig_src"):
                    mt["src"] = items[k]["af2ig_src"]
                out.append({"key": k, "model": mk, "status": "ok", "job_id": jid, "v": PARSE_V, "metrics": mt, "ts": dn.utcnow()})
                n += 1
            except Exception as e:
                print("reparse failed", k, mk, repr(e)[:120])
    append(out)
    print(f"reparsed {n} targets")


# ------------------------------------------------------------------ continuous watch
def now_hhmm():
    t = datetime.now()
    return (t.hour, t.minute)


def log(msg, fh=None):
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)


def tier_plan(st, per, designs, tier2_models):
    """(model -> [item]) still to launch: tier 1 = OF3 on every design (controls first, then newest);
    tier 2 = other models + the lysozyme decoy for controls, top ~35% of OF3 consensus, and all framework-diversified variants."""
    its = {i["id"]: i for i in design_items(decoy=True)}
    st, _ = canonical_view(st, designs)      # a sequence that was scored under any id is never scored again
    dups, _by = canon_ids(designs)
    dup_ids = {x for xs in dups.values() for x in xs}
    rp = rules_pass()           # tier 1 waits for S2's verdict and only takes designs that pass the VM abguard de novo gate
    # representative first: round-robin over (framework, hotspot set) groups, newest first inside each group
    groups = {}
    for d in sorted((d for d in designs if designs[d]["generator"] != "control" and rp.get(d) and d not in dup_ids),
                    key=lambda d: designs[d].get("created_utc", ""), reverse=True):
        groups.setdefault((designs[d].get("framework", ""), designs[d].get("hotspot_set", "")), []).append(d)
    newest = []
    while any(groups.values()):
        for k in sorted(groups):
            if groups[k]:
                newest.append(groups[k].pop(0))
    ctl = sorted(d for d in designs if designs[d]["generator"] == "control")
    seq = ctl + newest
    cand = [d for d in per if designs[d]["generator"] != "control" and rp.get(d) and d not in dup_ids and "of3_consensus_rank_score" in per[d]]
    k = max(3, math.ceil(0.35 * len(cand))) if cand else 0
    top = set(sorted(cand, key=lambda d: -per[d]["of3_consensus_rank_score"])[:k])
    t2 = set(ctl) | top

    def need(key, mk):
        r = st.get((key, mk))
        if r and r["status"] == "ok":
            return False
        return sum(1 for (kk, m), x in st.items() if kk == key and m == mk and x["status"] != "ok") < 2
    plan = {"of3": []}
    for d in seq:
        if f"{d}__her2" in its and need(f"{d}__her2", "of3"):
            plan["of3"].append(its[f"{d}__her2"])
    for mk in tier2_models:
        plan[mk] = []
    for d in seq:
        if d not in t2:
            continue
        if f"{d}__hewl" in its and need(f"{d}__hewl", "of3") and (f"{d}__her2", "of3") in st:
            plan["of3"].append(its[f"{d}__hewl"])
        for mk in tier2_models:
            if mk != "af2ig" and f"{d}__her2" in its and need(f"{d}__her2", mk):
                plan[mk].append(its[f"{d}__her2"])
    if "af2ig" in plan:      # tier-1 companion: every gate-passer with a design pose (cheap, deterministic) + controls with an OF3 pose
        base = [its[f"{d}__her2"] for d in seq if f"{d}__her2" in its and need(f"{d}__her2", "af2ig")]
        plan["af2ig"] = af2ig_items(base)
    return plan, len(top), len(t2)


def funnel_lines(per, designs, st, spent, cyc, plan_sz, n_top):
    rp = rules_pass()
    dups, _by = canon_ids(designs)
    dup_ids = {x for xs in dups.values() for x in xs}
    uniq = {d: v for d, v in designs.items() if d not in dup_ids}          # unique vhh_seq only
    cand = [d for d in per if designs[d]["generator"] != "control" and rp.get(d) and d not in dup_ids]
    gens = {}
    for d in uniq.values():
        if d["generator"] != "control":
            gens[d["generator"]] = gens.get(d["generator"], 0) + 1
    n_rule = sum(1 for d in uniq.values() if d["generator"] != "control" and rp.get(d["design_id"]))
    n_wait = sum(1 for d in uniq.values() if d["generator"] != "control" and d["design_id"] not in rp)
    n_ctl = sum(1 for d in designs.values() if d["generator"] == "control")
    sc = lambda mk: sum(1 for d in cand if f"{PREF[mk]}_iptm" in per[d])
    L = [f"cycle {cyc} | unique designs {sum(gens.values())} (+{len(dup_ids)} duplicate ids not scored) ({', '.join(f'{k} {v}' for k, v in sorted(gens.items())) or 'none yet'}) + {n_ctl} controls | "
         f"abguard-VM pass {n_rule} (awaiting verdict {n_wait}) | OF3 scored {sc('of3')} | tier-2 (OF3 top {n_top}) dde {sc('dde')} esm2 {sc('esm2')} | to launch {plan_sz} | spend ${spent:.2f}/{BUDGET:.0f}"]
    epi = [d for d in cand if per[d].get("of3_epi_overlap", 0) >= 0.5]
    p3 = [d for d in cand if per[d].get("of3_pass") == 1]
    p2 = [d for d in cand if per[d].get("n_models_pass", 0) >= 2]
    comp = [d for d in epi if per[d].get("of3_clash30", 0) >= 5 or per[d].get("of3_cdr_fv_frac4", 0) >= 0.15]
    L.append(f"funnel: OF3-scored {sc('of3')} -> epitope overlap >= 0.5: {len(epi)} -> overlapping trastuzumab Fv: {len(comp)} -> past OF3 gate: {len(p3)} -> 2+ models pass: {len(p2)}")
    gt = gates().get("models", {})
    ctl_ids = [d for d in per if designs[d]["generator"] == "control"]
    for mk in ("of3", "dde", "esm2"):
        p = PREF[mk]
        gate = gt.get(mk)
        met = gate["metric"] if gate else DEFAULT_METRIC.get(mk, ("ipae_tm", "high"))[0]
        col = f"{p}_{met}"
        tz = per.get("ctl-tz-fv", {}).get(col)
        dec = [per[d][col] for d in ctl_ids if d.startswith("ctl-fwrand") and col in per[d]]
        des = [per[d][col] for d in cand if col in per[d]]
        fmt = lambda x: "-" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.2f}"
        thr = fmt(gate["thr"]) if gate else "uncalibrated"
        L.append(f"gate {mk} [{met} {'>=' if not gate or gate['direction'] == 'high' else '<='} {thr}]: tz-fv {fmt(tz)} | decoys max {fmt(max(dec) if dec else None)} | "
                 f"designs n={len(des)} median {fmt(float(np.median(des)) if des else None)} past gate {sum(1 for d in cand if per[d].get(f'{p}_pass') == 1)}")
    # canary: do the gates separate the controls from the designs?
    ctl_n = sum(1 for d in designs.values() if d["generator"] == "control")
    got = {mk: sum(1 for d in ctl_ids if f"{PREF[mk]}_iptm" in per[d]) for mk in ("of3", "dde", "esm2")}
    fz = lambda x: "-" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.2f}"
    tz = per.get("ctl-tz-fv", {})
    fw = [per[d].get("of3_consensus_rank_score") for d in ctl_ids if d.startswith("ctl-fwrand")]
    fw = [x for x in fw if x is not None and not math.isnan(x)]
    hw = [1 for d in cand if "of3_hewl_iptm" in per[d]]
    L.append(f"canary: controls scored of3 {got['of3']}/{ctl_n} dde {got['dde']}/{ctl_n} esm2 {got['esm2']}/{ctl_n} | "
             f"tz-fv of3 iptm {fz(tz.get('of3_iptm'))} consensus pct {fz(tz.get('of3_consensus_rank_score'))} | "
             f"fwrand decoys consensus pct median {fz(float(np.median(fw)) if fw else None)} | designs with lysozyme decoy {len(hw)}")
    return L


_DONE_JOBS = set()


def adopt(ex, futs, inflight):
    """Resume jobs that were submitted earlier (by this stage, before a restart) and are not parsed yet."""
    st = load_state()
    items = {i["id"]: i for i in design_items(decoy=True) + calib_items()}
    n = 0
    for f in sorted(dn.JOBS.glob("s3-*.json")):
        jid, mk = json.loads(f.read_text())["job_id"], f.stem.split("-")[1]
        if mk not in MODELS or mk == "af2ig" or jid in _DONE_JOBS:
            continue
        try:
            j = json.loads(jf().call("GET", f"/v1/jobs/{jid}"))
        except Exception:
            continue
        ids = [t["id"] for t in j.get("targets", [])]
        if not ids or all((i, mk) in st and st[(i, mk)]["status"] == "ok" for i in ids):
            _DONE_JOBS.add(jid)
            continue
        its = [items[i] for i in ids if i in items]
        if not its or j.get("model") != MODELS[mk][0]:
            continue
        key = (mk, tuple(i["id"] for i in its))
        if key in futs:
            continue
        futs[key] = ex.submit(run_job_safe, mk, its, f.stem)
        inflight |= {(i["id"], mk) for i in its}
        n += 1
    return n


def watch(a):
    tier2 = [m for m in a.tier2.split(",") if m]
    jf0 = jf().credits()
    ex = cf.ThreadPoolExecutor(max_workers=a.workers)
    futs = {}      # (model, tuple(ids)) -> future
    inflight = set()
    cyc = 0
    log(f"watch start pid {os.getpid()}: tier1 of3, tier2 {tier2}, budget ${BUDGET}, credits ${jf0:.3f}")
    log(f"adopted {adopt(ex, futs, inflight)} earlier jobs that were still unparsed")
    _designs = {d["design_id"]: d for d in load_designs()}
    _, _stats = canonical_view(load_state(), _designs)
    log(f"dedupe by vhh_seq: of {_stats['of3_rows']} already-scored OF3 design rows, {_stats['dup_rows']} were duplicates of another scored id "
        f"({_stats['unique_seqs']} unique sequences); {_stats['ids_left_table']} scored ids are no longer in the designs table (their scores are "
        f"kept per sequence, not per id); the table now has {len(_designs)} ids / {len({d['vhh_seq'] for d in _designs.values()})} unique sequences")
    last_spent = 0.0
    while True:
        cyc += 1
        try:
            for k, f in list(futs.items()):
                if f.done():
                    try:
                        f.result()
                    except Exception as e:
                        log(f"job error {k[0]}: {repr(e)[:200]}")
                    inflight -= {(i, k[0]) for i in k[1]}
                    del futs[k]
            stop = STOP_FILE.exists() or now_hhmm() >= STOP_HHMM
            if cyc % 6 == 0:      # pick up jobs whose thread died (e.g. repeated network errors)
                adopt(ex, futs, inflight)
            try:
                last_spent = jf0 - jf().credits()
            except Exception as e:
                log(f"credits unavailable ({str(e)[:80]}), using last value")
            spent = last_spent
            st = load_state()
            per, designs = build(st)
            plan, n_top, n_t2 = tier_plan(st, per, designs, tier2)
            if not stop and spent < BUDGET - 0.5:
                # priority: tier 1 OF3, then tier-2 models, AF2-IG last; never more than --max-jobs jobs in flight so the
                # shared JapanFold account queue (HTTP 429 when full) stays free for the other stages
                caps = {"of3": 4, "af2ig": 2, "dde": 1, "esm2": 1, "of3s1": 1}    # per-model job slots, so AF2-IG on the design poses is not starved by OF3
                for mk in ["of3", "af2ig"] + [m for m in tier2 if m != "af2ig"]:
                    its = [i for i in plan.get(mk, []) if (i["id"], mk) not in inflight]
                    bs = BATCH_FOR.get(mk, BATCH)
                    for j in range(0, len(its), bs):
                        running = sum(1 for k in futs if k[0] == mk)
                        if len(futs) >= a.max_jobs or running >= caps.get(mk, 1):
                            break
                        b = its[j:j + bs]
                        key = (mk, tuple(i["id"] for i in b))
                        futs[key] = ex.submit(run_job_safe, mk, b)
                        inflight |= {(i["id"], mk) for i in b}
            assemble()
            per, designs = build(load_state())
            plan_sz = sum(len(v) for v in plan.values())
            for line in funnel_lines(per, designs, load_state(), spent, cyc, plan_sz, n_top):
                log(line)
            log(f"jobs in flight {len(futs)} (max {a.max_jobs})")
            if stop and not futs:
                log(f"stopping ({'STOP file' if STOP_FILE.exists() else 'time limit'}); final spend ${spent:.2f}")
                break
            if stop and futs:
                log(f"stop requested, waiting for {len(futs)} jobs in flight")
        except Exception as e:      # never let one bad cycle kill the loop
            log(f"cycle error: {repr(e)[:300]}")
        time.sleep(a.interval)
    ex.shutdown(wait=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["score", "watch", "calib", "assemble", "reparse"])
    ap.add_argument("--models", default="of3,dde,esm2")
    ap.add_argument("--tier2", default="dde,esm2,of3s1,af2ig", help="watch: models for tier 2 (comma list; add af2ig to enable it)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-jobs", type=int, default=9, help="watch: max JapanFold jobs in flight from this stage")
    ap.add_argument("--no-decoy", action="store_true")
    ap.add_argument("--interval", type=int, default=90, help="watch: seconds between cycles (<= 300)")
    ap.add_argument("--only", default="", help="comma separated design_id prefixes to restrict to")
    a = ap.parse_args()
    a.interval = min(a.interval, 300)
    models = [m for m in a.models.split(",") if m]
    bad = [m for m in models if m not in MODELS]
    if bad:
        sys.exit(f"unknown models {bad}; known {list(MODELS)}")
    pre = tuple(x for x in a.only.split(",") if x)
    if a.cmd == "assemble":
        print(f"assembled {assemble()} designs -> score_of3.csv / score_orth.csv")
    elif a.cmd == "reparse":
        reparse()
        print(f"assembled {assemble()} designs")
    elif a.cmd == "calib":
        score(models, calib_items(), a.workers)
    elif a.cmd == "watch":
        watch(a)
    else:
        its = design_items(decoy=not a.no_decoy)
        if pre:
            its = [i for i in its if i["id"].startswith(pre)]
        score(models, its, a.workers)
        print(f"assembled {assemble()} designs")


if __name__ == "__main__":
    main()
