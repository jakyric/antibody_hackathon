#!/usr/bin/env python
"""S1a: JapanFold rfd4 (RFD4-Proteina) VHH design against HER2 domain IV on a fixed framework. Resumable.

  denovo/.venv/bin/python denovo/s1_design/prep_inputs.py          # once: frameworks + placed inputs
  denovo/.venv/bin/python denovo/s1_design/run_design.py --reps 3 --workers 12     # wave 1
  denovo/.venv/bin/python denovo/s1_design/run_design.py --reps 10                 # wave 2 (re-uses wave-1 jobs)
  denovo/.venv/bin/python denovo/s1_design/run_design.py --dry-run                 # print the grid, cost nothing

One JOB = 8 designs (API max) for one grid cell: framework x hotspot set x CDR3 length range x pose x repeat.
Every job name is stable, so re-running never pays twice (dn_common.JapanFold.submit job cache, results/denovo/jobs/).
Outputs: results/denovo/designs_rfd4jf.csv (schema in CONTRACT.md), results/denovo/designs/<design_id>.cif,
         results/denovo/geom_rfd4jf.csv (per-design interface geometry; NOT a designs_* file on purpose).
"""
import argparse
import json
import math
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))
import dn_common as dn  # noqa: E402
import analyze as AN  # noqa: E402
import prep_inputs as PI  # noqa: E402

FWS = ["bcii10", "capla"]
HOTSETS = ["hs1_core", "hs2_broad", "hs3_alt_patch"]
C3 = PI.CDR3_RANGES                      # wave 1 (old): c3a 12-16, c3b 16-21; the designed block started at IMGT 105
CDR3_N = (15, 20)                        # wave 2 'c3n': designed block AFTER the fixed C104-A105-A106; total CDR3 incl. AA = 17-22
NEW_C3 = "c3n"
POSES = ["p0", "p1"]                     # p0 = trastuzumab-VH pose, p1 = perturbed pose (see make_pose)
RECIPE, NSTEPS, PER_JOB = "high/front_loaded", 60, 8
OUT = dn.OUT
CSV = OUT / "designs_rfd4jf.csv"
GEOM = OUT / "geom_rfd4jf.csv"
RUNS = HERE / "runs"
RUNS.mkdir(exist_ok=True)
LOCK = threading.Lock()
EPI = json.load(open(dn.TARGETS / "epitope.json"))


def hot_string(hs):
    """residue-level hotspot string for the API: A<n> for HER2 489-580, C<n> for 591-607."""
    nums = sorted({x["auth_seq"] for x in EPI["hotspot_sets"][hs]}, reverse=True)
    return ",".join(("A" if n <= 580 else "C") + str(n) for n in nums), nums


def make_pose(fw, pose, caa=False):
    """p0 = framework superposed on trastuzumab VH. pX = p0 rotated 10-14 deg about the CDR3-side pivot and shifted
    1.5 A away from HER2; skipped (falls back to p0) if the framework would clash with HER2 (CA-CA < 4 A)."""
    tag = "_caa" if (caa and fw == "capla") else ""      # capla needs the grafted A105/A106 (prep_caa.py); bcii10 already has them
    base = HERE / "inputs" / f"{fw}{tag}_on_her2.cif"
    if pose == "p0":
        return base
    path = HERE / "inputs" / f"{fw}{tag}_on_her2_{pose}.cif"
    if path.exists():
        return path
    seed = int(pose[1:]) * 7919 + (0 if fw == "bcii10" else 101)
    rng = np.random.default_rng(seed)
    ha = PI.read_atoms(base, "A", label_seq=False); hc = PI.read_atoms(base, "C", label_seq=False)
    vb = PI.read_atoms(base, "B", label_seq=False)
    her_ca = np.array([r["atoms"]["CA"] for r in ha + hc])
    epi = np.array([r["atoms"]["CA"] for r in ha + hc if r["num"] in AN.EPI_RES])
    vca = np.array([r["atoms"]["CA"] for r in vb])
    pivot = epi.mean(0)
    for attempt in range(50):
        ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
        ang = math.radians(rng.uniform(10, 14))
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        R = np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K
        away = vca.mean(0) - pivot; away /= np.linalg.norm(away)
        t = 1.5 * away
        new = [{"num": r["num"], "aa3": r["aa3"], "atoms": {a: tuple(R @ (np.array(x) - pivot) + pivot + t) for a, x in r["atoms"].items()}} for r in vb]
        nca = np.array([r["atoms"]["CA"] for r in new])
        dmin = np.linalg.norm(nca[:, None] - her_ca[None], axis=2).min()
        if dmin >= 4.0:
            PI.write_cif(path, [("A", ha), ("C", hc), ("B", new)], f"{fw}_{pose}")
            return path
    return base


def contig_for(fw_rec, fw, c3):
    L = PI.CDR_LEN[fw]
    if c3 == NEW_C3:
        fr = fw_rec["caa"]["contig_fr"]; lo3, hi3 = CDR3_N
    else:
        fr = fw_rec["contig_fr"]; lo3, hi3 = C3[c3]
    return (f"A489-580,/0,C591-607,/0,{fr[0]},{L['cdr1'][0]}-{L['cdr1'][1]},{fr[1]},{L['cdr2'][0]}-{L['cdr2'][1]},"
            f"{fr[2]},{lo3}-{hi3},{fr[3]}")


def grid(reps):
    cells = old_cells()
    jobs = []
    for ci, (fw, hs, c3, pose) in enumerate(cells):
        for rep in range(reps):
            seed = 1000 + ci * 100 + rep
            jobs.append({"name": f"s1a-{fw}-{hs}-{c3}-{pose}-r{rep}", "fw": fw, "hs": hs, "c3": c3, "pose": pose,
                         "rep": rep, "seed": seed})
    # interleave so that every cell progresses together and partial output is balanced
    jobs.sort(key=lambda j: (j["rep"], j["fw"], j["hs"], j["c3"], j["pose"]))
    return jobs


def load_rows():
    rows = dn.read_csv(CSV) if CSV.exists() else []
    geoms = dn.read_csv(GEOM) if GEOM.exists() else []
    return rows, geoms


STOP_AT = [None]
SLOTS = [threading.Semaphore(12)]   # max NEW jobs in flight on the shared account queue (set by --max-inflight)


def run_job(jf, job, fwrecs, state):
    if STOP_AT[0] and time.time() > STOP_AT[0] and dn.cached_job(job["name"]) is None:
        return {"job": job, "skipped": True}
    fw, hs = job["fw"], job["hs"]
    hot, hot_nums = hot_string(hs)
    rec = fwrecs[fw]
    caa = job["c3"] == NEW_C3
    fr_list = rec["caa"]["fr"] if caa else rec["fr"]
    contig = contig_for(rec, fw, job["c3"])
    inp = make_pose(fw, job["pose"], caa)
    body = {"protocol": "rfd4-binder", "name": job["name"], "structure": inp.read_text(), "contig": contig, "hotspots": hot,
            "params": {"num_designs": PER_JOB, "recipe": RECIPE, "nsteps": NSTEPS, "seed": job["seed"]}}
    held = dn.cached_job(job["name"]) is None
    if held:
        SLOTS[0].acquire()      # only jobs we actually submit count against the in-flight cap; already-submitted ones are just collected
    try:
        if held and STOP_AT[0] and time.time() > STOP_AT[0]:      # re-check after waiting for a slot: nothing new after --until
            return {"job": job, "skipped": True}
        for attempt in range(40):
            try:
                jid = jf.submit(job["name"], "/v1/designs", body)
                break
            except RuntimeError as e:
                msg = str(e)
                if "HTTP 429" in msg or "HTTP 503" in msg or ("HTTP 402" in msg and "still in progress" in msg):
                    # 402 here is NOT an empty wallet: the account's in-flight runs (S3's big OpenFold3 batches) reserve up to the balance
                    # (their max cost), so a new run waits until some finish. Back off and retry; a real out-of-credit 402 lacks this text.
                    time.sleep(min(60, 15 + 5 * attempt))
                    continue
                return {"job": job, "error": msg[:500]}
        else:
            return {"job": job, "error": "rate limited too long"}
        j = jf.wait(jid, max_s=2400, log=lambda *_: None)
    finally:
        if held:
            SLOTS[0].release()
    if j.get("status") != "succeeded":
        return {"job": job, "job_id": jid, "error": f"{j.get('status')}: {str(j.get('error'))[:400]}"}
    res = jf.results(jid)
    dest = RUNS / job["name"]
    jf.download(jid, dest)
    cost = 0.0
    try:
        cost = json.loads(jf.call("GET", f"/v1/jobs/{jid}/usage")).get("amount", 0.0)
    except Exception:
        pass
    rows, geoms = [], []
    for k, d in enumerate(res.get("designs", [])):
        cif = dest / d["structure"]
        did = f"{fw}-rfd4jf-s{job['seed']}-{k:03d}"
        g = AN.geometry(cif, fr_list, hot_nums)
        if not g.get("ok"):
            geoms.append({"design_id": did, "ok": 0, "why": g.get("why", ""), "job_id": jid})
            continue
        outcif = OUT / "designs" / f"{did}.cif"
        outcif.write_text(cif.read_text())
        rows.append({"design_id": did, "parent_id": "", "generator": "rfd4jf", "framework": fw, "seed": job["seed"],
                     "hotspot_set": hs, "vhh_seq": g["vhh_seq"], "cdr1": g["cdr1"], "cdr2": g["cdr2"],
                     "cdr3": (rec["caa"]["prefix"] if caa else "") + g["cdr3"],   # always the IMGT 105-117 style block, incl. the fixed AA
                     "cdr_scheme": "IMGT-ext 24-40/56-65/105-117", "design_cif": str(outcif.relative_to(dn.ROOT)), "job_id": jid,
                     "params_json": json.dumps({"contig": contig, "hotspots": hot, "recipe": RECIPE, "nsteps": NSTEPS,
                                                "c3_range": CDR3_N if caa else C3[job["c3"]], "c3_label": job["c3"],
                                                "caa_fixed": int(caa), "caa_grafted": bool(caa and rec["caa"].get("grafted")),
                                                "cdr3_designed_len": len(g["cdr3"]), "cdr3_total_len": len(g["cdr3"]) + (2 if caa else 0),
                                                "pose": job["pose"], "api_design": d["id"],
                                                "motif_rmsd": d.get("motif_rmsd"), "num_chains": d.get("num_chains")}),
                     "created_utc": dn.utcnow()})
        geoms.append({"design_id": did, "framework": fw, "hotspot_set": hs, "c3_range": job["c3"], "pose": job["pose"],
                      "motif_rmsd": d.get("motif_rmsd"), "job_id": jid, "caa_fixed": int(caa), "cdr3_designed_len": len(g["cdr3"]),
                      **{k2: v for k2, v in g.items() if k2 not in ("vhh_seq",)}})
    return {"job": job, "job_id": jid, "rows": rows, "geoms": geoms, "cost": cost}


def parse_hhmm(x):
    h, m = map(int, x.split(":"))
    t = time.localtime()
    return time.mktime((t.tm_year, t.tm_mon, t.tm_mday, h, m, 0, 0, 0, -1))


def old_cells():
    return [(fw, hs, c3, pose) for fw in FWS for hs in HOTSETS for c3 in C3 for pose in POSES]


def all_cells():
    """cells that NEW waves may use: both frameworks x 3 hotspot sets x c3n (CAA fixed, total CDR3 17-22) x 2 poses"""
    return [(fw, hs, NEW_C3, pose) for fw in FWS for hs in HOTSETS for pose in POSES]


def job_name(cell, rep):
    fw, hs, c3, pose = cell
    return f"{'s1a2' if c3 == NEW_C3 else 's1a'}-{fw}-{hs}-{c3}-{pose}-r{rep}"


def existing_reps(cell):
    fw, hs, c3, pose = cell
    pre = "s1a2" if c3 == NEW_C3 else "s1a"
    reps = [int(p.stem.rsplit("-r", 1)[1]) for p in dn.JOBS.glob(f"{pre}-{fw}-{hs}-{c3}-{pose}-r*.json")]
    return max(reps) + 1 if reps else 0


def recover_jobs(done_ids):
    """jobs that were submitted (cached id) but whose results never reached the CSV (process stopped/restarted): collect them, pay nothing.
    Handles both name styles: wave 1 `s1a-<fw>-<hs>-<c3a|c3b>-<pose>-r<n>` and wave 2 `s1a2-<fw>-<hs>-c3n-<pose>-r<n>`."""
    old, new = old_cells(), all_cells()
    out = []
    for f in sorted(list(dn.JOBS.glob("s1a-*.json")) + list(dn.JOBS.glob("s1a2-*.json"))):
        name = f.stem
        try:
            pre, fw, hs, c3, pose, rep = name.split("-")
            rep = int(rep[1:]); cell = (fw, hs, c3, pose)
            seed = (1000 + old.index(cell) * 100 + rep) if pre == "s1a" else (100000 + new.index(cell) * 100 + rep)
        except (ValueError, IndexError):
            continue
        if json.loads(f.read_text())["job_id"] in done_ids:
            continue
        out.append({"name": name, "fw": fw, "hs": hs, "c3": c3, "pose": pose, "rep": rep, "seed": seed})
    return out


def plan_jobs(n_jobs, log):
    import feedback as FB
    alloc, weights, why, basis = FB.plan(n_jobs, all_cells())
    log(f"plan basis: {basis}")
    jobs = []
    for cell, k in alloc.items():
        fw, hs, c3, pose = cell
        base = existing_reps(cell)
        ci = all_cells().index(cell)
        for rep in range(base, base + k):
            jobs.append({"name": job_name(cell, rep), "fw": fw, "hs": hs, "c3": c3, "pose": pose, "rep": rep,
                         "seed": 100000 + ci * 100 + rep})
    with open(HERE / "plan_log.jsonl", "a") as fh:
        fh.write(json.dumps({"ts": dn.utcnow(), "basis": basis, "alloc": {"|".join(k): v for k, v in alloc.items()},
                             "why": {"|".join(k): v for k, v in why.items()}}) + "\n")
    jobs.sort(key=lambda j: (j["rep"], j["fw"], j["hs"], j["c3"], j["pose"]))
    return jobs


def run_wave(jobs, jf, fwrecs, rows, geoms, workers, budget, tag=""):
    spent, nfail, nskip = 0.0, 0, 0
    t0 = time.time()
    ledger = HERE / "spend_ledger.jsonl"
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(run_job, jf, j, fwrecs, None) for j in jobs]
        for f in as_completed(futs):
            r = f.result()
            with LOCK:
                if r.get("skipped"):
                    nskip += 1
                    continue
                if "error" in r:
                    nfail += 1
                    print(f"FAILED {r['job']['name']}: {r['error'][:200]}", flush=True)
                    with open(HERE / "failed_jobs.jsonl", "a") as fh:
                        fh.write(json.dumps({"name": r["job"]["name"], "job_id": r.get("job_id"), "error": r["error"]}) + "\n")
                    continue
                rows.extend(r["rows"]); geoms.extend(r["geoms"])
                spent += r["cost"]
                dn.write_csv(CSV, rows, dn.DESIGN_COLS)      # flushed after EVERY finished job
                if geoms:
                    cols = sorted({k for g in geoms for k in g}, key=lambda k: (k not in ("design_id", "framework", "hotspot_set", "c3_range", "pose"), k))
                    dn.write_csv(GEOM, geoms, cols)
                with open(ledger, "a") as fh:
                    fh.write(json.dumps({"job": r["job"]["name"], "job_id": r["job_id"], "cost": r["cost"], "n": len(r["rows"])}) + "\n")
                print(f"{time.strftime('%H:%M:%S')} {tag}{r['job']['name']}: +{len(r['rows'])}; total {len(rows)}; wave spend ${spent:.3f}; fail {nfail}", flush=True)
                if sum(json.loads(x)["cost"] for x in open(ledger)) > budget:
                    print("budget reached; cancelling remaining jobs", flush=True)
                    STOP_AT[0] = time.time() - 1
    return spent, nfail, nskip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=3, help="uniform grid: repeats per cell (8 designs per repeat)")
    ap.add_argument("--auto", action="store_true", help="loop waves, re-planning cells from S3 scores (feedback.py)")
    ap.add_argument("--until", default=None, help="HH:MM local: submit nothing new after this time (in-flight jobs finish)")
    ap.add_argument("--jobs-per-wave", type=int, default=36)
    ap.add_argument("--workers", type=int, default=40, help="threads (cheap: they mostly wait); in-flight NEW jobs are capped by --max-inflight")
    ap.add_argument("--max-inflight", type=int, default=12, help="max new jobs queued/running on the account at once")
    ap.add_argument("--budget", type=float, default=14.0, help="total USD of rfd4 usage across all runs (spend_ledger.jsonl)")
    ap.add_argument("--only-fw", default=None); ap.add_argument("--only-hs", default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dn.wait_targets(60)
    fwrecs = {f: json.load(open(HERE / "frameworks" / f"{f}.json")) for f in FWS}
    SLOTS[0] = threading.Semaphore(a.max_inflight)
    if a.until:
        STOP_AT[0] = parse_hhmm(a.until)
    rows, geoms = load_rows()
    for fw in FWS:
        for pose in POSES:
            make_pose(fw, pose); make_pose(fw, pose, True)   # build perturbed-pose inputs once, before threads start
    if a.dry_run:
        jobs = grid(a.reps) if not a.auto else plan_jobs(a.jobs_per_wave, print)
        print(f"{len(jobs)} jobs ({len(jobs) * PER_JOB} designs)")
        for j in jobs[:4]:
            print(j["name"], contig_for(fwrecs[j["fw"]], j["fw"], j["c3"]), hot_string(j["hs"])[0])
        return
    jf = dn.JapanFold()
    if not a.auto:
        jobs = grid(a.reps)
        if a.only_fw:
            jobs = [j for j in jobs if j["fw"] == a.only_fw]
        if a.only_hs:
            jobs = [j for j in jobs if j["hs"] == a.only_hs]
        done = {r["job_id"] for r in rows}
        todo = [j for j in jobs if dn.cached_job(j["name"]) is None or dn.cached_job(j["name"]) not in done]
        print(f"grid {len(jobs)} jobs; to run {len(todo)}", flush=True)
        run_wave(todo, jf, fwrecs, rows, geoms, a.workers, a.budget)
    else:
        w = 0
        while not (STOP_AT[0] and time.time() > STOP_AT[0]):
            jobs = plan_jobs(a.jobs_per_wave, lambda m: print(m, flush=True))
            done = {r["job_id"] for r in rows} | {g.get("job_id") for g in geoms}
            jobs = recover_jobs(done) + [j for j in jobs if dn.cached_job(j["name"]) is None or dn.cached_job(j["name"]) not in done]
            if not jobs:
                print("nothing to run", flush=True); break
            print(f"wave {w}: {len(jobs)} jobs", flush=True)
            run_wave(jobs, jf, fwrecs, rows, geoms, a.workers, a.budget, tag=f"w{w} ")
            w += 1
    print(f"done: {len(rows)} designs, balance ${jf.credits():.2f}", flush=True)


if __name__ == "__main__":
    main()
