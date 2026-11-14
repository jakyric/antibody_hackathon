"""Extend cache/panel.csv with PLAbDab / PLAbDab-nano sequences (heavy-chain variable domains), deduplicated by Fv.
Sources (downloaded 2026-11-14 from opig.stats.ox.ac.uk/webapps/plabdab[-nano]/static/downloads):
  vhh_sequences.csv.gz       2.3k nanobodies (GenBank/patents/PDB), the 'plabdab_nano' arm of the abguard reference
  unpaired_sequences.csv.gz  all unpaired antibody chains (GenBank etc.); we keep chain == 'H'
  paired_sequences.csv.gz    paired antibodies; we keep heavy_sequence without ambiguity codes
Run: .venv/bin/python build_panel2.py   (keeps panel_v1.csv as backup; idempotent: re-running starts again from panel_v1.csv)
"""
import csv, gzip, shutil, sys, time
from multiprocessing import Pool
from pathlib import Path
import immunum

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
PL = CACHE / "plabdab"
STD = set("ACDEFGHIKLMNPQRSTVWY")
A = None

def init():
    global A
    A = immunum.Annotator(chains=["h"], scheme="imgt", min_confidence=0.4)

def work(item):
    src, cid, seq = item
    try:
        r = A.segment(seq)
        if r.error is not None:
            return None
        n = A.number(seq)
        fv = r.fr1 + r.cdr1 + r.fr2 + r.cdr2 + r.fr3 + r.cdr3 + r.fr4
        if len(fv) < 95:
            return None
        hall = "".join(n.numbering.get(p, "-") for p in ("42", "49", "50", "52"))
        vh = (src == "plabdab-nano") or (len(seq) <= 160 and sum(a in b for a, b in zip(hall, ("FY", "EQ", "R", "GFL"))) >= 3)
        return dict(id=cid, source=src, fv=fv, cdr3=r.cdr3, conf=round(n.confidence, 3), chain_len=len(seq), vhh_like=int(vh), n_chains=1)
    except Exception:
        return None

def items():
    with gzip.open(PL / "vhh_sequences.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            s = (r.get("sequence") or "").upper()
            if s and not set(s) - STD:
                yield ("plabdab-nano", "plnano:" + r["ID"], s)
    csv.field_size_limit(10**8)
    with gzip.open(PL / "unpaired_sequences.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            if r.get("chain") != "H":
                continue
            s = (r.get("GBSeq_sequence") or "").upper()
            if 90 <= len(s) <= 400 and not set(s) - STD:
                yield ("plabdab", "pl:" + r["GBSeq_accession-version"], s)
    with gzip.open(PL / "paired_sequences.csv.gz", "rt") as fh:
        for r in csv.DictReader(fh):
            s = (r.get("heavy_sequence") or "").upper()
            if 90 <= len(s) <= 400 and not set(s) - STD:
                yield ("plabdab", "plp:" + r["ID"], s)

if __name__ == "__main__":
    v1 = CACHE / "panel_v1.csv"
    if not v1.exists():
        shutil.copy(CACHE / "panel.csv", v1)
    base = list(csv.DictReader(open(v1)))
    seen = {r["fv"] for r in base}
    seen_seq, todo = set(), []
    for it in items():
        if it[2] in seen_seq:
            continue
        seen_seq.add(it[2]); todo.append(it)
    print("unique input sequences:", len(todo), flush=True)
    t0 = time.time(); out = []
    with Pool(16, initializer=init) as p:
        for k, r in enumerate(p.imap_unordered(work, todo, chunksize=500)):
            if r and r["fv"] not in seen:
                seen.add(r["fv"]); out.append(r)
            if k % 50000 == 0:
                print(f"  {k}/{len(todo)} {time.time()-t0:.0f}s kept {len(out)}", flush=True)
    rows = base + out
    with open(CACHE / "panel.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(base[0].keys())); w.writeheader(); w.writerows(rows)
    by = {}
    for r in rows:
        by[r["source"].split(":")[0]] = by.get(r["source"].split(":")[0], 0) + 1
    import json
    (CACHE / "panel_summary.json").write_text(json.dumps(dict(total=len(rows), by_source=by, vhh_like=sum(int(r["vhh_like"]) for r in rows)), indent=1))
    print("panel v2:", len(rows), by, f"{time.time()-t0:.0f}s")
