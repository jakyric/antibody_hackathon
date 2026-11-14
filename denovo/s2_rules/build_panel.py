"""Build the 'known antibodies' panel used for the whole-sequence identity gate.

Source 1: every protein chain in the wwPDB seqres file (length 100-300) that immunum numbers as an IGH variable domain
          (confidence >= 0.4). Deduplicated by variable-domain sequence. Includes Fab heavy chains, VH/VHH domains,
          scFvs (first/any heavy domain), therapeutics with PDB entries, nanobodies, HER2 antibodies with structures.
Source 2: UniProt reviewed human IGHV germline segments (and any camelid IGHV entries UniProt returns).
Source 3: hand-curated extra sequences in extra_known.fasta (literature HER2 VHHs etc.), each with its provenance.

Output: cache/panel.csv  (id, source, fv, cdr3, conf, chain_len, vhh_like)
Run:   .venv/bin/python build_panel.py
"""
import csv, gzip, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path
import immunum

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
A = immunum.Annotator(chains=["h"], scheme="imgt", min_confidence=0.4)


def fv_of(seq):
    r = A.segment(seq)
    if r.error is not None:
        return None
    n = A.number(seq)
    return ("".join([r.fr1, r.cdr1, r.fr2, r.cdr2, r.fr3, r.cdr3, r.fr4]), r.cdr3, n.confidence)


def pdb_panel(maxlen=300):
    seen, rows = {}, []
    t0, n = time.time(), 0
    cur_id = None
    with gzip.open(CACHE / "pdb_seqres.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split()
                cur_id, mol = parts[0], parts[1]
                cur_ok = mol == "mol:protein"
                continue
            if not cur_ok:
                continue
            s = line.strip()
            if not 100 <= len(s) <= maxlen:
                continue
            if set(s) - set("ACDEFGHIKLMNPQRSTVWY"):
                if set(s) - set("ACDEFGHIKLMNPQRSTVWYX"):
                    continue
            if s in seen:
                seen[s][1] += 1
                continue
            seen[s] = [cur_id, 1]
    print(f"unique PDB protein chains 100-{maxlen} aa: {len(seen)} ({time.time()-t0:.0f}s)", flush=True)
    by_fv = {}
    t0 = time.time()
    for i, (s, (cid, cnt)) in enumerate(seen.items()):
        # cheap prefilter: IG V domains have the canonical Cys..Trp-Gly-x-Gly signature
        if not re.search(r"W[GA].G", s) or s.count("C") < 2:
            continue
        r = fv_of(s)
        if r is None:
            continue
        fv, cdr3, conf = r
        if len(fv) < 95:
            continue
        if fv not in by_fv:
            by_fv[fv] = dict(id=cid, source="pdb", fv=fv, cdr3=cdr3, conf=round(conf, 3), chain_len=len(s),
                             vhh_like=int(len(s) <= 160), n_chains=cnt)
        else:
            by_fv[fv]["n_chains"] += cnt
            if len(s) <= 160:
                by_fv[fv]["vhh_like"] = 1
        if i % 20000 == 0:
            print(f"  numbered {i}/{len(seen)} ({time.time()-t0:.0f}s) panel so far {len(by_fv)}", flush=True)
    return list(by_fv.values())


def uniprot_ighv():
    rows = []
    for tax, label in [(9606, "human"), (9844, "llama"), (30538, "alpaca"), (9838, "camel")]:
        q = f"(gene:IGHV*) AND (organism_id:{tax})"
        url = "https://rest.uniprot.org/uniprotkb/search?format=fasta&size=500&query=" + urllib.parse.quote(q)
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            print("uniprot", label, "failed:", e)
            continue
        recs = [r for r in txt.split(">") if r.strip()]
        k = 0
        for rec in recs:
            head, *seq = rec.split("\n")
            s = "".join(seq).strip()
            r = fv_of(s) if len(s) >= 90 else None
            if r is None:
                continue
            fv, cdr3, conf = r
            acc = head.split("|")[1] if "|" in head else head.split()[0]
            rows.append(dict(id=f"uniprot:{acc}", source=f"uniprot-{label}", fv=fv, cdr3=cdr3, conf=round(conf, 3),
                             chain_len=len(s), vhh_like=0, n_chains=1))
            k += 1
        print(f"uniprot {label}: {k}/{len(recs)} numbered")
    return rows


def extra_fasta():
    f = HERE / "extra_known.fasta"
    rows = []
    if not f.exists():
        return rows
    for rec in f.read_text().split(">"):
        if not rec.strip():
            continue
        head, *seq = rec.strip().split("\n")
        s = "".join(seq).replace(" ", "")
        r = fv_of(s)
        if r is None:
            print("extra: cannot number", head[:40]); continue
        fv, cdr3, conf = r
        rows.append(dict(id=head.split()[0], source="literature:" + " ".join(head.split()[1:])[:60], fv=fv, cdr3=cdr3,
                         conf=round(conf, 3), chain_len=len(s), vhh_like=1, n_chains=1))
    print(f"extra_known.fasta: {len(rows)}")
    return rows


if __name__ == "__main__":
    rows = pdb_panel()
    print("PDB heavy-domain panel:", len(rows), "vhh_like:", sum(r["vhh_like"] for r in rows))
    rows += uniprot_ighv()
    rows += extra_fasta()
    with open(CACHE / "panel.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    by = {}
    for r in rows:
        by[r["source"].split(":")[0]] = by.get(r["source"].split(":")[0], 0) + 1
    (CACHE / "panel_summary.json").write_text(json.dumps(dict(total=len(rows), by_source=by,
                                                           vhh_like=sum(r["vhh_like"] for r in rows)), indent=1))
    print("panel written:", len(rows), by)
