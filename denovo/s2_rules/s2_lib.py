"""Shared helpers for stage S2 (numbering, identity panel, abguard-on-VM runner).

Numbering uses immunum==2.0.0, the exact pin of the abguard build on the team VM, so CDR boundaries and the
identity-to-trastuzumab number reproduce the VM tool. Identity is the abguard formula:
    identity = 1 - Levenshtein(Fv_a, Fv_b) / max(len(Fv_a), len(Fv_b))
on the numbered variable domain (FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4, tags/linkers/prefix/postfix removed).
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import immunum
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
sys.path.insert(0, str(HERE.parent))
import dn_common as dn  # noqa: E402

TRAS_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNT"
           "AYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")  # as in abguard.novelty (DrugBank DB00072)
TRAS_CDRH3_CORE = "WGGDGFYAMD"  # what the HER2 assay tables call CDRH3
ANN = immunum.Annotator(chains=["h", "k", "l"], scheme="imgt", min_confidence=0.1)  # abguard ANY_IG, floor 0.1
FR_REGIONS = ("fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4")
HALLMARK_POS = ("42", "49", "50", "52")  # IMGT; = Kabat 37, 44, 45, 47

import os
WORKERS = int(os.environ.get("S2_WORKERS", "8"))  # rapidfuzz threads; polite on a shared laptop
VM = os.environ.get("PIPELINE_VM_HOST", "")  # user@host of the shared VM that runs abguard; empty = VM abguard off
VM_KEY = os.environ.get("PIPELINE_VM_KEY", "")  # path of the ssh key
VM_DIR = "moritz/denovo_s2"
VM_ABGUARD = ".cache/uv/archive-v0/jsRWGsUwPztt-JwD/bin/abguard"  # newest build on the VM (2026-11-11 16:17)


def lev_identity(a: str, b: str) -> float:
    return 1.0 - Levenshtein.distance(a, b) / max(len(a), len(b))


def imgt_key(label: str):
    m = re.match(r"(\d+)([A-Z]?)", label)
    return int(m.group(1)), m.group(2)


@lru_cache(maxsize=None)
def analyse(seq: str) -> dict:
    """immunum numbering + segmentation. Returns {} fields with error text when the sequence cannot be numbered."""
    n = ANN.number(seq)
    s = ANN.segment(seq)
    out = dict(chain=n.chain, conf=float(n.confidence or 0.0), error=n.error_kind or n.error or "")
    if n.error_kind or s.error_kind:
        out.update(fv="", cdr1="", cdr2="", cdr3="", numbering={}, idx_label=[])
        return out
    fv = "".join(getattr(s, r) for r in FR_REGIONS)
    labels = list(n.numbering.keys())
    out.update(fv=fv, cdr1=s.cdr1, cdr2=s.cdr2, cdr3=s.cdr3, fr1=s.fr1, fr2=s.fr2, fr3=s.fr3, fr4=s.fr4,
               prefix=s.prefix, postfix=s.postfix, numbering=dict(n.numbering),
               idx_label=labels, q_start=n.query_start, q_end=n.query_end)
    return out


def imgt_region(label: str) -> str:
    """IMGT region of a position label (CDR1 27-38, CDR2 56-65, CDR3 105-117; the rest is framework)."""
    p = imgt_key(label)[0]
    if p <= 26:
        return "fr1"
    if p <= 38:
        return "cdr1"
    if p <= 55:
        return "fr2"
    if p <= 65:
        return "cdr2"
    if p <= 104:
        return "fr3"
    if p <= 117:
        return "cdr3"
    return "fr4"


TRAS_FV = analyse(TRAS_VH)["fv"]
TRAS_CDR3 = analyse(TRAS_VH)["cdr3"]  # IMGT 'SRWGGDGFYAMDY'


# ------------------------------------------------------------------ the known-antibody panel
class Panel:
    def __init__(self, path: Path | None = None):
        path = path or CACHE / "panel.csv"
        rows = list(csv.DictReader(open(path)))
        extra = HERE / "extra_known.fasta"
        self.ids = [r["id"] for r in rows]
        self.src = [r["source"] for r in rows]
        self.fv = [r["fv"] for r in rows]
        self.cdr3 = [r["cdr3"] for r in rows]
        self.vhh_like = np.array([int(r["vhh_like"]) for r in rows], dtype=bool)
        if TRAS_FV not in set(self.fv):  # trastuzumab VH is the reference of the abguard gate: always in the panel
            self.ids.append("trastuzumab_VH"); self.src.append("drugbank"); self.fv.append(TRAS_FV)
            self.cdr3.append(TRAS_CDR3); self.vhh_like = np.append(self.vhh_like, False)

    def __len__(self):
        return len(self.fv)

    def nearest(self, fvs: list[str], exclude_exact_self: bool = False, parent_fv: str | None = None):
        """For each query Fv: (best identity, id, source, best identity among VHH-like entries, id).
        `parent_fv`: ignore panel entries that are identical to the query itself only when exclude_exact_self."""
        M = process.cdist(fvs, self.fv, scorer=Levenshtein.normalized_similarity, workers=WORKERS)
        res = []
        for i in range(len(fvs)):
            row = M[i]
            if exclude_exact_self:
                row = np.where(row >= 0.999999, -1.0, row)
            j = int(np.argmax(row))
            jv = int(np.argmax(np.where(self.vhh_like, row, -1.0)))
            res.append(dict(identity_nearest_known=float(row[j]), nearest_known_id=self.ids[j], nearest_known_source=self.src[j],
                            identity_nearest_vhh=float(row[jv]), nearest_vhh_id=self.ids[jv]))
        return res

    def top(self, fv: str, k: int = 5):
        M = process.cdist([fv], self.fv, scorer=Levenshtein.normalized_similarity, workers=WORKERS)[0]
        idx = np.argsort(-M)[:k]
        return [(self.ids[i], self.src[i], round(float(M[i]), 3)) for i in idx]


# ------------------------------------------------------------------ abguard's own CDRH3 prior-art gate, replicated locally
_REF = None


def abguard_ref_cdrh3():
    global _REF
    if _REF is None:
        import pyarrow.parquet as pq
        t = pq.read_table(CACHE / "abguard_cdrh3.parquet").to_pandas()
        _REF = t
    return _REF


def cdrh3_min_dist(cdr3s: list[str]):
    """Min Levenshtein of each CDRH3 to the 610k CDRH3s that ship in abguard (+ trastuzumab). Chunked to bound memory."""
    ref = abguard_ref_cdrh3()
    refs = ref["cdrh3"].tolist()
    out = []
    for k in range(0, len(cdr3s), 16):
        chunk = cdr3s[k:k + 16]
        M = process.cdist(chunk, refs, scorer=Levenshtein.distance, workers=WORKERS, dtype=np.int32)
        for i, q in enumerate(chunk):
            j = int(np.argmin(M[i]))
            out.append((int(M[i, j]), refs[j], str(ref["source"].iloc[j]), str(ref["id"].iloc[j]) if "id" in ref else ""))
    return out


# ------------------------------------------------------------------ real abguard on the VM
def _ssh(cmd: str, timeout=300):
    return subprocess.run(["ssh", "-i", VM_KEY, "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", VM, cmd],
                          capture_output=True, text=True, timeout=timeout)


def abguard_vm(seqs: dict[str, str], use_cache: bool = True) -> dict[str, dict]:
    """Run the real abguard (`team-validate --all-liabilities`, bundled reference => novelty gate on) on the VM.
    seqs: {key: vhh sequence}. Returns {key: abguard row}. Results cached by sequence hash."""
    if not VM:
        raise RuntimeError("PIPELINE_VM_HOST is not set: VM abguard is off (use --no-vm)")
    cache_f = CACHE / "abguard_results.json"
    cache = json.loads(cache_f.read_text()) if (use_cache and cache_f.exists()) else {}
    h = {k: hashlib.sha1(s.encode()).hexdigest()[:12] for k, s in seqs.items()}
    todo = {h[k]: s for k, s in seqs.items() if h[k] not in cache}
    if todo:
        inp = CACHE / "abg_in.csv"
        with open(inp, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sequence_id", "kind", "heavy_sequence", "light_sequence", "design_description"])
            for hid, s in todo.items():
                w.writerow([hid, "vhh", s, "", "s2"])
        r = subprocess.run(["scp", "-q", "-i", VM_KEY, str(inp), f"{VM}:{VM_DIR}/abg_in.csv"], capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("scp up failed: " + r.stderr[:300])
        r = _ssh(f"cd ~/{VM_DIR} && ~/{VM_ABGUARD} team-validate abg_in.csv --all-liabilities -o abg_out.json 2>&1 | tail -3")
        if r.returncode:
            raise RuntimeError("abguard failed: " + (r.stdout + r.stderr)[-400:])
        r = subprocess.run(["scp", "-q", "-i", VM_KEY, f"{VM}:{VM_DIR}/abg_out.json", str(CACHE / "abg_out.json")],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("scp down failed: " + r.stderr[:300])
        for row in json.loads((CACHE / "abg_out.json").read_text()):
            row.pop("heavy_sequence", None)
            cache[row["sequence_id"]] = row
        cache_f.write_text(json.dumps(cache))
    return {k: cache[h[k]] for k in seqs}


# ------------------------------------------------------------------ watch-mode helpers
import datetime as _dt


def past_deadline(until_hhmm: str) -> bool:
    """True once local (CET on this machine) wall clock >= HH:MM today."""
    h, m = [int(x) for x in until_hhmm.split(":")]
    now = _dt.datetime.now()
    return (now.hour, now.minute) >= (h, m)


def stop_requested() -> bool:
    return (dn.DENOVO / "STOP").exists()
