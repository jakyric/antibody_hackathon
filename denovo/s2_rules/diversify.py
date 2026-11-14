"""S2 framework diversification: cut whole-Fv identity to every known antibody below the 70% novelty line while keeping the
design a believable VHH. Substitution only (no indels), CDR1-3 untouched, so the design backbone/CDR geometry stays valid.

Why: a fixed natural framework with brand-new CDRs stays ~75-90% identical to its parent, and the parent is a known antibody.
abguard on the VM classes a design 'de_novo' only if its whole Fv is < 70% identical to trastuzumab VH; the jury text also
says < 70% to ANY known antibody. We therefore drive identity_nearest_known (PDB heavy domains + human IGHV germlines +
trastuzumab) below --target (default 0.66, i.e. ~4 points of margin for sequences missing from our panel).

Method (greedy identity descent under naturalness + structure constraints):
  1. Mutable positions = IMGT framework positions minus a protected set (CDRs, Cys23/Cys104, Trp41, FR2 hallmarks 42/49/50/52,
     Trp118-Gly119-Gly121, CDR-flanking 26/39/55/66). Buried positions (rel. SASA < 0.20 in 3EAK) accept only
     conservative swaps inside their hydrophobic class.
  2. Allowed residues at a position = residues seen at that IMGT position in >= 1.5% of the panel (--min-freq) (PDB VHH/VH domains weighted
     x2 for VHH-like, + human germlines), never Cys, never creating an N-glyc sequon or new Met/Trp in the framework.
     This keeps every substitution a naturally occurring one (no 'unnatural' MPNN-style chimeras).
  3. Each step picks the substitution that lowers the mean of the top-10 nearest-known identities (and the top-1) most, with a
     bonus for (a) naturalness (log-odds of the new vs old residue) (b) ProteinMPNN vote on the VHH backbone (3EAK, hosted
     Nebius run, 96 samples T=0.3) (c) germline-likeness: residue found in a *different* human IGHV germline chosen per variant.
     Choice is sampled from the best 4 (softmax) so the K variants of one design differ. Stops at >= --min-mut substitutions
     AND identity_nearest_known <= --target (hard cap --max-mut).
  4. Optional --esm: Nebius ESM-2 pseudo-log-likelihood filter (variant must not lose more than a margin vs its parent).

Usage:
  .venv/bin/python denovo/s2_rules/diversify.py --in results/denovo/designs_rfd4jf.csv --k 6 --out results/denovo/designs_fwdiv.csv
  .venv/bin/python denovo/s2_rules/diversify.py --standin            # scratch test on stand-in designs (never writes designs_*.csv)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import s2_lib as L  # noqa: E402

dn = L.dn
CACHE = L.CACHE
AA20 = "ACDEFGHIKLMNPQRSTVWY"
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201, K=236, M=224, F=240, P=159,
              S=155, T=172, W=285, Y=263, V=174)  # Tien et al. 2013, theoretical
HYD_CLASS = [set("AVILM"), set("FYW"), set("ILVM")]
PROTECTED = {"23", "41", "42", "49", "50", "52", "104", "118", "119", "121", "26", "39", "55", "66"}


# ------------------------------------------------------------------ priors
def build_pssm(panel: L.Panel) -> dict:
    """Per-IMGT-position residue counts over the panel: counts[label] = {aa: weight}; VHH-like entries weigh 3, others 1."""
    f = CACHE / "pssm.json"
    if f.exists():
        return json.loads(f.read_text())
    cnt: dict[str, Counter] = defaultdict(Counter)
    for fv, vh in zip(panel.fv, panel.vhh_like):
        n = L.ANN.number(fv)  # not L.analyse(): its lru_cache would hold 180k numbering dicts
        if n.error_kind or n.error:
            continue
        w = 3.0 if vh else 1.0
        for lab, aa in n.numbering.items():
            cnt[lab][aa] += w
    out = {k: dict(v) for k, v in cnt.items()}
    f.write_text(json.dumps(out))
    return out


def freq(pssm: dict, lab: str, aa: str) -> float:
    c = pssm.get(lab, {})
    tot = sum(c.values()) or 1.0
    return c.get(aa, 0.0) / tot


def template_burial() -> dict[str, float]:
    """Relative SASA per IMGT label from the 3EAK crystal structure (NbBCII10-FGLA, 1.95 A), chain A, VHH alone."""
    f = CACHE / "burial_3eak.json"
    if f.exists():
        return json.loads(f.read_text())
    from Bio.PDB import PDBParser, ShrakeRupley
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    s = PDBParser(QUIET=True).get_structure("x", CACHE / "structs" / "3EAK.pdb")
    ch = s[0]["A"]
    ShrakeRupley().compute(ch, level="R")
    res = [r for r in ch if r.id[0] == " "]
    seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in res)
    a = L.analyse(seq)
    out = {}
    for i, lab in enumerate(a["idx_label"]):
        r = res[a["q_start"] + i]
        aa = seq[a["q_start"] + i]
        out[lab] = round(r.sasa / MAXASA[aa], 3)
    f.write_text(json.dumps(out))
    return out


def load_mpnn_votes(path: Path | None) -> dict[str, Counter]:
    """ProteinMPNN samples (hosted Nebius run on the 3EAK backbone) -> per-IMGT-label residue votes."""
    if path is None or not Path(path).exists():
        return {}
    d = json.loads(Path(path).read_text())
    recs = [r for r in d["mfasta"].split(">") if r.strip()]
    ref = recs[0].split("\n", 1)[1].replace("\n", "")
    a = L.analyse(ref)
    labels = a["idx_label"]
    votes: dict[str, Counter] = defaultdict(Counter)
    for r in recs[1:]:
        s = r.split("\n", 1)[1].replace("\n", "")
        if len(s) != len(ref):
            continue
        for i, lab in enumerate(labels):
            votes[lab][s[a["q_start"] + i]] += 1
    return votes


def germline_label_maps(panel: L.Panel) -> dict[str, dict[str, str]]:
    out = {}
    for i, src in enumerate(panel.src):
        if src.startswith("uniprot-human"):
            a = L.analyse(panel.fv[i])
            if a.get("numbering") and ("IGHV3" in panel.ids[i] or True):
                out[panel.ids[i]] = a["numbering"]
    return out


# ------------------------------------------------------------------ core
def sequon_or_bad(seq: str) -> bool:
    return bool(re.search(r"N[^P][ST]", seq))


class Diversifier:
    def __init__(self, panel: L.Panel | None = None, mpnn_json: Path | None = None, min_freq: float = 0.015, mask_p: float = 0.25):
        self.min_freq, self.mask_p = min_freq, mask_p
        self.panel = panel or L.Panel()
        self.pssm = build_pssm(self.panel)
        self.burial = template_burial()
        self.votes = load_mpnn_votes(mpnn_json)
        self.germ = germline_label_maps(self.panel)
        self.pfv = self.panel.fv

    # -- positions
    def mutable_positions(self, a: dict) -> list[tuple[int, str]]:
        """(index into Fv string, IMGT label) of mutable framework positions."""
        out = []
        for i, lab in enumerate(a["idx_label"]):
            if L.imgt_region(lab) in ("cdr1", "cdr2", "cdr3") or lab in PROTECTED:
                continue
            out.append((i, lab))
        return out

    def allowed(self, lab: str, aa: str) -> list[str]:
        buried = self.burial.get(lab, 0.5) < 0.20
        res = []
        for b in AA20:
            if b == aa or b in "CMW":
                continue
            if freq(self.pssm, lab, b) < self.min_freq:
                continue
            if buried and not any(aa in g and b in g for g in HYD_CLASS):
                continue
            res.append(b)
        return res

    def nat(self, lab: str, a: str, b: str) -> float:
        c = self.pssm.get(lab, {})
        return math.log((c.get(b, 0) + 1) / (c.get(a, 0) + 1))

    def diversify_one(self, seq: str, seed: int, target: float, min_mut: int, max_mut: int, germline_pool: list[str],
                      mode: str = "deep") -> dict | None:
        """mode 'deep': descend identity to the nearest known antibody (whole panel) to <= target.
        mode 'lite': the fewest natural substitutions (<= max_mut) that bring identity to trastuzumab VH to <= target."""
        lite = mode == "lite"
        rng = random.Random(seed)
        a = L.analyse(seq)
        fv0 = a["fv"]
        pre, post = a.get("prefix", ""), a.get("postfix", "")
        assert "".join(a["numbering"].values()) == fv0, "numbering/Fv mismatch"
        pos = self.mutable_positions(a)
        # per-variant random mask: ~25% of the mutable positions are off-limits for this variant, which is what makes the K
        # variants of one design differ from each other instead of converging on the same greedy optimum
        pos_idx = {lab: i for i, lab in pos if rng.random() > (0.12 if lite else self.mask_p)}
        gl_id = rng.choice(germline_pool) if germline_pool else None
        gl = self.germ.get(gl_id, {}) if gl_id else {}
        cur = list(fv0)
        muts: list[tuple[str, str, str]] = []  # (label, old, new)
        used = set()
        P_fv = self.pfv
        N = None
        for step in range(max_mut):
            curfv = "".join(cur)
            if lite:
                cur_top1 = cur_top10 = L.lev_identity(curfv, L.TRAS_FV)
                if cur_top1 <= target:
                    break
                N = [0]
            elif N is None or step % 4 == 0:
                M = L.process.cdist([curfv], P_fv, scorer=L.Levenshtein.normalized_similarity, workers=L.WORKERS)[0]
                N = np.argsort(-M)[:300]
                cur_top1 = float(M[N[0]])
                cur_top10 = float(np.mean(np.sort(M)[-10:]))
                if cur_top1 <= target and len(muts) >= min_mut:
                    break
            neigh = [L.TRAS_FV] if lite else [P_fv[j] for j in N]
            cands = []
            for lab, i in pos_idx.items():
                if i in used:
                    continue
                old = cur[i]
                for b in self.allowed(lab, fv0[i]):
                    if b == old:
                        continue
                    cands.append((i, lab, old, b))
            if not cands:
                break
            seqs = []
            for i, lab, old, b in cands:
                s = cur.copy(); s[i] = b
                seqs.append("".join(s))
            S = L.process.cdist(seqs, neigh, scorer=L.Levenshtein.normalized_similarity, workers=L.WORKERS)
            top1 = S.max(axis=1)
            top10 = np.sort(S, axis=1)[:, -10:].mean(axis=1)
            scores = []
            for k, (i, lab, old, b) in enumerate(cands):
                full = pre + seqs[k] + post
                if sequon_or_bad(full) and not sequon_or_bad(pre + "".join(cur) + post):
                    scores.append(-9); continue
                v = self.votes.get(lab)
                vote = (v.get(b, 0) / sum(v.values())) if v else 0.0
                if lite and cur_top1 - top1[k] <= 0:
                    scores.append(-9); continue  # lite only takes substitutions that actually lower identity to trastuzumab VH
                sc = ((cur_top10 - top10[k]) + (cur_top1 - top1[k])
                      + 0.003 * self.nat(lab, fv0[i], b)
                      + 0.004 * vote
                      + (0.004 if gl.get(lab) == b else 0.0))
                scores.append(sc)
            order = np.argsort(scores)[::-1][:8]
            w = np.array([math.exp((scores[o] - scores[order[0]]) / 0.004) for o in order])
            pick = int(rng.choices(list(order), weights=w)[0])
            i, lab, old, b = cands[pick]
            cur[i] = b
            used.add(i)
            muts.append((lab, old, b))
            cur_top1, cur_top10 = float(top1[pick]), float(top10[pick])
            if lite:
                continue
            if cur_top1 <= target and len(muts) >= min_mut:  # confirm against the whole panel, not just the 300 neighbours
                full = L.process.cdist(["".join(cur)], P_fv, scorer=L.Levenshtein.normalized_similarity, workers=L.WORKERS)[0]
                if float(full.max()) <= target:
                    break
        newfv = "".join(cur)
        new_seq = pre + newfv + post
        b = L.analyse(new_seq)
        if b["fv"] != newfv or (b["cdr1"], b["cdr2"], b["cdr3"]) != (a["cdr1"], a["cdr2"], a["cdr3"]) or b["chain"] != "H":
            return None  # numbering shifted: reject (never silently change a CDR)
        if not muts:
            return None
        return dict(seq=new_seq, muts=muts, germline=gl_id or "")

    def run_design(self, d: dict, k: int, targets: list[float], min_mut: int, max_mut: int, seed0: int = 0, min_ham: int = 8,
                   tier: str = "deep") -> list[dict]:
        seq = d["vhh_seq"]
        a = L.analyse(seq)
        if not a.get("fv"):
            return []
        pool = [g for g in self.germ if g.startswith("uniprot:")]
        # germline pool: human IGHV3 family members closest to the parent framework (natural targets), excluding its single nearest
        fr = a["fv"]
        sims = L.process.cdist([fr], [self.pfv[self.panel.ids.index(g)] for g in pool], scorer=L.Levenshtein.normalized_similarity, workers=L.WORKERS)[0]
        ranked = [g for _, g in sorted(zip(sims, pool), reverse=True)]
        gpool = ranked[1:13]
        out, seen = [], {seq}
        for v in range(k * 6):  # retries: variants must differ pairwise by >= min_ham positions
            if len(out) >= k:
                break
            tgt = targets[len(out) % len(targets)]
            r = self.diversify_one(seq, seed=int(hashlib.sha1(f'{d["design_id"]}|{seed0}|{tier}|{v}'.encode()).hexdigest()[:8], 16),
                                   target=tgt, min_mut=min_mut, max_mut=max_mut, germline_pool=gpool, mode=tier)
            if r is None or r["seq"] in seen:
                continue
            if any(sum(x != y for x, y in zip(r["seq"], o["seq"])) < min_ham for o in out):
                continue
            r["target"] = tgt
            seen.add(r["seq"])
            out.append(r)
        return out


def fv_identity(a: str, b: str) -> float:
    return L.lev_identity(a, b)


def make_rows(div: Diversifier, designs: list[dict], k: int, targets: list[float], min_mut: int, max_mut: int,
              tier: str = "deep"):
    """tier 'deep': targets = nearest-known identity targets (cycled), up to max_mut substitutions.
    tier 'lite': targets = [trastuzumab-VH identity target], max_mut = max substitutions (<= 8), min_mut ignored."""
    rows, report = [], []
    suffix = "fl" if tier == "lite" else "fd"
    for d in designs:
        pa = L.analyse(d["vhh_seq"])
        before = div.panel.nearest([pa["fv"]])[0]
        id_tr0 = L.lev_identity(pa["fv"], L.TRAS_FV)
        variants = div.run_design(d, k, targets, 1 if tier == "lite" else min_mut, max_mut, min_ham=3 if tier == "lite" else 8, tier=tier)
        for n, r in enumerate(variants, 1):
            vid = dn.safe_id(f'{d["design_id"][:40]}-{suffix}{n}')
            va = L.analyse(r["seq"])
            after = div.panel.nearest([va["fv"]])[0]
            rows.append(dict(design_id=vid, parent_id=d["design_id"], generator="fwdiv", framework=d.get("framework", ""),
                             seed=d.get("seed", ""), hotspot_set=d.get("hotspot_set", ""), vhh_seq=r["seq"], cdr1=va["cdr1"],
                             cdr2=va["cdr2"], cdr3=va["cdr3"], cdr_scheme="imgt", design_cif=d.get("design_cif", ""),
                             job_id="", params_json=json.dumps(dict(tier=tier, n_mut=len(r["muts"]), germline=r["germline"], target=r["target"],
                                                                    muts=",".join(f"{o}{l}{n_}" for l, o, n_ in r["muts"]))),
                             created_utc=dn.utcnow()))
            report.append(dict(design_id=vid, parent_id=d["design_id"], tier=tier, n_mut=len(r["muts"]), germline=r["germline"],
                               target=r["target"],
                               identity_to_parent=round(fv_identity(pa["fv"], va["fv"]), 3),
                               nearest_before=round(before["identity_nearest_known"], 3), nearest_before_id=before["nearest_known_id"],
                               nearest_after=round(after["identity_nearest_known"], 3), nearest_after_id=after["nearest_known_id"],
                               trast_before=round(id_tr0, 3), trast_after=round(L.lev_identity(va["fv"], L.TRAS_FV), 3),
                               ok=int(after["identity_nearest_known"] < 0.70),
                               lite_ok=int(L.lev_identity(va["fv"], L.TRAS_FV) <= 0.66), muts=",".join(f"{o}{l}{n_}" for l, o, n_ in r["muts"])))
    return rows, report



# ------------------------------------------------------------------ watch mode: parents come from S3's scores
STATE = HERE / "state_fwdiv.json"


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def select_parents(top_n: int = 40, per_gen: int = 3, skip_ident_below: float = 0.68, unscored_n: int = 0) -> list[tuple[dict, float, str]]:
    """Parents ranked by S3's OF3 consensus: top `top_n` overall + best `per_gen` of every generator, restricted to designs
    that pass every hard rule except identity. Returns [(design_row, score, why)] best first."""
    of3 = {r["design_id"]: r for r in dn.read_csv(dn.OUT / "score_of3.csv")} if (dn.OUT / "score_of3.csv").exists() else {}
    orth = {r["design_id"]: r for r in dn.read_csv(dn.OUT / "score_orth.csv")} if (dn.OUT / "score_orth.csv").exists() else {}
    designs = {d["design_id"]: d for d in dn.load_designs(include_controls=False)
               if d.get("generator") not in ("fwdiv", "control") and d.get("vhh_seq") and not d["design_id"].startswith("ctl-")}
    import rules as R
    rules = {r["design_id"]: r for r in R.from_csv_rows(dn.read_csv(dn.OUT / "rules.csv"))} if (dn.OUT / "rules.csv").exists() else {}
    cand, unscored = [], []
    for did, d in designs.items():
        o3, ob = of3.get(did), orth.get(did)
        if not o3:
            unscored.append(d)  # not scored yet: only used by the time-boxed fallback below
            continue
        sc = _f((ob or {}).get("consensus_rank_score"))
        if sc is None:
            sc = _f((ob or {}).get("of3_pct"))
        if sc is None:
            sc = _f(o3.get("of3_iptm"))
        if sc is None:
            continue
        r = rules.get(did) or R.run([d], use_vm=False)[0]
        ok = (int(r.get("numbering_ok", 0)) and int(r.get("canon_cys_ok", 0)) and not int(r.get("unpaired_cys", 1))
              and int(r.get("len_ok", 0)) and not int(r.get("invalid_aa", 1)) and int(r.get("cdrh3_novel_local", 0)))
        if not ok:
            continue
        idn = _f(r.get("identity_nearest_known"))
        if idn is not None and idn < skip_ident_below:
            continue  # already far enough from everything known: nothing to fix
        cand.append((d, sc))
    cand.sort(key=lambda x: -x[1])
    chosen, why = {}, {}
    if unscored_n and len([1 for _ in cand]) < unscored_n:
        # Fallback so the downstream stages (refold of variants, NanoMelt, ranking) are exercised on real rows even when S3
        # is behind: unscored designs that pass every hard rule except identity, one per CDRH3, fewest substitutions needed first.
        seen3, extra = set(), []
        for d in unscored:
            r = rules.get(d["design_id"]) or R.run([d], use_vm=False)[0]
            ok = (int(r.get("numbering_ok", 0)) and int(r.get("canon_cys_ok", 0)) and not int(r.get("unpaired_cys", 1))
                  and int(r.get("len_ok", 0)) and not int(r.get("invalid_aa", 1)) and int(r.get("cdrh3_novel_local", 0)))
            idn = _f(r.get("identity_nearest_known"))
            if not ok or (idn is not None and idn < skip_ident_below) or r.get("cdrh3") in seen3:
                continue
            seen3.add(r.get("cdrh3")); extra.append((d, -1.0 - float(r.get("identity_to_trastuzumab") or 1.0)))
        extra.sort(key=lambda x: -x[1])
        for d, sc in extra[: unscored_n - len(cand)]:
            chosen[d["design_id"]] = (d, sc); why[d["design_id"]] = "unscored-fallback"

    for d, sc in cand[:top_n]:
        chosen[d["design_id"]] = (d, sc); why[d["design_id"]] = f"top{top_n}"
    gens = {}
    for d, sc in cand:
        gens.setdefault(d.get("generator", ""), []).append((d, sc))
    for g, lst in gens.items():
        for d, sc in lst[:per_gen]:
            if d["design_id"] not in chosen:
                chosen[d["design_id"]] = (d, sc); why[d["design_id"]] = f"best-of-{g}"
    return sorted(((d, sc, why[d["design_id"]]) for d, sc in chosen.values()), key=lambda x: -x[1])


def repair_canonical(seq: str) -> tuple[str, list[str]] | None:
    """Restore the canonical VHH disulfide (Cys at IMGT 23 and 104) and Trp118 when a generator dropped them
    (BoltzGen's nanobody protocol never writes Cys). CDRs untouched. Returns (new_seq, [edits]) or None if not repairable."""
    a = L.analyse(seq)
    if a["chain"] != "H" or not a.get("numbering"):
        return None
    want = {"23": "C", "104": "C", "118": "W"}
    s = list(seq); edits = []
    for i, lab in enumerate(a["idx_label"]):
        if lab in want and a["numbering"][lab] != want[lab]:
            j = a["q_start"] + i
            edits.append(f"{s[j]}{j + 1}{want[lab]}")
            s[j] = want[lab]
    if not edits:
        return None
    new = "".join(s)
    b = L.analyse(new)
    if b["chain"] != "H" or (b["cdr1"], b["cdr2"], b["cdr3"]) != (a["cdr1"], a["cdr2"], a["cdr3"]) or new.count("C") % 2:
        return None
    return new, edits


def repair_pass(state: dict, rows: list[dict], known_seqs: set, cap: int = 60, log=print) -> list[dict]:
    """One '-cr' variant (canonical Cys/Trp restored) for every design that fails ONLY the canonical-disulfide rule."""
    import rules as R
    out = []
    done = state.setdefault("repaired", {})
    for d in dn.load_designs(include_controls=False):
        did = d["design_id"]
        if d.get("generator") in ("fwdiv", "control") or did.startswith("ctl-") or did in done or not d.get("vhh_seq"):
            continue
        if len(done) >= cap:
            break
        r = R.local_rules(d["vhh_seq"])
        if r["numbering_ok"] and not r["canon_cys_ok"] and not r["invalid_aa"]:
            res = repair_canonical(d["vhh_seq"])
            done[did] = bool(res)
            if not res or res[0] in known_seqs:
                continue
            new, edits = res
            va = L.analyse(new)
            vid = dn.safe_id(f"{did[:42]}-cr")
            known_seqs.add(new)
            out.append(dict(design_id=vid, parent_id=did, generator="fwdiv", framework=d.get("framework", ""), seed=d.get("seed", ""),
                            hotspot_set=d.get("hotspot_set", ""), vhh_seq=new, cdr1=va["cdr1"], cdr2=va["cdr2"], cdr3=va["cdr3"],
                            cdr_scheme="imgt", design_cif=d.get("design_cif", ""), job_id="",
                            params_json=json.dumps(dict(repair="canonical_Cys23_Cys104_Trp118", edits=",".join(edits))),
                            created_utc=dn.utcnow()))
            log(f"repair {did}: {','.join(edits)} -> {vid}")
    return out


def watch(a):
    import time
    log = lambda m: print("[fwdiv] " + str(m), flush=True)
    div = Diversifier(mpnn_json=Path(a.mpnn), min_freq=a.min_freq)
    targets = [float(x) for x in a.targets.split(",")]
    log(f"fwdiv watch start pid={__import__('os').getpid()} interval={a.interval}s until={a.until} cap={a.cap} k={a.k}")
    out_p, rep_p = Path(a.out), Path(a.report)
    while True:
        if L.stop_requested():
            log("STOP file found, exiting"); break
        if L.past_deadline(a.until):
            log(f"deadline {a.until} reached, exiting"); break
        try:
            state = json.loads(STATE.read_text()) if STATE.exists() else {"parents": {}, "n_variants": 0}
            rows = dn.read_csv(out_p) if out_p.exists() else []
            rep = dn.read_csv(rep_p) if rep_p.exists() else []
            known_seqs = {d["vhh_seq"] for d in dn.load_designs()}
            sel_all = select_parents(a.top_n, a.per_gen, unscored_n=(a.fallback_n if L.past_deadline(a.fallback_after) else 0))
            added = 0
            rp = repair_pass(state, rows, known_seqs, log=log)
            if rp:
                rows += rp
                dn.write_csv(out_p, rows, dn.DESIGN_COLS)
                STATE.write_text(json.dumps(state, indent=1))
                added += len(rp)
            # phase 1: every new parent gets its 'lite' variants first (cheap: a handful of substitutions vs trastuzumab VH);
            # phase 2: 'deep' variants (descent to <= target identity vs the nearest known antibody)
            for tier, kk, tg, mx in (("lite", a.k_lite, [a.lite_target], a.lite_max), ("deep", a.k_deep, targets, a.max_mut)):
                for d, sc, why in sel_all:
                    pid = d["design_id"]
                    if state["parents"].get(pid, {}).get(tier):
                        continue
                    if state["n_variants"] >= a.cap:
                        log(f"variant cap {a.cap} reached"); break
                    if L.stop_requested() or L.past_deadline(a.until):
                        break
                    nr, nrep = make_rows(div, [d], kk, tg, a.min_mut, mx, tier=tier)
                    nr = [r for r in nr if r["vhh_seq"] not in known_seqs]
                    keep = {r["design_id"] for r in nr}
                    for r in nr:
                        known_seqs.add(r["vhh_seq"])
                    for r in nrep:
                        if r["design_id"] in keep:
                            r["parent_score"] = round(sc, 4); r["why"] = why
                    rows += nr; rep += [r for r in nrep if r["design_id"] in keep]
                    st = state["parents"].setdefault(pid, dict(score=round(sc, 4), why=why))
                    st[tier] = len(nr) or -1  # -1: nothing to do/emit (e.g. parent already <= lite target)
                    state["n_variants"] += len(nr)
                    added += len(nr)
                    dn.write_csv(out_p, rows, dn.DESIGN_COLS)  # atomic after every parent: S3/S4 pick variants up immediately
                    dn.write_csv(rep_p, rep)
                    STATE.write_text(json.dumps(state, indent=1))
                    log(f"parent {pid} ({why}, score {sc:.3f}) {tier}: +{len(nr)} variants; total {state['n_variants']}")
            if not added:
                log(f"{dn.utcnow()} no new parents ({len(state['parents'])} done, {state['n_variants']} variants)")
        except Exception as e:
            import traceback
            log(f"{dn.utcnow()} cycle error: {type(e).__name__}: {str(e)[:300]} | {traceback.format_exc()[-300:]}")
        for _ in range(a.interval):
            if L.stop_requested() or L.past_deadline(a.until):
                break
            time.sleep(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", nargs="*", default=None)
    ap.add_argument("--watch", action="store_true", help="loop: diversify new top parents from results/denovo/score_of3.csv")
    ap.add_argument("--interval", type=int, default=240)
    ap.add_argument("--until", default="15:30", help="local (CET) HH:MM at which the watcher exits")
    ap.add_argument("--fallback-after", default="14:40", help="local HH:MM after which unscored designs may be diversified")
    ap.add_argument("--fallback-n", type=int, default=12, help="...until this many scored parents exist (never replaces scored parents)")
    ap.add_argument("--top-n", type=int, default=40)
    ap.add_argument("--per-gen", type=int, default=3)
    ap.add_argument("--cap", type=int, default=250, help="max total variants emitted by the watcher")
    ap.add_argument("--controls-test", default=None, help="comma-separated control design_ids to diversify as a scratch test (writes only to s2_rules/cache)")
    ap.add_argument("--standin", action="store_true", help="use stand-in designs; write only to s2_rules/cache")
    ap.add_argument("--ids-file", default=None, help="restrict to design_ids listed (one per line)")
    ap.add_argument("--k", type=int, default=6, help="(non-watch test modes) deep variants per parent")
    ap.add_argument("--k-lite", type=int, default=2)
    ap.add_argument("--k-deep", type=int, default=3)
    ap.add_argument("--lite-target", type=float, default=0.66, help="identity to trastuzumab VH the lite tier must reach")
    ap.add_argument("--lite-max", type=int, default=8, help="max substitutions in the lite tier")
    ap.add_argument("--targets", default="0.68,0.66,0.64", help="identity_nearest_known targets, cycled over the K variants")
    ap.add_argument("--min-mut", type=int, default=10)
    ap.add_argument("--max-mut", type=int, default=28)
    ap.add_argument("--max-variants", type=int, default=300)
    ap.add_argument("--min-freq", type=float, default=0.015)
    ap.add_argument("--mpnn", default=str(Path(__file__).resolve().parents[2] / "results/denovo_s2/mpnn-3eak/eeb61a5b-bf0f-51e5-b55d-f5fb78aae98c.json"))
    ap.add_argument("--out", default=str(dn.OUT / "designs_fwdiv.csv"))
    ap.add_argument("--report", default=str(HERE / "fwdiv_report.csv"))
    a = ap.parse_args()
    if a.watch:
        return watch(a)
    if a.controls_test:
        want = set(a.controls_test.split(","))
        designs = [d for d in dn.load_designs() if d["design_id"] in want]
        a.out = str(CACHE / "designs_fwdiv_controls_test.csv")
        a.report = str(CACHE / "fwdiv_report_controls_test.csv")
    elif a.standin:
        designs = dn.read_csv(HERE / "standins_all3.csv")
        a.out = str(CACHE / "designs_fwdiv_standin.csv")
        a.report = str(CACHE / "fwdiv_report_standin.csv")
    else:
        designs = [d for f in (a.inp or []) for d in dn.read_csv(f)] if a.inp else dn.load_designs(include_controls=False)
        designs = [d for d in designs if d.get("generator") not in ("fwdiv", "control", "standin") and d.get("vhh_seq")]
        if a.ids_file:
            keep = set(Path(a.ids_file).read_text().split())
            designs = [d for d in designs if d["design_id"] in keep]
        # prefer designs that pass the identity-independent hard rules
        import rules as R
        rr = {r["design_id"]: r for r in R.run(designs, use_vm=False)}
        designs = [d for d in designs if rr[d["design_id"]]["numbering_ok"] and not rr[d["design_id"]]["unpaired_cys"]
                   and rr[d["design_id"]]["canon_cys_ok"] and rr[d["design_id"]].get("cdrh3_novel_local", 0)]
        seen_cdr3, uniq = set(), []
        for d in designs:  # one design per CDRH3: variants of the same CDRH3 add little
            c = L.analyse(d["vhh_seq"])["cdr3"]
            if c not in seen_cdr3:
                seen_cdr3.add(c); uniq.append(d)
        designs = uniq[: max(1, a.max_variants // a.k)]
    div = Diversifier(mpnn_json=Path(a.mpnn), min_freq=a.min_freq)
    rows, rep = make_rows(div, designs, a.k_lite, [a.lite_target], a.min_mut, a.lite_max, tier="lite")
    r2, p2 = make_rows(div, designs, a.k_deep, [float(x) for x in a.targets.split(',')], a.min_mut, a.max_mut, tier="deep")
    rows += r2; rep += p2
    dn.write_csv(a.out, rows, dn.DESIGN_COLS)
    dn.write_csv(a.report, rep)
    ok = sum(r["ok"] for r in rep)
    print(f"{len(designs)} parents -> {len(rows)} variants ({ok} below 70% to nearest known) -> {a.out}")
    for r in rep[:12]:
        print(f'  {r["design_id"]:<34} mut={r["n_mut"]:>2} id_parent={r["identity_to_parent"]} nearest {r["nearest_before"]}({r["nearest_before_id"]}) -> {r["nearest_after"]}({r["nearest_after_id"]}) trast {r["trast_before"]}->{r["trast_after"]}')


if __name__ == "__main__":
    main()
