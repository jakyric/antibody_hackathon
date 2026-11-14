"""Parse the Nebius protenix-v2 result of `ptx-top` into results/denovo/score_ptx.csv (+ report.md).

    python3 ~/.claude/skills/nebius-scientific-ai/python3 denovo/s3d_protenix/finish.py            # fetches every finished run into results/ptx-runs/<label>
    baseline/.venv/bin/python denovo/s3d_protenix/parse_ptx.py [--src results/ptx-runs]

Per design: median over the samples of every metric, plus the best sample (highest ranking_score) and the spread.
Geometry (epitope overlap, trastuzumab competition proxy) comes from S3's `denovo/s3_score/geom.py`, imported read-only by
file path (module name `s3geom`, so it cannot collide with baseline/common.py or denovo/common.py).
Chain convention of the input: VHH complex = chain A VHH, chain B HER2; Fv control = A VH, B VL, C HER2.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dn_common as dn  # noqa: E402

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("s3geom", dn.DENOVO / "s3_score" / "geom.py")
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)

ap = argparse.ArgumentParser()
ap.add_argument("--src", default=str(dn.ROOT / "results" / "ptx-runs"))
a = ap.parse_args()
SRC = Path(a.src)

man = pd.read_csv(HERE / "manifest.csv").set_index("design_id")
designs = {r["design_id"]: r for r in dn.load_designs()}
ids = sorted(man.index, key=len, reverse=True)  # longest first so ctl-fwrand-01 is not matched as a prefix of something else


labmap = dict(zip(man["label"], man.index)) if "label" in man else {}


def owner(path: Path) -> str | None:
    rel = path.relative_to(SRC)
    if rel.parts and rel.parts[0] in labmap:  # results/ptx-runs/<label>/...: one Nebius run per complex
        return labmap[rel.parts[0]]
    s = str(rel)
    for i in ids:
        if i in s:
            return i
    return None


def sample_of(path: Path) -> int:
    m = re.search(r"sample[_-]?(\d+)", path.name) or re.search(r"sample[_-]?(\d+)", str(path))
    return int(m.group(1)) if m else 0


def scalars(d: dict, prefix="") -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[prefix + k] = float(v)
        elif isinstance(v, bool):
            out[prefix + k] = float(v)
    return out


cifs, confs, paes = {}, {}, {}
import hashlib  # noqa: E402

# Nebius stores every artifact under a UUID name; result.json lists the original path + sha256 of each file, so match by hash.
for rj in sorted(SRC.glob("*/result.json")):
    run = json.loads(rj.read_text())
    if run.get("status") != "succeeded":
        continue
    by_sha = {hashlib.sha256(f.read_bytes()).hexdigest(): f for f in rj.parent.iterdir() if f.is_file()}
    for ent in run.get("native_result", {}).get("files", []):
        meta = ent.get("file", {})
        loc = by_sha.get(meta.get("sha256"))
        if loc is None:
            continue
        orig = Path(meta["path"])
        parts = orig.parts
        did = parts[1] if len(parts) > 2 and parts[0] == "outputs" and parts[1] in man.index else None
        if did is None:
            continue
        k = (did, sample_of(orig))
        sem = ent.get("semantic_type", "")
        if sem == "protein-structure-mmcif/v1":
            cifs[k] = loc
        elif sem == "structure-upstream-summary-json/v1":
            confs[k] = json.loads(loc.read_text())
print(f"found {len(cifs)} structures, {len(confs)} confidence files, {len(paes)} pae files for {len({k[0] for k in cifs})} jobs")

rows = []
for did, grp in man.groupby(level=0):
    samples = sorted(k for k in cifs if k[0] == did)
    if not samples:
        continue
    vhh_seq = designs.get(did, {}).get("vhh_seq", "")
    cd = g.find_cdrs(vhh_seq, *(designs.get(did, {}).get(c, "") for c in ("cdr1", "cdr2", "cdr3"))) if man.loc[did, "kind"] != "fv" else {}
    per = []
    for k in samples:
        m = {"sample": k[1]}
        c = confs.get(k, {})
        m.update(scalars(c))
        # chain-pair ipTM VHH-HER2: matrix [[.,x],[x,.]] -> take the off-diagonal between binder chain A and HER2 (last chain)
        cp = c.get("chain_pair_iptm")
        if isinstance(cp, list) and len(cp) >= 2:
            m["chain_pair_iptm"] = float(cp[0][len(cp) - 1])
        for key in ("chain_iptm", "chain_ptm", "chain_plddt"):
            if isinstance(c.get(key), list) and c[key]:
                m[key + "_binder"] = float(c[key][0])
                m[key + "_her2"] = float(c[key][-1])
        try:
            at = g.load_cif(cifs[k])
            fv = man.loc[did, "kind"] == "fv"
            m.update(g.complex_metrics(at, cd, vhh_chain="A", her_chain="C" if fv else "B", extra_ab_chains=("B",) if fv else ()))
        except Exception as e:  # noqa: BLE001
            m["geom_error"] = str(e)[:80]
        pj = paes.get(k)
        if pj is not None:
            pae = np.asarray(pj.get("pae") or pj.get("token_pair_pae") or pj.get("predicted_aligned_error"), float)
            if pae.ndim == 2 and pae.shape[0] == pae.shape[1]:
                nb = len(vhh_seq) if man.loc[did, "kind"] != "fv" else pae.shape[0] - g.HER2_LEN
                blk = np.concatenate([pae[:nb, nb:].ravel(), pae[nb:, :nb].ravel()])
                m["min_ipae"], m["mean_ipae"] = float(blk.min()), float(blk.mean())
        per.append(m)
    df = pd.DataFrame(per)
    row = {"design_id": did, "is_control": int(man.loc[did, "is_control"]), "ptx_n_samples": len(df)}
    for c in df.columns:
        if c in ("sample", "geom_error", "comp_error") or not pd.api.types.is_numeric_dtype(df[c]):
            continue
        row["ptx_" + c] = float(df[c].median())
    if "ranking_score" in df:
        best = df.loc[df["ranking_score"].idxmax()]
        row["ptx_iptm_best"] = float(best.get("iptm", np.nan))
        row["ptx_best_sample"] = int(best["sample"])
    if "iptm" in df:
        row["ptx_iptm_sd"] = float(df["iptm"].std(ddof=0))
    row["stage_ts"] = dn.utcnow()
    rows.append(row)

T = pd.DataFrame(rows)
if T.empty:
    sys.exit("nothing parsed: check --src and file names")
pop = T[T.is_control == 0]


def pct(s, ref, hb=True):
    r = np.sort(ref.dropna().to_numpy(float))
    if len(r) < 2:
        return pd.Series(np.nan, index=s.index)
    p = np.searchsorted(r, s.to_numpy(float), side="right") / len(r)
    return pd.Series(p if hb else 1 - p, index=s.index)


parts = [pct(T["ptx_iptm"], pop["ptx_iptm"])] if "ptx_iptm" in T else []
if "ptx_min_ipae" in T:
    parts.append(pct(T["ptx_min_ipae"], pop["ptx_min_ipae"], hb=False))
T["ptx_score"] = pd.concat(parts, axis=1).mean(axis=1) if parts else np.nan
dn.write_csv(dn.OUT / "score_ptx.csv", T.to_dict("records"), list(T.columns))
print(f"score_ptx.csv: {len(T)} rows ({len(pop)} designs, {len(T) - len(pop)} controls)")

# ---- report: controls, designs vs decoys, correlation with the other models
L = ["# Protenix-v2 (Nebius) report", "", f"{len(T)} jobs parsed, sample median per design.", ""]
cols = [c for c in ("ptx_iptm", "ptx_ptm", "ptx_ranking_score", "ptx_chain_pair_iptm", "ptx_min_ipae", "ptx_epi_overlap",
                    "ptx_cdr_fv_frac4", "ptx_clash25") if c in T]
ctl = T[T.is_control == 1][["design_id"] + cols].round(3)
L += ["## Controls", "", ctl.to_string(index=False), ""]
dec = T[T.design_id.str.startswith("ctl-fwrand")]
L += ["## Designs vs decoys (median / max)", "",
      pd.DataFrame({"designs_median": pop[cols].median(), "designs_max": pop[cols].max(),
                    "decoys_median": dec[cols].median(), "decoys_max": dec[cols].max()}).round(3).to_string(), ""]
if "ptx_iptm" in T:
    tz = T[T.design_id == "ctl-tz-fv"]["ptx_iptm"]
    L.append(f"Designs with ptx_iptm above the best decoy: {int((pop.ptx_iptm > dec.ptx_iptm.max()).sum())} of {len(pop)}"
             + (f"; trastuzumab Fv control ptx_iptm {float(tz.iloc[0]):.3f}." if len(tz) else "."))
    L.append("")
try:
    from scipy.stats import spearmanr
    S3 = pd.read_csv(dn.OUT / "score_of3.csv")
    SO = pd.read_csv(dn.OUT / "score_orth.csv")
    M = T.merge(S3, on="design_id", how="left", suffixes=("", "_of3")).merge(SO, on="design_id", how="left", suffixes=("", "_orth"))
    M = M[M.is_control == 0]
    other = [c for c in M.columns if re.search(r"(of3_iptm$|of3_consensus|dde.*iptm$|esm2.*iptm$|consensus_rank_score$|of3_min_ipae$)", c)]
    out = []
    for c in other:
        ok = M[["ptx_iptm", c]].dropna()
        if len(ok) >= 8 and ok[c].nunique() > 1:
            out.append((c, len(ok), spearmanr(ok.ptx_iptm, ok[c]).statistic))
    L += ["## Spearman of ptx_iptm vs other models (same designs)", "",
          pd.DataFrame(out, columns=["column", "n", "spearman"]).round(3).to_string(index=False) if out else "(too few overlapping rows)", ""]
except Exception as e:  # noqa: BLE001
    L.append(f"(correlation skipped: {e})")
(HERE / "report.md").write_text("\n".join(L) + "\n")
print("\n".join(L[:40]))
