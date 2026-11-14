"""Build S3's controls and calibration sets (deterministic, seed 7).

Writes
  results/denovo/designs_control.csv   VHH-type controls (generator=control, ids ctl-*) for ALL stages
  denovo/s3_score/calib_set.csv        calibration complexes private to S3 (Fv binders vs non-binders, other-antigen
                                       VHH pairs). Deliberately NOT named designs_*.csv so other stages never read it.

Run: denovo/.venv/bin/python denovo/s3_score/make_calib.py
"""
import json
import random
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "denovo"))
import dn_common as dn  # noqa: E402

HERE = Path(__file__).resolve().parent
RNG = random.Random(7)
AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA19 = AA20.replace("C", "")
WT_H3 = "WGGDGFYAMD"
VH_TZ = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
         "SRWGGDGFYAMDYWGQGTLVTVSS")
VL_TZ = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
# NbBCII10-FGLA (3EAK), the framework chosen by S1 (denovo/s1_design/frameworks/bcii10.json)
BC_PRE = "QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
BC_POST = "WGQGTLVTVSS"
BC_CDR3_WT = "AAVRGYFMRLPSSHNFRY"
BC_CDR1, BC_CDR2 = "GGSEYSYSTFS", "IASMGGLT"


def ham(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 99


def row(did, seq, cdr3, gen="control", fw="bcii10", note="", **params):
    return {"design_id": did, "parent_id": "", "generator": gen, "framework": fw, "seed": 0, "hotspot_set": "",
            "vhh_seq": seq, "cdr1": BC_CDR1 if fw == "bcii10" else "", "cdr2": BC_CDR2 if fw == "bcii10" else "",
            "cdr3": cdr3, "cdr_scheme": "imgt" if fw == "bcii10" else "fv", "design_cif": "", "job_id": "",
            "params_json": json.dumps({"note": note, **params}), "created_utc": dn.utcnow()}


def controls():
    rows = [row("ctl-bcii10-wt", BC_PRE + BC_CDR3_WT + BC_POST, BC_CDR3_WT,
                note="NbBCII10-FGLA wild type: framework-only VHH, binds beta-lactamase not HER2 (negative)"),
            row("ctl-tz-fv", VH_TZ, "SRWGGDGFYAMDY", fw="trastuzumab", kind="fv", vl_seq=VL_TZ,
                note="trastuzumab Fv (VH+VL) as the full-signal positive; scored as a 3-chain complex, not a VHH")]
    for i in range(1, 13):
        n = RNG.randint(8, 14)
        h3 = "AA" + "".join(RNG.choice(AA19) for _ in range(n)) + "Y"
        rows.append(row(f"ctl-fwrand-{i:02d}", BC_PRE + h3 + BC_POST, h3,
                        note="BCII10-FGLA framework + random CDR3: no HER2-directed design (decoy / harness stand-in)"))
    dn.write_csv(dn.OUT / "designs_control.csv", rows, dn.DESIGN_COLS)
    return rows


def fv_calibration():
    spr = pd.read_csv(ROOT / "data" / "her2" / "trastuzumab_her2_spr_all_768.csv")
    spr = spr[spr.cdrh3.str.len().eq(10)].copy()
    spr["ed"] = [ham(s, WT_H3) for s in spr.cdrh3]
    spr = spr[(spr.ed >= 3) & (spr.ed <= 6)].drop_duplicates("cdrh3")
    bind = spr[(spr.binding_result == "quantified_binding") & (spr.kd_nM <= 10)]
    non = spr[spr.binding_result == "no_reportable_kd"]
    out = []
    for lab, df, tag in ((1, bind, "spr-bind"), (0, non, "spr-non")):
        df = df.sample(n=min(12, len(df)), random_state=7)
        for k, r in enumerate(df.itertuples(), 1):
            out.append({"calib_id": f"cal-{tag}-{k:02d}", "label": lab, "group": tag, "kind": "fv", "vhh_seq": r.vh_sequence,
                        "vl_seq": r.vl_sequence, "antigen_seq": "", "antigen_name": "her2",
                        "note": f"CDRH3 {r.cdrh3}, {r.ed} changes, KD {r.kd_nM if lab else 'none reported'} nM"})
    with zipfile.ZipFile(ROOT / "data" / "her2" / "tz_her2_facs_3classes_524346.csv.zip") as z:
        name = [n for n in z.namelist() if n.endswith(".csv")][0]
        facs = pd.read_csv(z.open(name))
    facs = facs[facs.edit_distance.between(4, 6)]
    for lab, cls, tag in ((1, "high", "facs-high"), (0, "low", "facs-low")):
        df = facs[facs.binding_class == cls].sample(n=16, random_state=7)
        for k, r in enumerate(df.itertuples(), 1):
            out.append({"calib_id": f"cal-{tag}-{k:02d}", "label": lab, "group": tag, "kind": "fv", "vhh_seq": r.vh_sequence,
                        "vl_seq": r.vl_sequence, "antigen_seq": "", "antigen_name": "her2",
                        "note": f"CDRH3 {r.cdrh3}, {r.edit_distance} changes, FACS {cls}"})
    return out


if __name__ == "__main__":
    c = controls()
    cal = fv_calibration()
    # keep any VHH-antigen pairs added by make_vhh_calib.py
    keep = []
    f = HERE / "calib_set.csv"
    if f.exists():
        keep = [r for r in dn.read_csv(f) if r["kind"] == "vhh"]
    dn.write_csv(f, cal + keep, ["calib_id", "label", "group", "kind", "vhh_seq", "vl_seq", "antigen_seq", "antigen_name", "note"])
    print(f"controls: {len(c)} -> results/denovo/designs_control.csv ; calibration complexes: {len(cal)} Fv + {len(keep)} VHH")
