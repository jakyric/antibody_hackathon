"""Batch TNP (Therapeutic Nanobody Profiler) without the per-sequence CLI overhead.

Run inside the micromamba env `tnp`:
    micromamba run -n tnp python denovo/s4_props/tnp_batch.py --in seqs.csv --out tnp.csv --pdb-dir pdbs [--ncores 8]

seqs.csv needs columns: design_id, vhh_seq.  Output columns: design_id, tnp_* metrics, tnp_flag_*, tnp_n_red/amber.

Why not `TNP --file --ncores`: that CLI launches NanoBodyBuilder2 as a subprocess per sequence with
`--n_threads 1` on CPU (about 4 min each).  Here ONE process loads NanoBodyBuilder2 once, models on the GPU
(about 2 s per sequence incl. OpenMM refinement), and a process pool only runs TNP's own downstream maths
(the numbering, the CDR3 compactness and the patch hydrophobicity/charge calculators), which are CPU-only.
The metric code is imported unchanged from tools/TNP (theraprofnano), so the numbers match `TNP --seq`
for the same model; the only difference is that the model is generated in-process.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import multiprocessing as mp
import os
import sys
import tempfile
import time
from contextlib import redirect_stdout
from pathlib import Path

H_SCALE = 0  # Kyte & Doolittle, TNP default


def _strip_h(pdb_path: str) -> None:
    from Bio.PDB import PDBIO, PDBParser, Select
    import re
    hyd = re.compile("[123 ]*H.*")

    class NoH(Select):
        def accept_atom(self, atom):
            return 0 if hyd.match(atom.get_id()) else 1

    s = PDBParser(QUIET=True).get_structure("nb", pdb_path)
    io_ = PDBIO()
    io_.set_structure(s)
    io_.save(pdb_path, NoH())


def model_all(rows, pdb_dir: Path, log):
    """Model every sequence with NanoBodyBuilder2 (GPU if available) -> pdb_dir/<id>.pdb (no hydrogens)."""
    import torch  # noqa: F401  (ImmuneBuilder picks the device itself)
    from ImmuneBuilder import NanoBodyBuilder2
    pdb_dir.mkdir(parents=True, exist_ok=True)
    nbb = NanoBodyBuilder2()
    failed = {}
    for r in rows:
        out = pdb_dir / f"{r['design_id']}.pdb"
        if out.exists() and out.stat().st_size > 1000:
            continue
        t = time.time()
        try:
            nb = nbb.predict({"H": r["vhh_seq"]})
            nb.save(str(out))
            _strip_h(str(out))
            log(f"  modelled {r['design_id']} in {time.time() - t:.1f}s")
        except Exception as e:  # noqa: BLE001
            failed[r["design_id"]] = f"nbb2: {type(e).__name__}: {e}"[:200]
            log(f"  NBB2 failed {r['design_id']}: {failed[r['design_id']]}")
    return failed


def _profile(args):
    """TNP downstream metrics for one (id, seq, pdb).  Mirrors theraprofnano.cli.profile_sequence."""
    did, seq, pdb = args
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "TNP" / "src"))
    res = {"design_id": did}
    try:
        with redirect_stdout(io.StringIO()):
            from theraprofnano.CDR_Profiler.CDR_Assigner import main as region_and_aa_dicts
            from theraprofnano.CDR_Profiler.CDR3_Conf_Assigner import main_compactness
            from theraprofnano.Hydrophobicity_and_Charge_Profiler.Hydrophobicity_and_Charge_Assigner import \
                CreateAnnotation
            from theraprofnano.thresholds import assign_flag
            from anarci import run_anarci
            with tempfile.TemporaryDirectory() as td:
                _, lens = region_and_aa_dicts(did + "_H", seq, "H", td, verbose=False)
            total = sum(lens["imgt"][k] for k in ("cdrh1", "cdrh2", "cdrh3"))
            l3 = lens["imgt"]["cdrh3"]
            rho = main_compactness(pdb, "imgt", verbose=False)
            comp = (l3 / rho) if rho is not None else float("nan")
            hc = CreateAnnotation(H_SCALE, 7.4, pdb, "IG", "imgt", verbose=False)[H_SCALE]
            try:
                num = run_anarci(seq, scheme="IMGT", allowed_species=None)[1][0][0][0]
                tetrad = "".join(n[1] for n in num if n[0][0] in (42, 49, 50, 52))
            except Exception:  # noqa: BLE001
                tetrad = ""
        vals = {"L": total, "L3": l3, "C": comp, "PSH": hc["Patch_Hydrophob_CDR"],
                "PPC": hc["Patch_Pos_Charge_CDR"], "PNC": hc["Patch_Neg_Charge_CDR"]}
        flags = {k: assign_flag(k, v) for k, v in vals.items()}
        res.update({"tnp_cdr_len": total, "tnp_cdr3_len": l3, "tnp_cdr3_compact": comp, "tnp_psh": vals["PSH"],
                    "tnp_ppc": vals["PPC"], "tnp_pnc": vals["PNC"], "tnp_tetrad": tetrad})
        res.update({f"tnp_flag_{k.lower()}": v for k, v in flags.items()})
        res["tnp_n_red"] = sum(v == "red" for v in flags.values())
        res["tnp_n_amber"] = sum(v == "amber" for v in flags.values())
        res["tnp_error"] = ""
    except Exception as e:  # noqa: BLE001
        res["tnp_error"] = f"{type(e).__name__}: {e}"[:200]
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--ncores", type=int, default=8)
    a = ap.parse_args()
    log = lambda m: print(m, flush=True)  # noqa: E731
    rows = list(csv.DictReader(open(a.inp)))
    pdb_dir = Path(a.pdb_dir)
    t0 = time.time()
    failed = model_all(rows, pdb_dir, log)
    log(f"modelling done in {time.time() - t0:.0f}s ({len(failed)} failed)")
    jobs = [(r["design_id"], r["vhh_seq"], str(pdb_dir / f"{r['design_id']}.pdb")) for r in rows
            if r["design_id"] not in failed]
    t1 = time.time()
    with mp.get_context("spawn").Pool(max(1, min(a.ncores, len(jobs) or 1))) as pool:
        out = pool.map(_profile, jobs, chunksize=1)
    out += [{"design_id": k, "tnp_error": v} for k, v in failed.items()]
    log(f"TNP metrics for {len(jobs)} seqs in {time.time() - t1:.0f}s")
    cols = ["design_id", "tnp_cdr_len", "tnp_cdr3_len", "tnp_cdr3_compact", "tnp_psh", "tnp_ppc", "tnp_pnc",
            "tnp_tetrad", "tnp_flag_l", "tnp_flag_l3", "tnp_flag_c", "tnp_flag_psh", "tnp_flag_ppc", "tnp_flag_pnc",
            "tnp_n_red", "tnp_n_amber", "tnp_error"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
