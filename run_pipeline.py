#!/usr/bin/env python3
"""Run the HER2 CDR-H3 + Vernier design pipeline.

    python run_pipeline.py --config config/pipeline.yaml                 # all stages
    python run_pipeline.py --config config/pipeline.yaml --stages contacts score rank
    python run_pipeline.py --list

Stages are decoupled by files in paths.out_dir, so any stage can be re-run on its own.
"""
import argparse
import importlib
import sys
import time

from pipeline.common import load_config

STAGES = [  # name, module, one-line description
    ("enumerate", "s1_enumerate", "BLOSUM45 variants of the best tested CDR-H3, novelty / liability filters"),
    ("panel", "s2_panel", "diverse CDR-H3 shortlist + BLOSUM45 Vernier-zone genotypes"),
    ("validate", "s3_validate", "organisers' abguard gate (novelty, liabilities, Cys parity)"),
    ("fold", "s4_fold", "JapanFold OpenDDE antibody-antigen folds of all sequences + controls"),
    ("contacts", "s5_contacts", "GetContacts: Fv-HER2 interface and CDR-H3-vs-rest, vs the crystal"),
    ("pose", "s5b_pose", "RMSD to the crystal after fitting the Fv framework: same HER2 epitope?"),
    ("score", "s6_score", "ipSAE and iLIS from PAE for designs passing the contact pre-filter"),
    ("rank", "s7_rank", "final filter, ranking and export (JSON / FASTA / CSV)"),
    ("figures", "s8_figures", "PyMOL structure figures (optional)"),
    ("plots", "s10_plots", "score distributions and ipSAE-iLIS plane with trastuzumab marked (optional)"),
    ("deck", "s9_deck", "PowerPoint summary from a template (optional)"),
]
DEFAULT = [s[0] for s in STAGES[:8]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--stages", nargs="+", choices=[s[0] for s in STAGES], help="default: enumerate .. rank")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for n, _, d in STAGES:
            print(f"{n:10s} {d}")
        return
    cfg = load_config(a.config)
    for name, mod, desc in STAGES:
        if name not in (a.stages or DEFAULT):
            continue
        print(f"[{name}] {desc}")
        t = time.time()
        importlib.import_module(f"pipeline.{mod}").run(cfg)
        print(f"[{name}] done in {time.time() - t:.0f}s\n")


if __name__ == "__main__":
    sys.exit(main())
