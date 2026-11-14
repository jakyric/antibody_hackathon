"""Stage 5: GetContacts on every predicted complex, two views, compared with the crystal (1N8Z).

  interface : Fv (H, L) vs HER2 (G)                    -> are the crystal Fv-HER2 contacts kept?
  intra     : CDR-H3 (H 99-108) vs rest of the Fv       -> does the loop still sit in the framework / Vernier zone?
The second view covers the R98 (H94) - D108 salt bridge, the Vernier contacts of positions 8-10 and the VL contacts.
Contacts are compared at residue-pair level. HER2 residues 581-590 are unresolved in the crystal and are ignored when
counting 'new' contacts.
"""
from __future__ import annotations

import os
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from .common import out, save_summary

POLAR = {"hbbb", "hbsb", "hbss", "hbls", "hblb", "sb", "pc", "ps", "ts", "hb"}
SEL_INTERFACE = ("chain H or chain L", "chain G")
SEL_INTRA = ("chain H and resid 99 to 108", "(chain H and not resid 99 to 108) or chain L")


def _gc(cfg: dict, pdb: Path, tsv: Path, sel: tuple[str, str]) -> None:
    if tsv.exists():
        return
    c = cfg["contacts"]
    py = str(Path(c["getcontacts_python"]).expanduser())
    script = str(Path(c["getcontacts_dir"]).expanduser() / "get_static_contacts.py")
    subprocess.run([py, script, "--structure", str(pdb), "--output", str(tsv), "--sele", sel[0], "--sele2", sel[1],
                    "--itypes", "all"], capture_output=True, check=False)


def read_interface(tsv: Path) -> dict:
    pairs = defaultdict(set)
    for line in open(tsv):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        a, b = p[2].split(":"), p[3].split(":")
        if a[0] == "G":
            a, b = b, a
        if b[0] == "G" and a[0] in "HL":
            pairs[(a[0], int(a[2]), int(b[2]))].add(p[1])
    return pairs


def read_intra(tsv: Path) -> dict:
    pairs = defaultdict(set)
    for line in open(tsv):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        a, b = p[2].split(":"), p[3].split(":")
        if not (a[0] == "H" and 99 <= int(a[2]) <= 108):
            a, b = b, a
        if a[0] == "H" and 99 <= int(a[2]) <= 108:
            pairs[(int(a[2]), b[0], int(b[2]))].add(p[1])
    return pairs


def run(cfg: dict) -> None:
    work = out(cfg, "contacts")
    pdbs, tsv_if, tsv_in = work / "pdb", work / "interface", work / "intra"
    for d in (pdbs, tsv_if, tsv_in):
        d.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PREP_IN=str(out(cfg, "structures")), PREP_OUT=str(pdbs),
               PREP_CRYSTAL=str(Path(cfg["paths"]["crystal"]).expanduser().resolve()))
    subprocess.run([cfg["contacts"]["pymol"], "-cq", str(Path(__file__).with_name("pymol_prep.py"))], env=env,
                   capture_output=True, check=True)
    files = sorted(pdbs.glob("*.pdb"))
    jobs = []
    with ThreadPoolExecutor(cfg["contacts"]["jobs"]) as ex:
        for f in files:
            jobs.append(ex.submit(_gc, cfg, f, tsv_if / f"{f.stem}.tsv", SEL_INTERFACE))
            jobs.append(ex.submit(_gc, cfg, f, tsv_in / f"{f.stem}.tsv", SEL_INTRA))
        for j in jobs:
            j.result()

    ref_if, ref_in = read_interface(tsv_if / "crystal.tsv"), read_intra(tsv_in / "crystal.tsv")
    ref_h3 = {k for k in ref_if if k[0] == "H" and 99 <= k[1] <= 108}
    ref_in_polar = {k for k, v in ref_in.items() if v & POLAR}
    print(f"  crystal: {len(ref_if)} Fv-HER2 pairs ({len(ref_h3)} with CDR-H3); {len(ref_in)} CDR-H3-vs-rest pairs ({len(ref_in_polar)} polar)")
    rows = []
    for f in files:
        if f.stem == "crystal":
            continue
        a, b = read_interface(tsv_if / f"{f.stem}.tsv"), read_intra(tsv_in / f"{f.stem}.tsv")
        ka, kb = set(a), set(b)
        resolved = lambda S: {k for k in S if not 581 <= k[2] <= 590}
        sb = b.get((108, "H", 98), set()) | b.get((108, "H", 98), set())
        rows.append(dict(
            id=f.stem, if_pairs=len(ka), if_recall=len(ka & set(ref_if)) / len(ref_if), if_missing=len(set(ref_if) - ka),
            if_h3_missing=len(ref_h3 - ka), if_new=len(resolved(ka) - set(ref_if)),
            in_pairs=len(kb), in_recall=len(kb & set(ref_in)) / len(ref_in), in_missing=len(set(ref_in) - kb),
            in_polar_recall=len(ref_in_polar & {k for k in kb if b[k] & POLAR}) / len(ref_in_polar),
            salt_bridge_r98_pos10="sb" in sb))
    df = pd.DataFrame(rows)
    df.to_csv(out(cfg, "contacts.csv"), index=False)
    ctrl = df[df.id.str.startswith("ctrl_")].set_index("id")
    print(ctrl[["if_recall", "if_h3_missing", "in_recall", "in_polar_recall", "salt_bridge_r98_pos10"]].round(2).to_string())
    save_summary(cfg, crystal_if_pairs=len(ref_if), crystal_h3_pairs=len(ref_h3), crystal_intra_pairs=len(ref_in))
