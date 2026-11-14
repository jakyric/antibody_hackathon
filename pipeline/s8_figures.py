"""Stage 8 (optional): structure figures with PyMOL (overview, Vernier close-up, predicted pose vs crystal)."""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

from PIL import Image, ImageChops

from .common import out


def _trim(path: Path, margin: int = 40) -> None:
    im = Image.open(path).convert("RGB")
    box = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
    if box:
        im.crop((max(0, box[0] - margin), max(0, box[1] - margin), min(im.width, box[2] + margin), min(im.height, box[3] + margin))).save(path)


def run(cfg: dict) -> None:
    top = json.load(open(out(cfg, "final") / "top_candidates.json"))["candidates"][0]["id"]
    figs = out(cfg, "figures")
    env = dict(os.environ, FIG_CRYSTAL=str(Path(cfg["paths"]["crystal"]).expanduser().resolve()),
               FIG_PRED=str(out(cfg, "structures") / f"{top}.pdb"), FIG_OUT=str(figs))
    subprocess.run([cfg["contacts"]["pymol"], "-cq", str(Path(__file__).with_name("pymol_figures.py"))], env=env, check=True, capture_output=True)
    for f in figs.glob("fig*.png"):
        _trim(f)
    print(f"  figures for {top} in {figs}")
