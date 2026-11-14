#!/usr/bin/env python3
"""Local layout check of the deck: renders the generated slides with headless Chromium, four per image.

    denovo/.venv/bin/python denovo/slides/preview.py        # writes denovo/slides/preview/sheet-N.png

Approximation of the Slides renderer (plain browser rules; x-shape arrows are not drawn). Used only to catch wraps and overflow.
"""
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = HERE / "artifact" / "project"
OUT = HERE / "preview"
OUT.mkdir(exist_ok=True)
idx = json.loads((PROJ / "deck.json").read_text())
css = (PROJ / "deck.css").read_text()
fonts = "".join(f'<link rel="stylesheet" href="{f["href"]}">' for f in idx["faces"].values() if "href" in f)
base = """
* { margin: 0; box-sizing: border-box }
body { font: 32px/1.2 sans-serif; background: #888 }
.wrap { width: 960px; height: 540px; overflow: hidden; position: relative }
.wrap > section { width: 1920px !important; height: 1080px !important; overflow: hidden !important; position: absolute; left: 0; top: 0;
  transform: scale(0.5); transform-origin: top left; isolation: isolate }
aside { display: none !important }
h1 { font-size: 96px; font-weight: 600 } h2 { font-size: 64px; font-weight: 600 } h3 { font-size: 44px; font-weight: 600 }
p, li, h1, h2, h3, td, th { overflow-wrap: break-word }
table { width: 100%; border-collapse: separate; border-spacing: 0; border: solid #0003; border-width: 2px 0 0 2px }
th, td { padding: .35em .6em; text-align: left; border: solid #0003; border-width: 0 2px 2px 0 }
x-shape { display: block; border-radius: 4px } x-shape[kind=arrow-right] { clip-path: polygon(0 25%, 60% 25%, 60% 0, 100% 50%, 60% 100%, 60% 75%, 0 75%) }
.grid { display: grid; grid-template-columns: 960px 960px; width: 1920px }
"""
ids = idx["order"]
for n in range(0, len(ids), 4):
    chunk = ids[n:n + 4]
    body = "".join(f'<div class="wrap">{(PROJ / "slides" / f"{i}.html").read_text()}</div>' for i in chunk)
    page = f"<!doctype html><meta charset=utf-8>{fonts}<style>{base}{css}</style><div class=grid>{body}</div>"
    h = OUT / f"sheet-{n // 4 + 1}.html"
    h.write_text(page)
    png = OUT / f"sheet-{n // 4 + 1}.png"
    r = subprocess.run(["chromium", "--headless=new", "--no-sandbox", "--disable-gpu", "--hide-scrollbars", "--virtual-time-budget=6000",
                        "--window-size=1920,1080", f"--screenshot={png}", f"file://{h}"], capture_output=True, text=True, timeout=120)
    print(png.name, "ok" if png.exists() else f"FAILED {r.stderr[-300:]}")
