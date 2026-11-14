"""Build a PyMOL session + PNGs of a few RFD4 designs on HER2 domain IV, with trastuzumab's Fv (1N8Z) superposed as a
transparent reference for the competition geometry.

Run from the repo root with the PyMOL conda env:
    micromamba run -p $PYMOL_ENV pymol -cq denovo/validation/make_pymol.py
Reads denovo/validation/pymol/show_top10.csv (design_id, framework, vhh_seq, cdr1, cdr2, cdr3, design_cif).
Design CIF layout (checked): chain A = HER2 489-580 (resi 1-92), chain B = designed VHH (resi 1..N), chain C = HER2 591-607 (resi 1-17).
Colours: HER2 grey, VHH framework light blue, CDR1 blue, CDR2 green, CDR3 red, Cys sticks yellow, trastuzumab epitope orange sticks,
trastuzumab Fv light orange (VL) / orange (VH), transparent.
"""
import csv
import json
import os

from pymol import cmd

ROOT = os.getcwd()
OUT = f"{ROOT}/denovo/validation/pymol"
rows = list(csv.DictReader(open(f"{OUT}/show_top10.csv")))
epi = [int(r["auth_seq"]) for r in json.load(open(f"{ROOT}/denovo/targets/epitope.json"))["residues"]]


def her2(obj, n):
    return f"({obj} and chain A and resi {n - 488})" if n <= 580 else f"({obj} and chain C and resi {n - 590})"


cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("cartoon_side_chain_helper", 1)
cmd.set("ray_shadow", 0)
cmd.set("depth_cue", 0)
cmd.set("cartoon_transparency", 0.0)

# reference: trastuzumab Fv bound to HER2 (1N8Z crystal; chains A = VL, B = VH, C/D = HER2)
cmd.load(f"{ROOT}/denovo/targets/ref_1n8z.cif", "tz_ref")
cmd.hide("everything", "tz_ref")
cmd.show("cartoon", "tz_ref and chain A+B")
cmd.color("wheat", "tz_ref and chain A")
cmd.color("orange", "tz_ref and chain B")
cmd.set("cartoon_transparency", 0.6, "tz_ref")
cmd.hide("everything", "tz_ref and chain C+D")

for i, r in enumerate(rows):
    obj = f"d{i + 1}"
    cmd.load(f"{ROOT}/{r['design_cif']}", obj)
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.color("gray70", f"{obj} and chain A+C")
    cmd.color("lightblue", f"{obj} and chain B")
    seq = r["vhh_seq"]
    notes = []
    for nm, col in (("cdr1", "tv_blue"), ("cdr2", "tv_green"), ("cdr3", "tv_red")):
        s = r[nm]
        k = seq.find(s)
        if k < 0:
            notes.append(f"{nm} not found in vhh_seq")
            continue
        cmd.color(col, f"{obj} and chain B and resi {k + 1}-{k + len(s)}")
        cmd.show("sticks", f"{obj} and chain B and resi {k + 1}-{k + len(s)} and not name H*")
        cmd.set("stick_radius", 0.12, obj)
    cmd.show("sticks", f"{obj} and chain B and resn CYS")
    cmd.color("yellow", f"{obj} and chain B and resn CYS and elem S")
    cmd.select(f"{obj}_epi", " or ".join(her2(obj, n) for n in epi))
    cmd.show("sticks", f"{obj}_epi and not name H*")
    cmd.color("orange", f"{obj}_epi and elem C")
    cmd.set("stick_radius", 0.14, f"{obj}_epi")
    if notes:
        print(r["design_id"], notes)

# superpose the 1N8Z HER2 onto design 1's HER2 once (all designs share the HER2 frame because HER2 is held rigid)
cmd.align("tz_ref and chain C+D and name CA", "d1 and chain A+C and name CA")
print("tz_ref aligned on d1 HER2")


def snap(tag, obj, views=2):
    cmd.disable("all")
    cmd.enable("tz_ref")
    cmd.enable(obj)
    cmd.enable(f"{obj}_epi")
    cmd.zoom(f"{obj} and chain B or {obj}_epi", 4)
    cmd.orient(f"{obj} and chain B or {obj}_epi")
    for v in range(1):
        if v == 1:
            cmd.turn("y", 90)
        cmd.ray(1300, 950)
        cmd.png(f"{OUT}/{tag}_v{v + 1}.png", dpi=150)


for i, r in enumerate(rows):
    snap(r["design_id"], f"d{i + 1}")
    cmd.scene(f"scene_{r['design_id']}", "store")

# VHH-only view with CDR colours for fold plausibility (design 1)
cmd.disable("all")
cmd.enable("d1")
cmd.hide("cartoon", "d1 and chain A+C")
cmd.hide("sticks", "d1_epi")
cmd.orient("d1 and chain B")
cmd.ray(1100, 850)
cmd.png(f"{OUT}/{rows[0]['design_id']}_vhh_only.png", dpi=150)
cmd.show("cartoon", "d1 and chain A+C")
cmd.show("sticks", "d1_epi")

cmd.disable("all")
cmd.enable("tz_ref")
cmd.enable("d1")
cmd.save(f"{OUT}/rfd4_top10.pse")
print("saved", f"{OUT}/rfd4_top10.pse")
