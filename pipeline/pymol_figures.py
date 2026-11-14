"""Run with PyMOL (see s8_figures): renders the three structure figures.
Env: FIG_CRYSTAL (1N8Z cif), FIG_PRED (predicted pdb of the top design: A=VH, B=VL, C=HER2), FIG_OUT (dir)."""
import json
import os

import numpy as np
from pymol import cmd

crystal, pred, outd = os.environ["FIG_CRYSTAL"], os.environ["FIG_PRED"], os.environ["FIG_OUT"]
os.makedirs(outd, exist_ok=True)
cmd.set("ray_shadows", 0); cmd.set("antialias", 2); cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)
cmd.set("label_font_id", 7); cmd.set("label_size", 20); cmd.set("label_color", "black"); cmd.set("label_outline_color", "white")
cmd.set("label_position", [0.0, 0.0, 6.0]); cmd.set("stick_radius", 0.22); cmd.set("dash_radius", 0.05); cmd.set("dash_gap", 0.4)
cmd.set("depth_cue", 0)
for n, c in (("c_vh", [0.80, 0.82, 0.85]), ("c_vl", [0.65, 0.69, 0.75]), ("c_her2", [0.62, 0.85, 0.88]), ("c_var", [0.95, 0.55, 0.10]),
             ("c_cont", [0.98, 0.80, 0.10]), ("c_vern", [0.80, 0.12, 0.55]), ("c_pred", [0.95, 0.55, 0.10])):
    cmd.set_color(n, c)
H3_VAR, H3_CONT = "100+101+102+106+107+108", "99+103+104+105"           # CDR-H3 positions 2,3,4,8,9,10 / 1,5,6,7 (B99-B108)
VERNIER = [("B", 98, "H94"), ("B", 97, "H93"), ("B", 47, "H47"), ("A", 36, "L36"), ("A", 46, "L46"), ("A", 49, "L49")]


def style():
    cmd.hide("everything", "cx"); cmd.show("cartoon", "cx")
    cmd.color("c_vh", "cx and chain B"); cmd.color("c_vl", "cx and chain A"); cmd.color("c_her2", "cx and chain C+D")
    cmd.color("lightorange", "cx and chain B and resi 99-108")
    for sel, col in ((H3_VAR, "c_var"), (H3_CONT, "c_cont")):
        cmd.show("sticks", f"cx and chain B and resi {sel}"); cmd.color(col, f"cx and chain B and resi {sel} and elem C")
    for ch, r, _ in VERNIER:
        cmd.show("sticks", f"cx and chain {ch} and resi {r}"); cmd.color("c_vern", f"cx and chain {ch} and resi {r} and elem C")
    cmd.util.cnc("cx")


cmd.load(crystal, "cx"); cmd.remove("cx and (hydro or solvent)"); style()
cmd.orient("cx and (chain A+B)"); cmd.zoom("cx and (chain A+B)", 2); cmd.turn("y", -10)
cmd.png(f"{outd}/fig1_overview.png", 1600, 1100, dpi=200, ray=1)
cmd.hide("cartoon", "cx and chain C+D")
for p in (2, 3, 4, 8, 9, 10):
    cmd.label(f"cx and chain B and resi {98 + p} and name CA", f"'H3-{p}'")
for ch, r, name in VERNIER:
    cmd.label(f"cx and chain {ch} and resi {r} and name CA", f"'{name}'")
side = "not name N+CA+C+O"
cmd.distance("vh3", "(cx and chain B and resi 97+98+47 and " + side + ") or (cx and chain A and resi 36+46+49 and " + side + ")",
             f"cx and chain B and resi {H3_VAR} and " + side, 4.0, mode=0)
cmd.color("c_vern", "vh3"); cmd.hide("labels", "vh3")
cmd.orient("cx and chain B and resi 97-108 or (cx and chain A and resi 36+46+49)")
cmd.png(f"{outd}/fig2_closeup.png", 1600, 1100, dpi=200, ray=1)

cmd.delete("all"); cmd.load(crystal, "cx"); cmd.remove("cx and (hydro or solvent)"); cmd.load(pred, "p")
fit = cmd.align("p and chain A and not resi 99-108 and name CA", "cx and chain B and not resi 99-108 and name CA")


def xyz(sel):
    return np.array([a.coord for a in cmd.get_model(sel).atom])


rm = lambda a, b: float(np.sqrt(((xyz(a) - xyz(b)) ** 2).sum(1).mean()))
stats = {"framework_fit_rmsd": round(fit[0], 2),
         "cdrh3_ca_rmsd": round(rm("p and chain A and resi 99-108 and name CA", "cx and chain B and resi 99-108 and name CA"), 2),
         "vl_ca_rmsd": round(rm("p and chain B and name CA", "cx and chain A and name CA"), 2)}
json.dump(stats, open(f"{outd}/pose_rmsd.json", "w"))
cmd.hide("everything", "all"); cmd.show("cartoon", "all")
cmd.color("c_vh", "cx and chain B"); cmd.color("c_vl", "cx and chain A"); cmd.color("c_her2", "cx and chain C+D"); cmd.set("cartoon_transparency", 0.55, "cx")
cmd.color("grey50", "p and chain A"); cmd.color("grey70", "p and chain B"); cmd.color("c_pred", "p and chain C"); cmd.color("c_vern", "p and chain A and resi 99-108")
cmd.show("sticks", "p and chain A and resi 99-108"); cmd.color("c_pred", "p and chain A and resi 99-108 and elem C")
cmd.show("sticks", "cx and chain B and resi 99-108"); cmd.color("tv_green", "cx and chain B and resi 99-108 and elem C")
cmd.set("cartoon_transparency", 0.0, "p"); cmd.util.cnc("all")
cmd.zoom("cx and chain B and resi 96-110 or (cx and chain C+D and resi 489-580)", -4)
cmd.png(f"{outd}/fig3_pose.png", 1600, 1100, dpi=200, ray=1)
print("RES", json.dumps(stats))
