"""Run with PyMOL:  PREP_IN=<dir of pdb> PREP_OUT=<dir> pymol -cq pipeline/pymol_prep.py
Adds hydrogens (GetContacts needs them for hydrogen bonds) and standardises chains: H (VH), L (VL), G (HER2).
JapanFold predictions use A = VH, B = VL, C = HER2 domain IV numbered 1-119; HER2 is renumbered onto the 1N8Z
numbering (489-580 and 591-607; the 10 residues in between are not resolved in the crystal and get 581-590).
The crystal (PREP_CRYSTAL=<cif>) uses B = VH, A = VL, C+D = HER2 and keeps its own numbering."""
import glob
import os

from pymol import cmd


def her2_num(i):
    return 488 + i if i <= 92 else (581 + (i - 93) if i <= 102 else 591 + (i - 103))


def save(src, dst, crystal=False):
    cmd.delete("all")
    cmd.load(src, "m")
    cmd.remove("m and (solvent or hydro)")
    if crystal:
        cmd.alter("m and chain B", "chain='H'")
        cmd.alter("m and chain A", "chain='L'")
        cmd.alter("m and chain C+D", "chain='G'")
    else:
        cmd.alter("m and chain A", "chain='H'")
        cmd.alter("m and chain B", "chain='L'")
        cmd.alter("m and chain C", "resi=str(her2num(int(resi)))", space={"her2num": her2_num})
        cmd.alter("m and chain C", "chain='G'")
    cmd.sort()
    cmd.h_add("m")
    cmd.save(dst, "m")


out_dir = os.environ["PREP_OUT"]
os.makedirs(out_dir, exist_ok=True)
if os.environ.get("PREP_CRYSTAL"):
    save(os.environ["PREP_CRYSTAL"], f"{out_dir}/crystal.pdb", crystal=True)
for f in sorted(glob.glob(os.environ.get("PREP_IN", "") + "/*.pdb")):
    dst = f"{out_dir}/{os.path.basename(f)}"
    if not os.path.exists(dst):
        save(f, dst)
print("RES prepared", len(glob.glob(out_dir + "/*.pdb")))
