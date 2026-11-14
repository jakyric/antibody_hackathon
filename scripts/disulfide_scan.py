"""Look for sites for a new VH-VL disulfide bridge in trastuzumab.

Run:  python scripts/disulfide_scan.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from abx.thermo import (HEAVY_KEEP, KNOWN_ENGINEERED_PAIR, cdr_positions_heavy,
                        cdr_positions_light, disulfide_candidates,
                        existing_disulfides, in_variable_domain, load_residues)

res = load_residues()
print(f"Residues read from 1N8Z, chain A (VL) and B (VH): {len(res)}")

print("\nDISULFIDE BRIDGES THAT ALREADY EXIST — check that the method works")
print("-" * 62)
for a, b, d in existing_disulfides(res):
    print(f"   {a} – {b}     Cβ distance {d:.2f} Å")
print("   (if the method finds the known bridges, it will probably find the new one too)")

excl = ({("B", i) for i in cdr_positions_heavy()}
        | {("A", i) for i in cdr_positions_light()}
        | {("B", i) for i in HEAVY_KEEP})
alle = disulfide_candidates(res, exclude=excl)
# An scFv is VH + linker + VL. A pair in CL or CH1 does not exist there.
hits = [h for h in alle if in_variable_domain(h["a"]) and in_variable_domain(h["b"])]
print(f"\nOf {len(alle)} geometrically possible pairs, {len(hits)} lie in the")
print("variable domains. The remainder are in CL/CH1 and do not exist in an scFv.")

print(f"\nCANDIDATES for a NEW bridge between VL and VH: {len(hits)}")
print("The CDR residues are kept out — they bind and must be left alone.")
print("-" * 62)
print(f'{"VL":>7} {"VH":>7} {"Cβ":>7} {"Cα":>7}   assessment')


def is_known(h) -> bool:
    return (h["a"].resi, h["b"].resi) == KNOWN_ENGINEERED_PAIR


for h in hits[:12]:
    a, b = h["a"], h["b"]
    note = "very typical" if h["deviation"] < 0.3 else (
        "typical" if h["deviation"] < 0.6 else "at the edge")
    if is_known(h):
        note += "  <- built and measured by others, on 4D5"
    print(f'{str(a):>7} {str(b):>7} {h["cb_dist"]:>7.2f} {h["ca_dist"]:>7.2f}   {note}')

# Geometry alone cannot know which pair has been tried in a laboratory.
# We can. Kabat VH44 + VL100 is the most widely used engineered bridge between
# VH and VL, and Wörn & Plückthun (Biochemistry 1999) built exactly that one on 4D5,
# the antibody trastuzumab was humanized from: clearly more resistance to
# denaturation, binding retained, but the production yield dropped sharply.
# That is why we rank that pair first, and our own geometry finding second.
known = next((h for h in hits if is_known(h)), None)
best = known or (hits[0] if hits else None)

if best:
    a, b = best["a"], best["b"]
    print()
    print("RECOMMENDED PAIR" + ("  (known from the literature)" if known else "  (geometry only)"))
    print("-" * 62)
    print(f"   {a.resn}{a.resi} in the light chain  ->  CYS")
    print(f"   {b.resn}{b.resi} in the heavy chain  ->  CYS")
    print(f"   Cβ–Cβ {best['cb_dist']:.2f} Å, Cα–Cα {best['ca_dist']:.2f} Å")
    if known:
        print("   Kabat VH44 + VL100. The same bridge has been built on 4D5 (Wörn &")
        print("   Plückthun 1999) and on other Fvs (Reiter/Brinkmann/Pastan).")
        print("   Reported: clearly more stable, binding retained, lower yield.")
    print()
    print("   Novelty: both lie in the framework, so requirements A and B are untouched.")
    print("   Cysteines: two new ones, so the parity stays even. The hard filter")
    print("   against unpaired cysteines is still satisfied.")
    print()
    print("   THIS IS GEOMETRY PLUS LITERATURE, NOT A MEASUREMENT ON OUR DESIGN.")
    print("   It must be folded, and preferably measured, before anyone trusts it.")
