"""Build the shared HER2 domain IV target files from the 1N8Z crystal (trastuzumab Fv + HER2 domain IV).

Run:  denovo/.venv/bin/python denovo/targets/make_targets.py
Writes (in denovo/targets/): ref_1n8z.cif, her2_domIV.cif, her2_domIV.fasta, her2_domIV.a3m, numbering_map.json,
epitope.json, README.md, READY.

Crystal chains in data/her2/her2_fv_1n8z.cif: A = VL, B = VH, C = HER2 489-580, D = HER2 591-607 (581-590 unresolved).
"""
import json
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "optimization" / "data" / "her2" / "her2_fv_1n8z.cif"
OUT = Path(__file__).resolve().parent
A3M = ROOT / "optimization" / "baseline" / "refold_inputs" / "her2_domIV.a3m"
CUT = 4.5
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}

lines = SRC.read_text().splitlines()
hdr = [l.strip().split(".")[1] for l in lines if l.startswith("_atom_site.")]
atoms = []
for l in lines:
    if l.startswith(("ATOM", "HETATM")):
        t = l.split()
        d = dict(zip(hdr, t))
        if d["type_symbol"] == "H":
            continue
        atoms.append(dict(chain=d["label_asym_id"], label_seq=int(d["label_seq_id"]), auth_seq=int(d["auth_seq_id"]),
                          resname=d["label_comp_id"], atom=d["label_atom_id"].strip('"'),
                          xyz=np.array([float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])]), line=l))
print(len(atoms), "heavy atoms; chains", sorted({a["chain"] for a in atoms}))

her2 = [a for a in atoms if a["chain"] in ("C", "D")]
fv = [a for a in atoms if a["chain"] in ("A", "B")]
H = np.array([a["xyz"] for a in her2])
F = np.array([a["xyz"] for a in fv])
D = np.linalg.norm(H[:, None, :] - F[None, :, :], axis=2)

per = defaultdict(lambda: dict(n=0, dmin=9e9, atom=None, fv_chain=defaultdict(int), fv_res=set(), atoms=defaultdict(int)))
for i, j in zip(*np.where(D <= CUT)):
    a, b = her2[i], fv[j]
    r = per[a["auth_seq"]]
    r["resname"] = a["resname"]
    r["n"] += 1
    r["fv_chain"][b["chain"]] += 1
    r["fv_res"].add((b["chain"], b["auth_seq"]))
    r["atoms"][a["atom"]] += 1
    if D[i, j] < r["dmin"]:
        r["dmin"], r["atom"] = float(D[i, j]), a["atom"]

# --- clean domain IV target: one chain "A", contiguous label_seq_id over modelled residues
seen, order = set(), []
for a in her2:
    k = (a["chain"], a["auth_seq"])
    if k not in seen:
        seen.add(k)
        order.append(k)
order.sort(key=lambda k: k[1])
label_of = {k[1]: i + 1 for i, k in enumerate(order)}
seq = "".join(THREE[next(a["resname"] for a in her2 if a["auth_seq"] == k[1])] for k in order)
numbering = {str(label_of[k[1]]): k[1] for k in order}
(OUT / "numbering_map.json").write_text(json.dumps({"label_to_auth": numbering, "note": "HER2 numbering; 581-590 unresolved"}, indent=1))
(OUT / "her2_domIV.fasta").write_text(">her2_domIV_modelled_1N8Z\n" + seq + "\n")

cif = ["data_her2_domIV", "#", "loop_", "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
       "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
       "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x",
       "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
       "_atom_site.auth_seq_id", "_atom_site.auth_comp_id", "_atom_site.auth_asym_id", "_atom_site.auth_atom_id",
       "_atom_site.pdbx_PDB_model_num"]
n = 0
for a in sorted(her2, key=lambda a: (a["auth_seq"],)):
    n += 1
    el = a["atom"][0]
    x, y, z = a["xyz"]
    cif.append(f"ATOM {n} {el} {a['atom']} . {a['resname']} A 1 {label_of[a['auth_seq']]} ? {x:.3f} {y:.3f} {z:.3f} "
                f"1.00 0.00 {a['auth_seq']} {a['resname']} A {a['atom']} 1")
(OUT / "her2_domIV.cif").write_text("\n".join(cif) + "\n")
shutil.copy(SRC, OUT / "ref_1n8z.cif")
shutil.copy(A3M, OUT / "her2_domIV.a3m")

# --- epitope + hotspot sets
ep = []
for auth, r in sorted(per.items()):
    ep.append(dict(auth_seq=auth, label_seq=label_of[auth], resname=r["resname"], aa=THREE[r["resname"]],
                   n_contacts=r["n"], min_dist=round(r["dmin"], 2), closest_atom=r["atom"],
                   fv_chains=dict(r["fv_chain"]), frac_VH=round(r["fv_chain"].get("B", 0) / r["n"], 2)))
ep_sorted = sorted(ep, key=lambda e: -e["n_contacts"])
ca = {a["auth_seq"]: a["xyz"] for a in her2 if a["atom"] == "CA"}
res_ids = [e["auth_seq"] for e in ep]
P = np.array([ca[i] for i in res_ids])

# 2-means on epitope CA to split the epitope into two patches (deterministic init: farthest pair)
d = np.linalg.norm(P[:, None] - P[None], axis=2)
i0, j0 = np.unravel_index(np.argmax(d), d.shape)
c = np.array([P[i0], P[j0]])
for _ in range(50):
    lab = np.argmin(np.linalg.norm(P[:, None] - c[None], axis=2), axis=1)
    c = np.array([P[lab == k].mean(0) for k in range(2)])
patch = {res_ids[i]: int(lab[i]) for i in range(len(res_ids))}


def hot(e):
    return dict(auth_seq=e["auth_seq"], label_seq=e["label_seq"], resname=e["resname"], atom=e["closest_atom"])


by_patch = {k: sorted([e for e in ep if patch[e["auth_seq"]] == k], key=lambda e: -e["n_contacts"]) for k in (0, 1)}
big = max(by_patch, key=lambda k: sum(e["n_contacts"] for e in by_patch[k]))
small = 1 - big
sets = {
    "hs1_core": [hot(e) for e in by_patch[big][:6]],           # most-contacted residues of the dominant patch
    "hs2_broad": [hot(e) for e in ep_sorted[:12]],             # twelve most-contacted residues overall
    "hs3_alt_patch": [hot(e) for e in by_patch[small][:6]],    # the other patch (diversity / flank)
}
doc = dict(source="1N8Z", cutoff_A=CUT, n_epitope_res=len(ep), n_contact_pairs=int(sum(e["n_contacts"] for e in ep)),
           residues=ep, patches={str(k): [e["auth_seq"] for e in by_patch[k]] for k in by_patch},
           hotspot_sets=sets, numbering="auth_seq = HER2 numbering used in the paper; label_seq = contiguous 1..N in her2_domIV.cif")
(OUT / "epitope.json").write_text(json.dumps(doc, indent=1))

(OUT / "README.md").write_text(f"""# Targets
* `her2_domIV.cif`: HER2 domain IV as ONE chain `A`, {len(order)} modelled residues (489-580, 591-607 of 1N8Z; 581-590 unresolved).
  `label_seq_id` is contiguous 1..{len(order)}; `auth_seq_id` keeps HER2 numbering. Map in `numbering_map.json`.
* `her2_domIV.fasta`: the {len(seq)} modelled residues. (The full 119 aa domain IV sequence for folding is in `her2_domIV.a3m` query.)
* `her2_domIV.a3m`: existing MSA for domain IV (119 aa, 489-607) used by the earlier refolds.
* `ref_1n8z.cif`: the crystal complex, chains A=VL, B=VH, C+D=HER2.
* `epitope.json`: HER2 residues within {CUT} A of trastuzumab Fv heavy atoms ({len(ep)} residues, from the crystal, not
  from literature lists) and three hotspot sets: hs1_core (6 most-contacted of the dominant patch), hs2_broad (12 most-contacted),
  hs3_alt_patch (6 of the other patch).
Regenerate with `denovo/.venv/bin/python denovo/targets/make_targets.py`.
""")
(OUT / "READY").write_text("ok\n")
print("epitope residues:", len(ep), "patches:", {k: [e['auth_seq'] for e in v] for k, v in by_patch.items()})
for k, v in sets.items():
    print(k, [(h["resname"] + str(h["auth_seq"]), h["atom"]) for h in v])
print("seq", len(seq), seq)
