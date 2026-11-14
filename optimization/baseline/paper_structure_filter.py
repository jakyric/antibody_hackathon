"""Paper-faithful BLOSUM library design, phase 2: antibody structure filter and final selection.

  paper_structure_filter.py prepare   write one H/L FASTA per sampled design (nebius-runs/paper/inputs/)
  (then run baseline/abb3_nebius.sh to model them with ABodyBuilder3 on Nebius)
  paper_structure_filter.py analyze   apply the paper's filter and pick the final designs

The paper modelled designs with ABodyBuilder2 and removed those whose CDRH3 RMSD to the wild-type model was above
3.5 A or whose total CDR RMSD was above 1.5 A, then sampled equally across edit distances. We use ABodyBuilder3
(the successor from the same group) and the same two thresholds. RMSDs are over backbone N, CA, C atoms after
superposing the framework (non-CDR) C-alpha atoms of VH + VL, with IMGT CDR definitions from immunum.
"""
import glob
import sys

import immunum
import numpy as np
import pandas as pd

import common as C
import structure_tools as S

VH_T = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
        "SR{H3}YWGQGTLVTVSS")
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
IN = C.ROOT / "nebius-runs" / "paper" / "inputs"
RES = C.ROOT / "results" / "paper"
CDRH3_MAX, CDR_MAX = 3.5, 1.5


def cmd_prepare():
    mod = pd.read_csv(C.OUT / "paper_modelling_set.csv")
    IN.mkdir(parents=True, exist_ok=True)
    for r in mod.itertuples():
        (IN / f"{r.id}.fasta").write_text(f">H\n{VH_T.format(H3=r.cdrh3)}\n>L\n{VL}\n")
    print(f"wrote {len(mod)} FASTA files to {IN} (wt.fasta already exists)")


def regions(vh, vl):
    """Residue index sets (0-based within each chain): CDR residues and framework residues."""
    out = {}
    for name, seq, chains in (("H", vh, ["H"]), ("L", vl, ["K", "L"])):
        seg = immunum.Annotator(chains=chains, scheme="imgt", min_confidence=0.0).segment(seq)
        parts = [("prefix", seg.prefix), ("fr1", seg.fr1), ("cdr1", seg.cdr1), ("fr2", seg.fr2), ("cdr2", seg.cdr2),
                 ("fr3", seg.fr3), ("cdr3", seg.cdr3), ("fr4", seg.fr4)]
        pos, d = 0, {}
        for k, s in parts:
            d[k] = list(range(pos, pos + len(s)))
            pos += len(s)
        out[name] = d
    return out


def load_chain(pdb_path):
    """Return {chain: DataFrame(res_index, atom, xyz)} with residues numbered 0.. in file order."""
    rows = []
    for line in open(pdb_path):
        if line.startswith("ATOM"):
            rows.append((line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip(),
                         float(line[30:38]), float(line[38:46]), float(line[46:54])))
    df = pd.DataFrame(rows, columns=["chain", "res", "resname", "atom", "x", "y", "z"])
    out = {}
    for ch, g in df.groupby("chain"):
        order = {r: i for i, r in enumerate(sorted(g.res.unique()))}
        g = g.assign(i=g.res.map(order))
        out[ch] = g
    return out


def backbone(g, idxs):
    d = g[g.i.isin(idxs) & g.atom.isin(["N", "CA", "C"])].sort_values(["i", "atom"])
    return d[["x", "y", "z"]].to_numpy()


def ca(g, idxs):
    d = g[g.i.isin(idxs) & (g.atom == "CA")].sort_values("i")
    return d[["x", "y", "z"]].to_numpy()


def compare(wt_chains, ds_chains, reg):
    fw_h = [i for k in ("fr1", "fr2", "fr3", "fr4") for i in reg["H"][k]]
    fw_l = [i for k in ("fr1", "fr2", "fr3", "fr4") for i in reg["L"][k]]
    P = np.vstack([ca(ds_chains["H"], fw_h), ca(ds_chains["L"], fw_l)])
    Q = np.vstack([ca(wt_chains["H"], fw_h), ca(wt_chains["L"], fw_l)])
    R, t = S.kabsch(P, Q)

    def rm(chain, idxs):
        a = backbone(ds_chains[chain], idxs) @ R.T + t
        b = backbone(wt_chains[chain], idxs)
        return S.rmsd(a, b), len(b)
    h3, _ = rm("H", reg["H"]["cdr3"])
    allcdr_a, allcdr_b = [], []
    for ch in ("H", "L"):
        idx = [i for k in ("cdr1", "cdr2", "cdr3") for i in reg[ch][k]]
        allcdr_a.append(backbone(ds_chains[ch], idx) @ R.T + t)
        allcdr_b.append(backbone(wt_chains[ch], idx))
    tot = S.rmsd(np.vstack(allcdr_a), np.vstack(allcdr_b))
    return h3, tot


def cmd_analyze():
    mod = pd.read_csv(C.OUT / "paper_modelling_set.csv")
    reg = regions(VH_T.format(H3=C.WT), VL)
    wtp = glob.glob(str(RES / "abb3_wt" / "*.pdb"))[0]
    wt = load_chain(wtp)
    rows = []
    for r in mod.itertuples():
        pdbs = glob.glob(str(RES / "abb3" / r.id / "*.pdb"))
        if not pdbs:
            rows.append((r.id, np.nan, np.nan))
            continue
        h3, tot = compare(wt, load_chain(pdbs[0]), reg)
        rows.append((r.id, h3, tot))
    res = pd.DataFrame(rows, columns=["id", "abb3_cdrh3_rmsd", "abb3_cdr_rmsd"])
    mod = mod.merge(res, on="id")
    mod["structure_ok"] = (mod.abb3_cdrh3_rmsd <= CDRH3_MAX) & (mod.abb3_cdr_rmsd <= CDR_MAX)
    print(f"modelled {int(mod.abb3_cdrh3_rmsd.notna().sum())} of {len(mod)}; passing the paper's thresholds "
          f"(CDRH3 <= {CDRH3_MAX} A and total CDR <= {CDR_MAX} A): {int(mod.structure_ok.sum())}")
    print(mod[["id", "cdrh3", "edit_distance", "abb3_cdrh3_rmsd", "abb3_cdr_rmsd", "structure_ok"]].round(2).to_string(index=False))
    mod.to_csv(C.OUT / "paper_structure_filter.csv", index=False)
    # final: up to 2 passing designs per edit-distance bucket (paper: equal sampling across edit distances)
    ok = mod[mod.structure_ok]
    final = ok.groupby("edit_distance", group_keys=False).head(2).sort_values(["edit_distance", "id"]).reset_index(drop=True)
    final["n_changed"] = final.edit_distance
    cols = ["cdrh3", "changed_positions", "n_changed", "blosum_logp", "p_mean", "p_cnn", "p_hgb", "liabilities",
            "edit_distance", "abb3_cdrh3_rmsd", "abb3_cdr_rmsd"]
    final[cols].round(4).to_csv(C.OUT / "blosum_candidates_paper.csv", index=False)
    print(f"\nfinal selection: {len(final)} designs ({final.edit_distance.value_counts().sort_index().to_dict()})")
    print("saved baseline/outputs/blosum_candidates_paper.csv")


if __name__ == "__main__":
    {"prepare": cmd_prepare, "analyze": cmd_analyze}[sys.argv[1]]()
