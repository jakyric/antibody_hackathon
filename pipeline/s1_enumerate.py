"""Stage 1: enumerate CDR-H3 variants of the best tested loop with the BLOSUM-45 recipe.

Positions: non-contact CDR-H3 positions (config design.positions). Substitution probabilities follow the
paper's reverse-BLOSUM recipe, P(b|a) ~ f_b * exp(lambda * BLOSUM45(a, b)), computed from the *base loop* residue,
keeping the m most likely per position. Candidates change k positions of the base; they must
  - contain no Cys and no more sequence liabilities than trastuzumab,
  - be >= novelty_min_distance (Hamming) from every sequence in the provided data (abguard's rule).
Output: enum.csv sorted by BLOSUM joint probability.
"""
from __future__ import annotations

import csv
import itertools
import math
import re

from .common import (AAS, H3_COUNTS, WT_H3, blosum, fit_ridge, known_sequences, mutations, out, ridge_score,
                     save_summary)

LIABILITY_MOTIFS = ["NG", "NS", "NT", "NN", "DG", "DS", "DD", "DT", "DH", "DN", "DP", "DQ", "KE", "KD", "EK", "ED", "M", "W"]


def liabilities(seq: str) -> int:
    n = sum(len(re.findall(f"(?={m})", seq)) for m in LIABILITY_MOTIFS)
    n += len(re.findall(r"N[^P][ST]", seq)) + sum(seq.count(x) for x in ("RGD", "RYD", "LDV", "GPR"))
    return n


def substitution_table(base: str, positions: list[int], m: int, lam: float) -> dict[int, list[tuple[str, float]]]:
    B = blosum("BLOSUM45")
    tot = sum(H3_COUNTS.values())
    bg = {k: v / tot for k, v in H3_COUNTS.items()}
    table = {}
    for p in positions:
        a = base[p - 1]
        w = {b: bg[b] * math.exp(lam * B[a][b]) for b in AAS if b != a}
        z = sum(w.values())
        table[p] = sorted(((b, v / z) for b, v in w.items()), key=lambda kv: -kv[1])[:m]
    return table


def far_from_known(seq: str, known: set[str], min_dist: int) -> bool:
    """True if seq differs from every known sequence at >= min_dist positions (min_dist <= 3)."""
    if seq in known:
        return False
    for i in range(10):
        for a in AAS:
            if a != seq[i] and seq[:i] + a + seq[i + 1:] in known:
                return False
    if min_dist >= 3:
        for i, j in itertools.combinations(range(10), 2):
            for a in AAS:
                if a == seq[i]:
                    continue
                t = seq[:i] + a + seq[i + 1:]
                for b in AAS:
                    if b != seq[j] and t[:j] + b + t[j + 1:] in known:
                        return False
    return True


def run(cfg: dict) -> None:
    d = cfg["design"]
    base = d["base_h3"]
    table = substitution_table(base, d["positions"], d["top_substitutions"], d["blosum_lambda"])
    for p, subs in table.items():
        print(f"  pos {p} {base[p - 1]} -> " + " ".join(f"{b}:{v:.2f}" for b, v in subs))
    coef = fit_ridge(cfg)
    known = known_sequences(cfg)
    wt_liab = liabilities(WT_H3)

    raw = []
    for k in d["changes_from_base"]:
        for ps in itertools.combinations(d["positions"], k):
            for combo in itertools.product(*[table[p] for p in ps]):
                s = list(base)
                lp = 0.0
                for p, (b, pr) in zip(ps, combo):
                    s[p - 1] = b
                    lp += math.log(pr)
                raw.append(("".join(s), lp))
    no_cys = [r for r in raw if "C" not in r[0]]
    ok_liab = [r for r in no_cys if liabilities(r[0]) <= wt_liab]
    final = [r for r in ok_liab if far_from_known(r[0], known, d["novelty_min_distance"])]
    print(f"  raw {len(raw)} -> no Cys {len(no_cys)} -> liabilities <= WT {len(ok_liab)} -> distance >= {d['novelty_min_distance']}: {len(final)}")

    with open(out(cfg, "enum.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cdrh3", "blosum_logp", "ridge_sum", "liabilities", "mutations_vs_wt", "mutations_vs_base"])
        for s, lp in sorted(final, key=lambda r: -r[1]):
            w.writerow([s, f"{lp:.3f}", f"{ridge_score(s, coef):.3f}", liabilities(s),
                        ";".join(mutations(s, WT_H3)), ";".join(mutations(s, base))])
    save_summary(cfg, n_raw=len(raw), n_no_cys=len(no_cys), n_liab=len(ok_liab), n_enum=len(final))
