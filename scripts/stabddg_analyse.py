"""What the existing StaB-ddG run actually shows -- and what it cannot show.

READ THIS BEFORE CITING THE NUMBERS

  `results/stabddg/` holds one completed StaB-ddG run, 60 variants. It is a
  VALIDATION run, not a scoring run: `scripts/stabddg_prep.py` picked 60
  trastuzumab CDR-H3 variants that already have a measured KD in the provided
  SPR table, spanning the whole KD range, so that the model's output could be
  checked against a known answer before we trusted it on anything.

  Two facts about that run decide what it can be used for:

  1. The column is `binding_ddg_kcal_mol`, with partner A = {antibody light,
     antibody heavy} and partner B = {HER2}. That is the free energy of the
     ANTIBODY-ANTIGEN INTERFACE. It is not a folding free energy and not a
     melting temperature. It speaks to the binding axis, not the
     thermostability axis.

  2. Every mutation in the run is on chain B, residue numbers 99-108: the
     CDR-H3 10-mer. The light chain (chain A) is untouched in all 60 variants.
     The run therefore contains no information whatsoever about VL R66G.

  So: this run cannot support a thermostability claim for any design, and it
  cannot say anything about the stabilising mutation. What it CAN do is tell
  us whether StaB-ddG's binding ddG tracks measured affinity well enough to
  rank loops with. This script answers that one question.

  Sign convention: a POSITIVE binding ddG means the mutation destabilises the
  interface, i.e. weaker binding, i.e. a HIGHER KD. So if the model works we
  expect a POSITIVE correlation between ddG and KD.

USAGE
  python scripts/stabddg_analyse.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def spearman(x: list[float], y: list[float]) -> float:
    def rank(v: list[float]) -> list[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def bootstrap_ci(x: list[float], y: list[float], n_boot: int = 2000, seed: int = 0):
    import random

    rng = random.Random(seed)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        vals.append(spearman([x[i] for i in idx], [y[i] for i in idx]))
    vals.sort()
    return vals[int(0.025 * n_boot)], vals[int(0.975 * n_boot)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--scores",
        type=Path,
        default=ROOT / "results/stabddg/0a7e70d5-01ad-56f4-92e7-f44e9a12159c.csv",
    )
    ap.add_argument(
        "--truth", type=Path, default=ROOT / "nebius-runs/stabddg-validate/truth.json"
    )
    ap.add_argument(
        "--params",
        type=Path,
        default=ROOT / "results/stabddg/225bade5-bf39-5cb0-af46-7996aa5d4c71.json",
    )
    ap.add_argument("--out", type=Path, default=ROOT / "results/stabddg/ANALYSIS.md")
    args = ap.parse_args(argv)

    scores = {r["variant_id"]: r for r in csv.DictReader(args.scores.open())}
    truth = {t["id"]: t for t in json.loads(args.truth.read_text())}
    params = json.loads(args.params.read_text())["parameters"]

    # --- what the run covers, stated from the run itself, not from memory ---
    chains = sorted({m["chain"] for v in params["variants"] for m in v["mutations"]})
    positions = sorted({m["residue_number"] for v in params["variants"] for m in v["mutations"]})
    light_chain_mutations = sum(
        1 for v in params["variants"] for m in v["mutations"] if m["chain"] == "A"
    )

    ids = [i for i in scores if i in truth]
    ddg = [float(scores[i]["binding_ddg_kcal_mol"]) for i in ids]
    err = [float(scores[i]["mc_standard_error_kcal_mol"]) for i in ids]
    kd = [truth[i]["kd_nM"] for i in ids]
    logkd = [math.log10(k) for k in kd]
    nmut = [truth[i]["n_mut"] for i in ids]

    rho = spearman(ddg, logkd)
    lo, hi = bootstrap_ci(ddg, logkd)
    rho_nmut = spearman(nmut, logkd)
    rho_ddg_nmut = spearman(ddg, nmut)

    # THE DECISIVE TEST.
    #   All ten submitted designs carry exactly 4 CDR-H3 mutations. Within our
    #   candidate set, mutation count is a constant. So a model whose apparent
    #   skill comes from counting mutations has zero skill where we need it.
    #   Restrict to a narrow band of mutation counts and ask again.
    bands = {}
    for label, lo_n, hi_n in [("exactly 4", 4, 4), ("4-5", 4, 5), ("3-6", 3, 6)]:
        sel = [k for k in range(len(ids)) if lo_n <= nmut[k] <= hi_n]
        if len(sel) >= 8:
            r = spearman([ddg[k] for k in sel], [logkd[k] for k in sel])
            clo, chi = bootstrap_ci([ddg[k] for k in sel], [logkd[k] for k in sel])
            bands[label] = (len(sel), r, clo, chi)
        else:
            bands[label] = (len(sel), float("nan"), float("nan"), float("nan"))

    # Does it separate the best binders from the worst?
    pairs = sorted(zip(kd, ddg))
    tight = [d for k, d in pairs[:15]]
    weak = [d for k, d in pairs[-15:]]

    # Is the spread even real, given Monte Carlo noise?
    spread = max(ddg) - min(ddg)
    median_err = sorted(err)[len(err) // 2]

    lines = [
        "# StaB-ddG: what the existing run supports",
        "",
        f"Source: `{args.scores.relative_to(ROOT)}`, operation `op-b3f8656a538edcc0`,",
        f"{len(ids)} variants with both a model score and a measured KD.",
        "",
        "## 1. What the run measured",
        "",
        "| | |",
        "| --- | --- |",
        "| Output column | `binding_ddg_kcal_mol` |",
        f"| Partners | A = {params['partner_a']} (antibody) vs B = {params['partner_b']} (HER2) |",
        f"| Chains mutated | {', '.join(chains)} |",
        f"| Residue numbers mutated | {min(positions)}-{max(positions)} (the CDR-H3 10-mer) |",
        f"| Light-chain mutations in the run | **{light_chain_mutations}** |",
        f"| Monte Carlo samples | {params.get('mc_samples')} |",
        "",
        "**This is an interface free energy, not a folding free energy.** It is a",
        "binding-axis quantity. It is not a melting temperature and does not",
        "become one.",
        "",
        f"**The light chain is untouched in all {len(params['variants'])} variants.**",
        "VL R66G does not appear anywhere in this run. Nothing here can support,",
        "or refute, a claim about the stabilising mutation.",
        "",
        "## 2. Does it track measured affinity?",
        "",
        "Sign convention: positive ddG = destabilised interface = weaker binding =",
        "higher KD. A working model gives a positive correlation.",
        "",
        "| Comparison | Spearman |",
        "| --- | ---: |",
        f"| StaB-ddG vs log10 KD | **{rho:+.3f}**  (95% CI [{lo:+.2f}, {hi:+.2f}], n={len(ids)}) |",
        f"| mutation count vs log10 KD | {rho_nmut:+.3f} |",
        f"| StaB-ddG vs mutation count | {rho_ddg_nmut:+.3f} |",
        "",
        "### The decisive test: does it work at fixed mutation count?",
        "",
        "All ten submitted designs carry exactly 4 CDR-H3 mutations, so mutation",
        "count is a constant within our set. Any skill that comes from counting",
        "mutations is skill we cannot use. Restricting the panel to a narrow band:",
        "",
        "| Mutation count | n | Spearman vs log10 KD | 95% CI |",
        "| --- | ---: | ---: | --- |",
    ] + [
        f"| {lab} | {n} | {r:+.3f} | [{c0:+.2f}, {c1:+.2f}] |"
        for lab, (n, r, c0, c1) in bands.items()
    ] + [
        "",
        f"Median ddG of the 15 tightest binders: {sum(tight)/len(tight):+.3f} kcal/mol",
        f"Median ddG of the 15 weakest binders: {sum(weak)/len(weak):+.3f} kcal/mol",
        "",
        f"Score spread across all variants: {spread:.2f} kcal/mol.",
        f"Median Monte Carlo standard error on a single score: {median_err:.2f} kcal/mol.",
        "",
    ]

    verdict = []
    if lo > 0:
        verdict.append(
            "The correlation with measured KD is positive and its 95% interval "
            "excludes zero. StaB-ddG carries real affinity signal on this panel."
        )
    elif hi < 0:
        verdict.append(
            "The correlation runs the WRONG WAY and the interval excludes zero: "
            "the model's ddG anti-predicts measured affinity here. Do not rank with it."
        )
    else:
        verdict.append(
            "The 95% interval includes zero. On this panel StaB-ddG does not "
            "demonstrably track measured affinity, so it cannot be used to rank "
            "designs by binding. That is a negative result, and it is worth "
            "reporting as one."
        )
    if abs(rho_ddg_nmut) > 0.4:
        verdict.append(
            f"But the score correlates with mutation count at {rho_ddg_nmut:+.2f}, and "
            f"mutation count alone predicts log KD at {rho_nmut:+.2f} -- BETTER than "
            f"StaB-ddG's {rho:+.2f}. On this panel the model does not beat simply "
            "counting how many residues changed."
        )
    n4, r4, c40, c41 = bands["4-5"]
    if n4 >= 8 and not (c40 > 0 or c41 < 0):
        verdict.append(
            f"At the mutation count our designs actually sit at (4-5, n={n4}), the "
            f"correlation is {r4:+.2f} with 95% CI [{c40:+.2f}, {c41:+.2f}] -- it "
            "includes zero. We have no demonstrated ability to rank designs that "
            "all carry 4 mutations. Running StaB-ddG on our ten would produce "
            "numbers we could not justify ranking on."
        )
    elif n4 >= 8 and c40 > 0:
        verdict.append(
            f"At 4-5 mutations (n={n4}) the correlation holds at {r4:+.2f}, CI "
            f"[{c40:+.2f}, {c41:+.2f}]. The signal survives at the mutation count "
            "our designs sit at, so StaB-ddG is worth running on the ten."
        )
    lines += ["## 3. Verdict", ""] + [f"{v}\n" for v in verdict]
    lines += [
        "## 4. What this means for the submission",
        "",
        "- It is **not** evidence of a thermostability improvement for any design.",
        "  It is a binding-axis measurement and must not be quoted on the",
        "  thermostability axis.",
        "- It says **nothing about VL R66G**, which is the only change in our",
        "  designs that is claimed to affect thermostability.",
        "- Spending more compute on StaB-ddG for the CDR-H3 loops would refine an",
        "  axis we already rank on by other means, and would still leave the",
        "  thermostability claim resting entirely on the published literature.",
        "",
    ]

    args.out.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
