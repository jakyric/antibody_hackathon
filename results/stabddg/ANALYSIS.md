# StaB-ddG: what the existing run supports

Source: `results/stabddg/0a7e70d5-01ad-56f4-92e7-f44e9a12159c.csv`, operation `op-b3f8656a538edcc0`,
60 variants with both a model score and a measured KD.

## 1. What the run measured

| | |
| --- | --- |
| Output column | `binding_ddg_kcal_mol` |
| Partners | A = ['A', 'B'] (antibody) vs B = ['C'] (HER2) |
| Chains mutated | B |
| Residue numbers mutated | 99-108 (the CDR-H3 10-mer) |
| Light-chain mutations in the run | **0** |
| Monte Carlo samples | 20 |

**This is an interface free energy, not a folding free energy.** It is a
binding-axis quantity. It is not a melting temperature and does not
become one.

**The light chain is untouched in all 60 variants.**
VL R66G does not appear anywhere in this run. Nothing here can support,
or refute, a claim about the stabilising mutation.

## 2. Does it track measured affinity?

Sign convention: positive ddG = destabilised interface = weaker binding =
higher KD. A working model gives a positive correlation.

| Comparison | Spearman |
| --- | ---: |
| StaB-ddG vs log10 KD | **+0.378**  (95% CI [+0.14, +0.59], n=60) |
| mutation count vs log10 KD | +0.544 |
| StaB-ddG vs mutation count | +0.508 |

### The decisive test: does it work at fixed mutation count?

All ten submitted designs carry exactly 4 CDR-H3 mutations, so mutation
count is a constant within our set. Any skill that comes from counting
mutations is skill we cannot use. Restricting the panel to a narrow band:

| Mutation count | n | Spearman vs log10 KD | 95% CI |
| --- | ---: | ---: | --- |
| exactly 4 | 15 | +0.132 | [-0.41, +0.65] |
| 4-5 | 25 | +0.138 | [-0.26, +0.53] |
| 3-6 | 40 | +0.307 | [-0.00, +0.56] |

Median ddG of the 15 tightest binders: +0.421 kcal/mol
Median ddG of the 15 weakest binders: +0.890 kcal/mol

Score spread across all variants: 2.82 kcal/mol.
Median Monte Carlo standard error on a single score: 0.22 kcal/mol.

## 3. Verdict

The correlation with measured KD is positive and its 95% interval excludes zero. StaB-ddG carries real affinity signal on this panel.

But the score correlates with mutation count at +0.51, and mutation count alone predicts log KD at +0.54 -- BETTER than StaB-ddG's +0.38. On this panel the model does not beat simply counting how many residues changed.

At the mutation count our designs actually sit at (4-5, n=25), the correlation is +0.14 with 95% CI [-0.26, +0.53] -- it includes zero. We have no demonstrated ability to rank designs that all carry 4 mutations. Running StaB-ddG on our ten would produce numbers we could not justify ranking on.

## 4. What this means for the submission

- It is **not** evidence of a thermostability improvement for any design.
  It is a binding-axis measurement and must not be quoted on the
  thermostability axis.
- It says **nothing about VL R66G**, which is the only change in our
  designs that is claimed to affect thermostability.
- Spending more compute on StaB-ddG for the CDR-H3 loops would refine an
  axis we already rank on by other means, and would still leave the
  thermostability claim resting entirely on the published literature.
