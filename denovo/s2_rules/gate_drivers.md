# What drives the abguard "identity to trastuzumab" gate for a VHH

Scope: the track decision of the VM abguard (build 0.1.0, 2026-11-11 16:17, immunum 2.0.0). Every number below is **measured** by
`gate_drivers_analysis.py` on the designs in `results/denovo/designs_*.csv` at 14:18 CET (136 RFD4 designs on the BCII10 framework, 64 on the
caplacizumab-like framework, 115 BoltzGen; per-design data in `gate_drivers_data.csv`, summary in `gate_drivers_stats.json`) unless marked
*read* (from the abguard source, `abguard/novelty.py`) or *inferred*. Our local replica reproduces the VM value **exactly** (max |difference| 0.0000
over 430 designs with identical sequences).

## 1. Exactly how the number is computed *(read + measured)*

```
identity = 1 - Levenshtein(Fv_design, Fv_trastuzumab) / max(len(Fv_design), len(Fv_trastuzumab))      # novelty.py: _identity_to_trastuzumab
class    = 'novel'   if identity >= 0.70          (IDENTITY_THRESHOLD = 0.70, ">=")                  # engineered trastuzumab variant, needs >= 3 CDRH3 mutations
           'de_novo' if identity <  0.70                                                             # needs >= 2 CDRH3 mutations from trastuzumab and from every reference CDRH3
```

* **Alignment type:** the *global* unit-cost edit distance (substitution, insertion, deletion all cost 1; `polars-distance` `dist_str.levenshtein`) between two
  strings. It is not a position-by-position comparison and not a local alignment: indels are allowed anywhere and cost the same as a mismatch.
* **Which residues count:** the whole numbered variable domain = `fr1 + cdr1 + fr2 + cdr2 + fr3 + cdr3 + fr4` as segmented by immunum (IMGT, chain H, confidence >= 0.1).
  Tags, linkers and residues before IMGT 1 / after IMGT 128 (`prefix`/`postfix`) are removed. CDRs count exactly like framework; there is no CDR masking or weighting.
* **Reference:** trastuzumab VH (DrugBank DB00072), `EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS`,
  numbered the same way: **120 residues** (CDR1 8, CDR2 8, CDR3 13 = `SRWGGDGFYAMDY`).
* **Denominator:** `max(len(design Fv), 120)`. A VHH has no light chain, so abguard's pooled formula has only the heavy terms (the VL is not charged).
* **Threshold:** a design is `de_novo` iff `Levenshtein distance > 0.30 * max(len, 120)`; for a 125-residue Fv that means 38 or more edits.
  `identity_to_trastuzumab` in `rules.csv` and `abguard_identity_trast` agree to 4 decimals.

## 2. Why length moves it *(measured)*

Write `m = L - d` (length minus edit distance, the "matched-equivalent" count). Then `identity = m / L` for `L >= 120`. Across the 136 RFD4 designs on the BCII10
framework **`m` is almost constant: 87.2 residues, sd 1.8, while `L` varies with sd 3.2** (range 119-132). So identity is, to first order, `87.2 / L`, and the 70 % line
sits at **L = 87.2 / 0.70 = 124.6**. Where the 87 comes from (Levenshtein alignment of the BCII10 designs, means over designs): 79.0 framework residues identical to trastuzumab's
(BCII10 shares 80 of 91 framework positions with trastuzumab at the same IMGT label, 0.88; the caplacizumab-like framework 75/91, 0.82) + 10.7 CDR residues that happen to align
with an identical residue in trastuzumab's CDRs, minus insertion/deletion bookkeeping. Every extra residue (any CDR, any content) adds 1 to the denominator and 1 edit, so it
lowers identity by about `m/L^2 = 0.0056`. Measured slope: **-0.004 per Fv residue** (regression over the 136 designs), **-0.0036 per CDR3 residue** (sweep: BCII10 framework + its own CDR1/CDR2
+ random CDR3 of 14 -> 24 residues: 0.6965 -> 0.6602).

| BCII10 framework, RFD4 (n = 136) | Fv length 119-120 | 121-122 | 123-124 | 125-126 | 127-128 | >= 129 |
|---|---|---|---|---|---|---|
| designs | 9 | 22 | 23 | 27 | 31 | 24 |
| `de_novo` (identity < 0.70) | 0 | 0 | 6 | 19 | 24 | 22 |
| median identity | 0.708 | 0.712 | 0.707 | 0.691 | 0.685 | 0.673 |

| BCII10, by immunum CDR3 length | <= 13 | 14-15 | 16-17 | 18-19 | >= 20 |
|---|---|---|---|---|---|
| designs | 23 | 31 | 39 | 22 | 21 |
| `de_novo` | 1 | 11 | 21 | 19 | 19 |

The pattern agrees with the length table the orchestrator measured on the first 96 designs (total length <= 122: 0 %, >= 128: ~90 %); CDR3 lengths differ by about 2 residues between S1's contig
block and immunum's IMGT CDR3 (immunum includes the flanking `AA`/`DY`-type residues), so the bins are not the same. CDR3 length is the part of `L` that varies most (r(L, CDR3) = 0.87; CDR1 0.44, CDR2 0.28).
The **framework decides where the line is**: the caplacizumab-like framework has `m = 81.0 +- 2.0`, so the crossing is at L = 115.7 and 63 of 64 designs (98 %) are `de_novo` at any
length (Fv 116-128). With the BCII10 framework (a humanised VH3-like VHH) the line sits in the middle of the length distribution, which is why the class "flips on single residues".

Worked examples (BCII10 framework): `bcii10-rfd4jf-s1400-000` has L = 119 (CDR 8/9/12): 79 framework + 12 CDR residues identical, 26 substitutions, 2 deletions, 3 insertions, d = 31,
identity 0.742 -> `novel`. `bcii10-rfd4jf-s1301-004`, L = 132 (CDR3 21): same 79 framework residues + 9 CDR matches, 30 substitutions, **14 deletions** (the extra CDR3 residues have no partner), d = 46, identity
0.6515 -> `de_novo`.

## 3. Content or only length? *(measured)*

* Replace every CDR residue of a design by `X` (a letter that can never match trastuzumab): identity then correlates with `L` at **r = -0.99**; with random amino acids at **r = -0.94**.
  Length alone therefore sets the "content-free" identity.
* Real designs vs the same designs with random CDR residues: **+0.042**; vs CDRs that cannot match: **+0.075**. Designed CDRs resemble trastuzumab's CDRs a little more than random sequence does,
  because Y, G, S, D, A recur at aligned positions.
* Explained variance of identity (BCII10 RFD4 designs, sd of identity = 0.020): Fv length alone **R^2 = 0.54**; Fv length + number of aligned CDR residues identical to trastuzumab's **R^2 = 0.84**
  (r(identity, CDR matches) = 0.54). So length is the main knob, CDR content is a real second-order effect worth about +-0.02 identity, which is exactly the width of the uncertain zone around 0.70.
* Content only becomes a first-order effect when it *copies trastuzumab*: the BCII10 framework with trastuzumab's own CDR3 (13 residues) gives 0.797, with a 22-residue CDR3 that starts with it still 0.742.
  (BCII10 wild type, natural CDRs, L = 128: 0.680 `de_novo`. BCII10 + random CDR3 of 8-13 residues: 0.693-0.698, flat; *inferred* reason: below 13 residues each missing CDR3 residue is a deletion edit while the denominator is already pinned at 120.)

## 4. "New CDRs land below 70 %": a length effect? *(measured; for the BCII10 framework, mostly yes)*

For the BCII10 framework the statement is **not true by itself**: 65 of 136 new-CDR designs (48 %) are `novel` (identity >= 0.70), and *none* of the 31 designs with Fv length <= 122 is `de_novo`, whatever their
CDR content. What pushes designs below 0.70 is Fv length >= ~125 (a long CDR3), i.e. the denominator effect of section 2, plus a few percentage points of luck in how many CDR residues coincide with trastuzumab's. The
wild-type BCII10, which is a *known* VHH, is `de_novo` (0.680) purely because its CDR3 is 18 residues long. With the caplacizumab-like framework the statement is true (98 %), but because that framework is less
similar to trastuzumab VH (82 % framework identity vs 88 %), not because the CDRs are new.

## 5. `identity_nearest_known` is a different quantity *(measured)*

Same formula, but the maximum over 183,567 known heavy domains instead of trastuzumab alone. For these designs the framework is the parent's, so the nearest entries are BCII10 relatives (BCII10 designs: `9vq6_A` for 80/136,
`9nh7_E` 23, `3eak_A` 15), and `identity_nearest_known` = 1 - (edits in the redesigned CDRs and their indels) / L = **0.807 +- 0.022**, range 0.76-0.87. It is driven by how much of the parent scaffold plus CDR1/CDR2 pattern is
retained, not by similarity to trastuzumab (r with trastuzumab identity 0.61, with Fv length -0.73: longer CDR3s move *away* from the 18-residue parent loop). **It does not enter the VM abguard decision at all**; it
estimates the jury sentence "under 70 % identity to any known antibody" (PDF), which this abguard build does not implement. With the framework intact, reaching < 0.70 on that quantity is impossible without
changing about 20 framework residues (`decisions.md` section 4), so for the intact-framework strategy the only gate we can verify and optimise is the trastuzumab-identity one; the nearest-known number remains a
disclosed, unverifiable risk (0 of 144 BCII10 and 0 of 72 caplacizumab-like designs are below 0.70, 14:20 CET; maxima 0.87 and 0.85).

## 6. What this means for selection (framework intact)

* Target **identity <= 0.68** (margin of 0.02), not just < 0.70: 29 of 136 BCII10 designs satisfy it, 42 more sit in 0.68-0.70 where an immunum boundary shift or a single residue decides the class, 52 are in 0.70-0.72 and 13 above
  0.72. Caplacizumab-like: 51 of 64 are <= 0.68.
* For BCII10 designs, Fv length 125-126 gives 70 % `de_novo`, 127-128 77 %, >= 129 92 %; length <= 122 gives 0 %. By immunum CDR3 length: <= 13 -> 4 %, 14-15 -> 35 %, 16-17 -> 54 %, 18-19 -> 86 %, >= 20 -> 90 %. Asking S1 for CDR3
  lengths 18-20 would raise the pass rate for free, but long CDR3s are the least reliably predicted loops (check S3's `of3_cdr3_plddt`) and every extra residue is a modelling risk.
* The class is a property of the final sequence, so check `denovo_vm_pass` in `rules.csv` for every candidate (the watcher keeps it current); never infer it from the generator.
