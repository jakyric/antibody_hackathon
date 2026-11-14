# S2: rule harness and framework diversification

> **STATUS (14:15 CET, decision of the user and organisers): the framework stays intact; framework diversification is NOT RUN.**
> `designs_fwdiv.csv` was moved to `results/denovo/parked_fwdiv.csv` (outside the `designs_*.csv` glob, so S3/S4 do not spend on it) and the
> `diversify.py --watch` process was stopped. `diversify.py` stays in the repo as a documented fallback only. **`rules.py --watch` is the gate** (`denovo_vm_pass`,
> `identity_to_trastuzumab`); how that number behaves is explained in [`gate_drivers.md`](gate_drivers.md).
>
> If diversification is ever forced: mutate CDR-adjacent framework positions first, never IMGT 42 nor the C-A-A at the FR3/CDR3 junction (Cys104-Ala105-Ala106),
> keep the camelid/human FR2 hallmarks (IMGT 42/49/50/52) and every Cys; use the `lite` tier (<= 8 natural substitutions, goal identity to trastuzumab VH <= 0.66) before `deep`.
> The code protects IMGT 23/41/42/49/50/52/104/118/119/121 and 26/39/55/66 and all CDR positions already; the "CDR-adjacent first" ordering is a recommendation, not implemented.

Two jobs, both resumable and both running as watchers during the hackathon:

1. `rules.py` checks every design against the hackathon's hard rules and the novelty rule, with the **real abguard on the team VM**
   plus a local replica, and writes `results/denovo/rules.csv` (one row per `design_id`, `True`/`False` booleans, numeric identities).
2. `diversify.py` takes the best scored designs (from S3) and writes framework-diversified variants to
   `results/denovo/designs_fwdiv.csv` (generator `fwdiv`, `parent_id` set) so S3 refolds them and S4 profiles them as first-class designs.

Why this stage exists: whatever generator we use, VHH frameworks stay near the parent's (RFD4 and BoltzGen keep or re-use a fixed
scaffold). Measured on the first 96 real RFD4 designs (BCII10 framework, new CDR1-3): **0/96 are below 70 % whole-Fv identity to the nearest
known antibody** (median 0.80, range 0.76-0.87), and only 51/96 are classed `de_novo` by the VM abguard (identity to trastuzumab VH < 70 %,
median 0.696). See `decisions.md` for the evidence and what abguard really checks.

## Setup (once; everything already done on this machine)

```bash
cd denovo/s2_rules
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python immunum==2.0.0 rapidfuzz numpy pandas biopython requests pyarrow   # immunum 2.0.0 = the abguard pin
.venv/bin/python build_panel.py      # PDB seqres (wwPDB, 14 Nov 2026) heavy domains + UniProt human IGHV germlines  (cache/pdb_seqres.txt.gz is downloaded first:
                                     #   curl -o cache/pdb_seqres.txt.gz https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz)
# PLAbDab downloads (opig.stats.ox.ac.uk/webapps/plabdab[-nano]/static/downloads/{vhh,unpaired,paired}_sequences.csv.gz) -> cache/plabdab/
.venv/bin/python build_panel2.py     # + PLAbDab(-nano) heavy domains; result cache/panel.csv (183,567 distinct Fv), backup cache/panel_v1.csv
```

abguard reference CDRH3s (610,109) are copied from the VM wheel: `cache/abguard_cdrh3.parquet`
(`~/.cache/uv/archive-v0/Oobkm_d4jbq4a0Sw/abguard/reference_data/`). Real abguard runs on the VM:
`~/.cache/uv/archive-v0/jsRWGsUwPztt-JwD/bin/abguard` (build of 2026-11-11 16:17, version 0.1.0, immunum 2.0.0), work dir `~/moritz/denovo_s2`.

## Run

```bash
PY=denovo/.venv/bin/python
$PY denovo/s2_rules/rules.py                        # one-shot over all results/denovo/designs_*.csv -> results/denovo/rules.csv
$PY denovo/s2_rules/rules.py --no-vm                # same without ssh (abguard_* columns empty)
$PY denovo/s2_rules/rules.py --watch --interval 120 --until 15:45          # incremental watcher (new/changed design_ids only)
# NOT RUN (see status): $PY denovo/s2_rules/diversify.py --watch --interval 240 --until 15:30   # parents from results/denovo/score_of3.csv, lite + deep tiers
$PY denovo/s2_rules/diversify.py --standin                                  # scratch test on stand-in VHHs (writes only to cache/)
$PY denovo/s2_rules/diversify.py --in results/denovo/designs_rfd4jf.csv --ids-file ids.txt --out cache/x.csv --report cache/x_rep.csv   # scratch on chosen designs
$PY denovo/s2_rules/fwdiv_eval.py                   # which settings keep the parent's OF3 signal -> fwdiv_eval.md
$PY denovo/s2_rules/rules_summary.py                # per-generator table of rule results
$PY denovo/s2_rules/make_table.py denovo/s2_rules/fwdiv_report.csv > denovo/s2_rules/before_after.md
bash denovo/s2_rules/run_s2.sh status|stop|stop-all    # (`start` also launches the diversify watcher: do not use it now; start rules only, see below)
$PY denovo/s2_rules/rules.py --watch --interval 120 --until 15:45 >> denovo/s2_rules/watch.log 2>&1 &   # the only watcher that should run
# watchers (nohup, log denovo/s2_rules/watch.log); `touch denovo/STOP` also stops both
```

Both watchers exit on `denovo/STOP`, or at 15:45 (rules) / 15:30 (diversify) local time, write every output atomically (tmp + rename),
survive a failed cycle (log + retry), and only process new/changed rows.

## rules.csv columns (all booleans are `True`/`False`)

| group | columns |
|---|---|
| local hard rules | `len_ok` (<= 250), `invalid_aa`, `cys_count`, `noncanon_cys`, `unpaired_cys` (abguard's definition: odd count), `canon_cys_ok` (Cys at IMGT 23 and 104), `numbering_ok` (immunum 2.0.0, chain H, confidence >= 0.1), `imm_conf`, `cdrh3` (IMGT), `hallmarks` (IMGT 42/49/50/52), `cdrh3_dist_trast` (min Levenshtein to the IMGT 13-mer and the 10-mer core of trastuzumab's CDRH3) |
| abguard replica | `identity_to_trastuzumab` (abguard formula: 1 - Levenshtein/max length on the numbered Fv; **< 0.70 = de_novo class**), `novelty_class_local`, `cdrh3_min_dist_ref` / `cdrh3_nearest_ref(_src)` (min distance to abguard's 610k known CDRH3s), `cdrh3_novel_local` (>= 2 mutations de novo, >= 3 if class `novel`) |
| known-antibody panel | `identity_nearest_known` (whole Fv, 183,567 known heavy domains), `nearest_known_id`, `nearest_known_source`, `identity_nearest_vhh`, `ident_lt70`, `ident_lt68` |
| real abguard (VM) | `abguard_pass` (`is_valid`, bundled reference, novelty gate on), `abguard_errors`, `abguard_class` (`de_novo`/`novel`), `abguard_identity_trast`, `abguard_cdrh3_min_dist`, `abguard_nearest`, `abguard_n_high_liab`, `abguard_liab_tags` (LAP liabilities) |
| verdicts | `basic_pass` (local hard rules + CDRH3 novelty), `denovo_vm_pass` (abguard_pass and class `de_novo`), `denovo_strict_pass` (`denovo_vm_pass` and `ident_lt70`), `gap_vm_ok_ident_ge70` (passes the VM tool but sits >= 70 % to something known: **the gap to watch**) |

`identity_nearest_known` is `1.0` for a sequence immunum cannot number (worst case, with `numbering_ok = False`).

## diversify.py: two tiers per parent

Parents come from `score_of3.csv` / `score_orth.csv` (top 40 by `consensus_rank_score`, falling back to `of3_pct`/`of3_iptm`, plus the best 3 of
each generator), restricted to designs that pass every hard rule except identity and whose nearest-known identity is still >= 0.68. From
14:40 an unscored fallback (12 parents max, only while fewer than 12 scored parents exist) keeps the downstream stages busy. Cap: 250 variants.

| tier | ids | what | goal |
|---|---|---|---|
| `lite` | `<parent>-fl1`, `-fl2` | the fewest natural substitutions (<= 8) that bring whole-Fv identity to trastuzumab VH <= 0.66 | clears the abguard `de_novo` class with margin; keeps everything else |
| `deep` | `<parent>-fd1..3` | greedy descent to nearest-known identity <= 0.68 / 0.66 / 0.64 (one target per variant) | the "under 70 % to ANY known antibody" hedge, verifiable only against our panel |
| `cr` | `<parent>-cr` | restore the canonical Cys23 / Cys104 (and Trp118) in designs that lack them (BoltzGen never writes Cys) | rescue; CDRs untouched |

Common rules: substitution only (no indels, geometry of the parent stays valid), CDR1-3 untouched (checked after re-numbering),
Cys23/Cys104, Trp41, the FR2 hallmarks (IMGT 42/49/50/52), W118-G119-G121 and the CDR-flanking positions 26/39/55/66 protected; buried positions
(relative SASA < 0.20 in 3EAK) take only conservative hydrophobic swaps; a residue is allowed at a position only if >= 1.5 % of the panel
carries it there (natural residues only); never Cys, new Met/Trp or an N-glycosylation sequon. Choice among the best 8 candidates per step is sampled
(softmax) with a per-variant random 25 % position mask, a ProteinMPNN vote (96 samples on the 3EAK backbone, hosted on Nebius), and a bonus for residues of a
different human IGHV germline, so the variants of one parent differ pairwise by >= 8 positions (lite: >= 3). Everything is seeded (sha1 of design id, tier, retry).
`fwdiv_report.csv` records substitutions (IMGT label), identity to parent, trastuzumab / nearest-known identity before and after, tier and target.
`fwdiv_eval.py` joins S3's scores of variants and parents and tabulates which tier / substitution count / target keeps the parent's OpenFold3 signal.

Optional ESM-2 650M check (`esm_filter.py submit|collect`, Nebius, masked marginals of the substitution set in the parent context): informational only,
see `decisions.md`.
