# S4: yield, melting temperature and developability of the designed VHHs

> **STATUS 14 Nov 14:35: monoyield is PAUSED.**  User rule: nothing starts on the shared VM (GPU or otherwise) without the user's explicit
> permission.  `my_yield` / `my_mdratio` / `my_dtype` are therefore **partial**: 50 rows only (7 reference nanobodies, 14 controls, 29 RFD4 designs
> from the first candidate batch at 14:28-14:30); every later design has them blank.  Everything else (liab, nb9, NanoMelt on CPU, TNP) runs locally
> and continuously.  The VM code paths are disabled unless `--allow-vm` (props.py) / `S4_ALLOW_VM=1` (run_vm.sh) is given, which only the user may decide.

Reads every design in `results/denovo/designs_*.csv` (plus the reference nanobodies in `baselines.csv`), writes one row per
`design_id` to `results/denovo/props.csv`.  Incremental, resumable, safe to run while other stages are still producing designs.

| Question | Tool | Where it runs | Columns |
|---|---|---|---|
| Will it melt late? | **NanoMelt** (Ramon et al., CC BY-NC-SA 4.0, non-commercial) | local env `nanomelt`; mirrored on VM | `nm_tm` (apparent Tm, C) |
| Is it a developability outlier? | **TNP**, Therapeutic Nanobody Profiler (Gordon & Deane, Oxford) | local env `tnp` (GPU); mirrored on VM | `tnp_cdr_len tnp_cdr3_len tnp_cdr3_compact tnp_psh tnp_ppc tnp_pnc tnp_tetrad tnp_flag_* tnp_n_red tnp_n_amber` |
| Will it express as a monomer? | **monoyield-predictor** ESMC-6B nanobody models (Tsuboyama lab) | VM only (needs 26 GB fp32 / 13 GB bf16); **paused, partial** | `my_yield(_sd)` log10 mg/L, `my_mdratio(_sd)` log10 monomer/dimer, `my_dtype` |
| Cheap liabilities | `liab.py` (pure Python + ANARCI numbering) | local; mirrored on VM | `liab_*` counts, `liab_f_*` flags, `liab_n_flags`, `liab_flags` |
| Yield prior from CDR2/CDR3 residues | `nanobody_rules.py`: the nine-rule (and six-rule) score of the monoyield paper, Fig. 5h | local; mirrored on VM | `nb9_score` (0-9), `nb9_applicable`, `nb6_score`, `nb6_applicable`, `nb9_flags`, `nb_has_cys_met`, `nb_cdr2_core`, `nb_cdr3_loop` |

Housekeeping columns: `design_id, is_control, generator, framework, vhh_seq, stage_ts, tools_done` (`tools_done` lists the tools whose value really exists for that row, e.g. `liab;nb9;nanomelt`; count a row as finished for a tool only if it is listed).  `is_control = 1` for `ctl-*` rows and the
reference nanobodies.  Blank = not computed (yet).  Rows are keyed by `(design_id, vhh_seq)`: a changed sequence under the same id is
recomputed from scratch.

## Run

```bash
# one pass now, local only (the default; monoyield and every VM path are off unless --allow-vm, which needs the user's permission)
python3 denovo/s4_props/props.py --tools liab nb9 nanomelt tnp --where local
# continuous (what the pipeline uses): every <= 5 min, stops 15:45 local or when denovo/STOP exists
nohup python3 denovo/s4_props/props.py --watch 240 --tools liab nb9 nanomelt tnp --where local >> denovo/s4_props/watch.log 2>&1 &   # pid in denovo/s4_props/watch.pid
touch denovo/STOP                                  # graceful stop (after the current pass)
# everything on the VM from a CSV (design_id,vhh_seq), outputs to a local folder -- ONLY with the user's permission
S4_ALLOW_VM=1 bash denovo/s4_props/run_vm.sh in.csv out_dir [nanomelt,tnp,liab,nb9,monoyield] [auto|fp32|bf16]
# ad-hoc sequences, separate output file
python3 denovo/s4_props/props.py --in-csv x.csv --no-designs --out /tmp/p.csv --where local --tools liab nanomelt tnp
```

`--where auto|local|vm|both`: `both` also computes liab/NanoMelt/TNP on the VM into `*_vm` columns (cross-machine check).
`--max-tnp N` caps TNP work per pass.  `--tier off` runs every tool on every design.

### Cost tiers (why not everything on everything)
* liab + nb9 + NanoMelt: every design (NanoMelt about 0.6 s per sequence on 8 CPU threads, in batches of 32).
* TNP + monoyield: controls, reference baselines, and designs that pass S2's gate (`denovo_vm_pass` in `rules.csv`) AND are in the top
  35% of those gate-passers by `of3_consensus_rank_score` (S3, `score_of3.csv`, higher = better).  Rows with generator `fwdiv`,
  `boltzgen`, `boltzgen-cys` are ignored (cancelled); their old rows may remain in `props.csv`.  A monoyield call costs about 8-10 min whatever the batch size
  (ESMC-6B is reloaded for each of the 10 fold checkpoints), so there is one batched VM call per cycle, at most one per 10 min.

## How the pieces work

* `tnp_batch.py`: one process, NanoBodyBuilder2 on the GPU (about 2-3 s per VHH) then TNP's own metric code from `tools/TNP`
  imported unchanged.  The stock `TNP --file --ncores` takes about 4 min per sequence (CPU, one process per sequence).
* `liab.py`: IMGT CDRs by ANARCI; N-glyc (N-x-S/T, x not P), NG/NS deamidation, DG isomerisation, DP cleavage, Met/Trp in CDRs, free /
  unpaired Cys (odd count or canonical IMGT 23/104 missing), pI and net charge at pH 7.4 (Biopython), Kyte-Doolittle CDR GRAVY, max
  5-mer window, longest V/I/L/M/F/W run, K+R and D+E in CDRs.  Flags: unpaired_cys, nglyc, ng_cdr, dg_cdr, met_trp_cdr (M+W >= 3),
  hydrophobic (5-mer KD >= 2.5 or run >= 5), pI (< 5.5 or > 9.5), cdr_pos (K+R >= 6).
* `vm_install.sh push|install|status`: builds the VM env `s4` (micromamba under `~/moritz/s4`, nothing else touched).
  `remote_s4.sh` is what runs on the VM; `run_vm.sh` is the local driver (rsync, scp, ssh, scp back).
* `nanobody_rules.py`: CDR2 = chars [1:7] of an 8-residue IMGT CDR2 (positions 55-60; other lengths: CDR2 rules not applicable); CDR3 loop = IMGT
  105-117 minus the first two residues (after framework C-A-A, before W-G-x-G); rules 1-4 CDR2 pos 1/4/5/6, rule 5 C1 in Y/W, rule 6 C3 in F/H/L/Y,
  rule 7 N1 not in Y/W/V/F/L/I, rule 8 N2 in L/T, rule 9 C2 not in F/L/Y; six-rule = rules 1-4, C1 = Y, C3.  `--selftest` checks the parent
  hNbBCII10-FGLA = 6/9 (fails N1=V, N2=R, CDR2 pos6=L).  Caveat from the paper: the added rules were derived on the same dataset (33% of nine-rule
  designs > ~3 mg/L vs 29% for the six-rule scheme), one fixed framework, CDR3 of 5-13 residues; ours are longer.  A small prior, not a gate.
* `baselines.csv`: reference nanobodies with the PDB entry each sequence was verified against in RCSB (3EAK, 3DWT, 7EOW, 3K1K, 3P0G, 4KRL)
  and the trastuzumab VH as a labelled non-nanobody sanity row.

## Install state (14 Nov 14:00-14:10 CET)

| Tool | Local | VM (`~/moritz/s4`) |
|---|---|---|
| NanoMelt | env `nanomelt`, proven, run on the CPU only (user decision; GPU ran out of memory) | env `s4`, installed and proven (14:05), Tm identical to local (0.00 C on 7 sequences); not in use now |
| TNP | env `tnp`, proven | env `s4`, proven (patch metrics model-dependent, see decisions D1) |
| liab | env `tnp` | env `s4`, outputs identical |
| monoyield | not possible (8 GB GPU) | `~/moritz/monoyield-predictor` venv, proven (fp32 and bf16, max difference 0.010 log10 mg/L on 7 sequences); **paused, 50 rows done** |

## Known limits
See `decisions.md` D4.  In short: NanoMelt MAE about 4-7 C and out of distribution for de novo CDR3s; trastuzumab (Fab/IgG) versus VHH Tm is
not like for like; TNP flags are not Tm predictors and its patch metrics move by +-15% between two NanoBodyBuilder2 models of the same
sequence; monoyield is an E. coli monomer-yield proxy (Adaptyv expresses cell-free); no SoluProt/NetSolP/Aggrescan3D (not attempted).
