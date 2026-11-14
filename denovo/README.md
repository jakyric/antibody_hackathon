# De novo VHH pipeline against HER2 domain IV

A continuously running, **filter-first** pipeline that designs new single-domain antibodies (VHHs) against the trastuzumab epitope on HER2 domain IV, and
chooses a final list of 20 by Pareto dominance over binding evidence, epitope proxy, predicted Tm, predicted yield and developability. Built for the
Oslo AI x Antibodies hackathon (14 Nov 2026), where designs are expressed and measured by Adaptyv Bio against trastuzumab (SPR affinity, epitope competition,
thermostability). Every stage is resumable, writes partial results as it goes, and is safe to stop at any moment.

**Scores are filters, not binding evidence.** Nothing in this repository has been measured in the lab (see Limitations).

```mermaid
flowchart LR
  T["S0 targets<br/>HER2 domain IV + 19-residue epitope<br/>from the 1N8Z crystal"] --> G["S1 design<br/>RFD4-Proteina on JapanFold<br/>intact framework, CDR1-3 designed"]
  G -->|designs_rfd4jf.csv| R["S2 hard rules + novelty gate<br/>abguard class, canonical Cys, nine rules"]
  R -->|gate passers| S["S3 scoring<br/>OpenFold3, then OpenDDE, ESMFold-2, AF2-IG"]
  S --> P["S3d Protenix-v2 (Nebius)<br/>top designs only"]
  G --> D["S4 properties<br/>NanoMelt Tm, TNP, monoyield, liabilities"]
  S --> D
  R --> L["S5 leaderboard (rank.py)<br/>weighted composite"]
  S --> L
  P --> L
  D --> L
  L --> A["S6 Pareto (pareto.py)<br/>final top 20"]
  S -. feedback.py .-> G
```

## Quickstart

```bash
export JAPANFOLD_API_KEY=...          # required for design and scoring; never commit or print it
bash denovo/setup.sh                  # venv, key/target checks, no-cost smoke test   (--check = read-only, --tools, --panel optional)
bash denovo/run_pipeline.sh dry-run   # the resolved plan and the exact commands; starts nothing
bash denovo/run_pipeline.sh start     # background loops, logs in denovo/.run/
bash denovo/run_pipeline.sh status    # alive loops, funnel, spend
bash denovo/run_pipeline.sh final     # stop, drain, rank + Pareto, snapshots into results_denovo/, safety scan
```

Settings are environment variables with safe defaults (deadline `PIPELINE_UNTIL`, spend caps, tier-2 scorers, VM opt-in); see the header of
`run_pipeline.sh`. With Claude Code, the skill `.claude/skills/denovo-vhh-pipeline` runs this for you and explains the outputs.

## What each stage does

| Stage | Folder | What | Key output (in `results/denovo/`) |
|---|---|---|---|
| S0 targets | `targets/` | HER2 domain IV as one chain, 19 epitope residues within 4.5 A of trastuzumab, three hotspot sets, MSA | `targets/*` (committed) |
| S1 design | `s1_design/` | JapanFold `rfd4` (RFD4-Proteina, `rfd4-binder`) on two fixed VHH frameworks (hNbBCII10-FGLA, caplacizumab-like), CDR1-3 designed, hotspots on the epitope, two poses; adaptive waves steered by S3 | `designs_rfd4jf.csv` |
| S2 rules | `s2_rules/` | length, Cys pairing, canonical Cys23/Cys104, numbering, abguard track class (identity to trastuzumab VH < 70%), local replica and nearest-known identity; `gate_drivers.md` explains the gate | `rules.csv` |
| S3 scoring | `s3_score/` | OpenFold3 on JapanFold for every gate passer; OpenDDE-abag, ESMFold-2, second seed, AF2-IG for the top third; geometry of the pose vs trastuzumab; calibration against decoys and cognate/shuffled VHH pairs | `score_of3.csv`, `score_orth.csv` |
| S3d Protenix | `s3d_protenix/` | Protenix-v2 on Nebius (no MSA) for the top designs | `score_ptx.csv` |
| S4 properties | `s4_props/` | NanoMelt Tm (CPU), TNP flags, monoyield (optional, shared-VM GPU), sequence liabilities, nine CDR2/CDR3 rules | `props.csv` |
| S5 leaderboard | `rank.py` | weighted composite with explicit priors; gates; `rank_reason` per design | `leaderboard.csv` |
| S6 Pareto | `pareto.py` | non-dominated sorting over bind / epitope / Tm / yield / developability, weakest-objective tie-break, CDR3 and framework diversity | `leaderboard_pareto.csv`, **`submission_top20_pareto.csv`** |
| overview | `slides/` | rebuilds a 13-slide overview deck from the result files | |

The contract between stages (file names, columns, budgets, timeline, decisions) is in [`CONTRACT.md`](CONTRACT.md). Each stage folder has a README and a
`decisions.md` with the evidence behind its choices, including what was rejected.

## Design decisions in one paragraph each

* **Intact framework, new CDR1-3.** The organisers' novelty rule is met by fully new CDRs; no framework mutation is run (a diversification fallback is kept in
  `s2_rules/diversify.py`, not used). The VM abguard gate depends mostly on total length: see `s2_rules/gate_drivers.md`, which is why generation steers to longer CDR3s.
* **Gate first, then score.** A design that fails the hard rules or the abguard class cannot be submitted, so it is never scored (this halves the scoring cost).
* **Independent evidence, calibrated on decoys.** OpenFold3 plus orthogonal predictors, with gates set at the 95th percentile of decoys (shuffled VHH-antigen pairs,
  lysozyme and framework-only decoys), on metrics fixed before looking at the designs.
* **Pareto, not a weighted sum.** A design is tiered by how many of the judged properties beat trastuzumab; the weakest objective is what loses the tier.
* **Canonical disulfide is a hard gate.** Designs without Cys23/Cys104 lose Tm (BoltzGen designs, which exclude Cys, had a median predicted Tm of 60.9 C vs 68.2 C once restored).

## Tried and dropped

* **BoltzGen `nanobody-anything` (JapanFold).** The hosted protocol designs with Cys excluded, so every design lacked the canonical disulfide; stopped to focus on RFD4.
* **RFD4-Proteina on a shared VM (Apheris Foundry image).** Ruled out by the team (shared GPU).
* **Framework diversification (S2).** Built and tested, then cancelled: the organisers said fully new CDRs land under 70%, and the gate analysis shows the measured effect is length.
* **Templated OpenFold3 on the VM.** Stopped before it ran; Nebius and JapanFold OpenFold3 accept no template (only Chai-1 on Nebius does).

## Results snapshot

`results_denovo/` (created by `run_pipeline.sh final` or the packaging sync) holds the leaderboard, the Pareto list, the funnel and the per-design tables of the last run,
with a README giving the timestamp and row counts. No structures or job ids are stored.

## Limitations (read before trusting a number)

* **No design has been measured.** OpenFold3 confidence scores rank real VHH-antigen pairs poorly (Smorodina et al. 2026; reproduced on our calibration set: AUC about 0.5-0.65 for
  cognate vs shuffled pairs). On the designs, OpenFold3 sees almost no confident interface (median ipTM about 0.1, decoys 0.14), partly because our recipe is harder than upstream's
  (MSA on, no template; the upstream RFD4-Proteina oracle is templated and single-sequence, which JapanFold and Nebius OpenFold3 cannot do).
* **Every design shares one rigid binding pose** (the framework is placed where trastuzumab's VH sits), so epitope overlap is built in; the pose-based proxies are weak.
* **The known-VHH panel is incomplete** (PDB, PLAbDab and UniProt germlines, no OAS or patents), so nearest-known identity can only be too low. The organisers' abguard build compares
  the track identity only against trastuzumab and checks prior art only on CDRH3.
* **Predictors are filters.** NanoMelt has a mean error of about 4-7 C and a non-commercial licence; TNP thresholds come from clinical nanobodies; monoyield was trained on one scaffold.
* **Developability and expression predictors the organisers recommend (SoluProt, NetSolP, Aggrescan3D) are not included.**
* **Expect an A- or B-tier outcome, not S:** trastuzumab's sub-nanomolar affinity is the hardest criterion for a de novo VHH to match.

## Reproducibility notes

Python 3.10+, `denovo/requirements.txt` (immunum pinned to 2.0.0 to match abguard), optional micromamba envs in `denovo/env/`. All paths are repo-relative. Access to JapanFold
(`JAPANFOLD_API_KEY`) and, for Protenix-v2, Nebius (`NEBIUS_CREDS` plus the `nebius-scientific-ai` skill from the hackathon kit) comes from the environment; the optional VM is
configured by `PIPELINE_VM_HOST` and `PIPELINE_VM_KEY` and is **off by default**. The repository is public: run `python3 denovo/scan_public.py .` before any commit.
