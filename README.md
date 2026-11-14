# De novo VHH design against HER2

How I would run a de novo campaign on the Nebius Scientific AI Platform, as
runnable code rather than a diagram. Standard library only, no install step.

```bash
# see the whole plan without spending a GPU-second
python3 -m denovo --forge ~/.claude/skills/nebius-scientific-ai/forge.py --dry-run

# run it
export NEBIUS_CREDS=/path/to/credentials.private.json
python3 -m denovo --forge ~/.claude/skills/nebius-scientific-ai/forge.py
```

## The argument

Generating antibody sequences is the easy part — AntiFold returns twelve in
110 seconds. Three things make the difference between a list of sequences and
a panel someone would act on, and all three are cheap:

**1. A control, scored the same way.** Trastuzumab is co-folded against the
same target, through the same model, before any candidate is. A candidate
scoring 0.74 means nothing on its own; 0.74 against a known binder's 0.77,
with 60% of its contact count, means something. Contact counts in the report
are normalised against it, so a column reads as "fraction of a known binder".

**2. The interface is measured, not reported.** Boltz2's output has no ipTM
(the platform's own recipes say so explicitly). Its confidence and pTM
describe the *model*, not the *interface* — a confidently-folded monomer
sitting next to a target scores well on both. So the interface is computed
from the returned coordinates: heavy-atom contacts across chains within 4.5 Å,
interface residue count, and the overlap with the trastuzumab epitope. That
epitope is itself measured, from the 1N8Z crystal structure, rather than
quoted from memory as a list of residue numbers.

**3. Triage before the GPU, not after.** Sixteen teams share the platform, so
the scarce resource is wall-clock. A sequence carrying a free cysteine or a
new N-glycosylation sequon in CDR3 is not a candidate a drug discovery team
would carry forward, whatever its predicted interface looks like. Those are
found locally in milliseconds, so they never consume a 75-second co-fold slot.
The rejected list is kept and reported — a candidate list without its
rejections cannot be audited.

## The pipeline

| stage | model | why this one | first call |
|---|---|---|---|
| `target` | — | 1N8Z gives the target sequence, a positive control, and the epitope in one download | — |
| `msa` | `msa-search-pdb70` | one alignment, reused by every later co-fold | 40 s |
| `baseline` | `boltz2` | target alone, then trastuzumab + target: the yardstick | 80 s |
| `generate` | `antifold` | antibody-aware CDR redesign: keeps the framework, proposes new loops | 110 s |
| `triage` | — | local liability scan; drops what no team would progress | instant |
| `build` | `nanobodybuilder2` | antibody-specific model, far better than a general folder for a VHH | 130 s |
| `cofold` | `boltz2` | candidate against target, the actual question | 80 s each |
| `score` | — | contacts, interface size, epitope overlap from the coordinates | instant |
| `confirm` | `openfold3` | second opinion; two models agreeing beats one high number | 165 s |
| `report` | — | ranked CSV + markdown, run id on every row | instant |

### Why AntiFold is the default generator

Three routes are implemented; `--route` picks one.

- `antifold` **(default)** — inverse folding with an antibody prior. It keeps a
  VHH framework and redesigns CDR1/2/3, so what comes out is still a nanobody.
  110 s.
- `rfd4` — `rfd4-proteina` for a binder backbone guided to a site, then
  `proteinmpnn` for a sequence. No antibody prior: the output is a general
  protein binder, not a VHH, and it needs the sequence-design step after it.
  ~5 min + 100 s, and it needs an *annotated* mmCIF marking the site.
- `mosaic` — de novo binder from the target sequence alone. The most literally
  de novo of the three, and the slowest at ~22 min. It has no tested recipe,
  so its input entry is read from the model's own
  `input_artifact_contract` rather than guessed at.

Default is `antifold` because an antibody prior is worth more here than
generality: the deliverable is an antibody, and 110 seconds against 22 minutes
buys more iterations in an afternoon. `mosaic` is worth starting early in
parallel if you want a non-antibody comparison — queue it first, collect later.

The framework donor is 1ZVY chain A, a VHH. Only its framework is used; the
CDRs are replaced, so the scaffold's original target (lysozyme) is irrelevant.

## Operational choices that cost people time

These are enforced in `forge_client.py` rather than left to whoever runs it:

- **Nothing is submitted twice.** Every call carries a label. A label with an
  existing run id in forge.py's receipts is *followed*, not resubmitted.
  Resubmitting does not get an answer sooner; it means waiting twice.
- **A storage failure resumes.** `private storage is temporarily unavailable`
  is a platform hiccup: the code calls `resume` on the same run id. It never
  resubmits on that error.
- **Long stages queue in parallel.** Co-folds and builds are submitted with
  `--no-wait`, then collected, so six candidates cost one wait rather than six.
- **The wait is announced before it starts**, per model, from the measured
  first-call times.
- **OpenFold2 is never called.** It failed to get an endpoint up on 12
  November, under load and idle, burning 55 minutes per attempt. It is in
  `config.AVOID`.
- **Everything is resumable.** Run ids live in forge.py's `nebius-runs/`
  receipts and campaign state in `results/campaign.json`, so you can close the
  laptop and pick up with `--stages score,confirm,report`.

## What this does not do

- **It does not predict affinity, and nothing here ranks candidates by
  binding strength.** A co-fold is not a binding assay. Two models agreeing
  that a complex is plausible is weak evidence of binding and no evidence of
  how tightly. The composite score is a documented triage ordering, with its
  components printed beside it so you can disagree with the weighting.
- **CDR1 and CDR2 are only labelled when AntiFold returns its region CSV.**
  Otherwise the liability scan treats them as framework and under-counts. Real
  IMGT numbering for arbitrary sequences needs `anarcii`.
- **The control calibrates the scoring, not the design.** Trastuzumab is an
  IgG Fab; the candidates are single-domain VHHs. A smaller paratope makes a
  smaller interface before any question of quality enters.
- **Low epitope overlap is not failure.** It may be a different and useful
  epitope. The column is information, not a verdict.
- **Not for clinical use.**

## Files

| file | role |
|---|---|
| `config.py` | model ids, lanes and timings copied from forge.py's catalog; campaign sizes |
| `forge_client.py` | idempotent wrapper around forge.py: labels, receipts, resume, contract reading |
| `bio.py` | PDB/mmCIF parsing, FASTA, tar, interchain contacts, epitope overlap |
| `triage.py` | CDR3 locator, liability motifs weighted by region, pI / charge / GRAVY |
| `pipeline.py` | the stages and the campaign state |
| `report.py` | ranked CSV and the markdown note, including the rejections |
| `__main__.py` | CLI |
| `tests/test_denovo.py` | 24 tests: geometry against hand-placed coordinates, CDR3 against a real VHH, offline CLI dry-run |

## Tests

```bash
python3 -m unittest discover -s tests
```

No network and no credentials needed. The contact counter is tested against
atoms placed by hand at known distances, including the 4.5 Å boundary, so a
wrong count fails a test instead of quietly producing a wrong ranking.

One finding worth recording: the CDR3 locator first used a regex anchor on the
framework-3 `[YF]C` motif, which **misses real VHHs** — the nanobody in the
platform's own `RECIPES.md` ends framework 3 with `...TAIYTC`. It now anchors
on the conserved cysteine itself, scanning for the first one that leaves a
plausible loop length, which also makes an extra cysteine inside CDR3 unable
to truncate the loop. Both cases are tests.

## Keys

`forge.py` finds credentials via `--creds`, then `$NEBIUS_CREDS`. This package
passes the path through and never reads, prints or copies the file. The
submission repo is public, so before the first commit:

```bash
echo "*.private.json" >> .gitignore
```

`nebius-runs/` and `results/` hold run ids and outputs, no keys. They are
worth committing: they are how a score can be traced back to the run that
produced it.
