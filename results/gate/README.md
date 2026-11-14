# Entry gate: all 20 designs checked against the organisers' official checker

Run on 14 November 2026. This is the organisers' own validator,
not our local reimplementation in `src/abx/gate.py`.

```
official-checker team-validate gate/designs_only_paired.csv
official-checker novelty-report  gate/designs_only_paired.csv -o gate/novelty_report.json
```

Using the reference bundled with the official checker (`reference_data`).

## Result

**All 20 designs pass, in both submission formats.**

| Table | Rows | Pass |
| --- | --- | --- |
| `submission_paired.csv` (20 designs + trastuzumab control) | 21 | 20 |
| `submission_scfv.csv` (20 designs + trastuzumab control) | 21 | 20 |
| `designs_only_paired.csv` (the 20 designs) | 20 | **20** |

The single failure in the first two tables is trastuzumab itself, included
deliberately as a negative control:

```
[novelty] tz_reference: CDRH3 SRWGGDGFYAMDY is 0 mutation(s) from
          SRWGGDGFYAMDY in trastuzumab; a novel design needs at least 3
```

That is the gate working. A run where trastuzumab passes would mean the
novelty check is not switched on, and the result would be worthless.

## What the run taught us that the task description did not

**1. The enforced rule is 3 mutations, not 4.**
`required_cdrh3_mutations: 3` for every design. The written rules are
ambiguous ("at least 2" in one place, "strictly more than 2" in another, plus
a separate "under 70% CDRH3 identity"), and we had read the identity clause as
implying 4. The official checker enforces 3. Our designs carry 4, so they satisfy both the
enforced rule and the strictest reading — but the margin is ours, not required.

**2. CDRH3 is 13 residues here, not 10.**
The official checker numbers trastuzumab's CDRH3 as `SRWGGDGFYAMDY` (13), while
`src/abx/` works on the 10-mer core `WGGDGFYAMD` that the Chinery screen
mutates. Both describe the same designs; only the denominator of any identity
fraction changes. On 13 residues, 4 mutations is 69.2% identity, just under
the 70% line. On 10 residues it would be 60%. Do not quote an identity
percentage without saying which length it is over.

**3. The bundled reference does not contain the provided HER2 data.**
For all 20 designs the nearest known CDRH3 is trastuzumab itself, at
distance 4:

```json
{"cdrh3_nearest_source": "trastuzumab", "cdrh3_min_distance": 4,
 "required_cdrh3_mutations": 3, "novelty_class": "novel"}
```

If the 524,346 provided variants were in the reference, the nearest neighbour
would be one of them. Our own gate puts the nearest provided sequence at
**distance 3** — passing, but with zero margin.

**This risk is now closed.** See the section below: we built that reference
and re-ran. All 20 designs pass.

`identity_to_trastuzumab` is 0.9824 for all 20 — 4 mutations across the
226-residue paired V region. That is global identity, and it is not what the
novelty rule is measured on; the rule is measured on CDRH3.

---

# Second run: checked against the provided HER2 data

The rule that actually constrains us is distance from the 525,063 provided
CDRH3s, and the bundled reference does not contain them. So we built one.

```
python scripts/build_reference.py --out provided_ref.csv
official-checker reference-ingest provided.duckdb provided_ref.csv
official-checker reference-build  provided.duckdb
official-checker team-validate <table> --reference-db provided.duckdb
```

`scripts/build_reference.py` grafts each provided 10-residue loop into
trastuzumab's VH, because `reference-ingest` reads the CDR3 off a heavy chain.
It keeps every provided sequence within 4 mutations of any design, plus 400
random controls: 4,579 sequences. Sequences further than 4 away cannot change
any design's minimum distance, so this is a faithful test of the constraint,
not a shortcut around it.

## Every design sits at exactly 3

```
tz-opt-01 ... tz-opt-20: nearest provided sequence is 3 mutation(s) away
```

All twenty. The requirement is 3. **Zero margin, uniformly.** That is not a
coincidence -- the generator was built to sit as close to trastuzumab as the
rules allow, and this is what that looks like from the other side.

## Result: all 20 pass

With 20 designs and 3 planted controls in one table, `--drop` removed exactly
the three controls and kept all twenty designs:

```
[paired] heavy_cysteines_are_paired: 1 row(s)
[novelty] control_d1_0_0: CDRH3 SRCDTTACTLWIY is 1 mutation(s) from
          SRADTTACTLWIY in provided_ref; a novel design needs at least 3
[novelty] control_d1_2_0: CDRH3 SRAAAGRFYTKDY is 1 mutation(s) from
          SRCAAGRFYTKDY in provided_ref; a novel design needs at least 3
```

## Why the controls matter

A gate that passes everything looks exactly like a gate that is switched off.
`scripts/reference_control.py` builds designs one mutation from a sequence
that IS in the reference and far from trastuzumab, so the only thing that can
reject them is the provided-data comparison. They are rejected, naming the
provided sequence they collide with. The filter is demonstrably on, and the
twenty passed the same run that caught them.

## Gotcha: `novelty-report` does not report the reference comparison

This nearly produced a false alarm. With `--reference-db` pointing at the
provided data, `novelty-report` still prints:

```json
{"cdrh3_min_distance": 4, "cdrh3_nearest_source": "trastuzumab"}
```

and it prints the same thing for a control sitting 1 mutation from a reference
sequence (distance 10 from trastuzumab, reported as `nearest: trastuzumab`).

`cdrh3_min_distance` and `cdrh3_nearest_source` are the **trastuzumab**
comparison. They do not reflect the reference database. Only `team-validate`
enforces the provided-data check, and only its error lines name the colliding
sequence. Do not read a novelty-report as evidence about prior art.

## Files

| File | What it is |
| --- | --- |
| `submission_paired.csv` | 20 designs as VH/VL pairs, plus the trastuzumab control |
| `submission_scfv.csv` | the same designs as single-chain scFv, plus the control |
| `designs_only_paired.csv` | the 20 designs alone, the table that passes cleanly |
| `novelty_report.json` | per-design novelty fields (trastuzumab comparison only) |
| `novelty_vs_provided.json` | the same, run with `--reference-db`: identical, which is the gotcha above |
| `provided_ref.csv` | the 4,579 provided sequences the reference was built from |
| `control.csv` | planted designs that must fail, and do |

Regenerate with `scripts/build_submission.py`, `scripts/build_reference.py`
and `scripts/reference_control.py`.

---

# Third run: the final ten, with exit codes

The ten rows actually submitted (`results/SUBMISSION.csv`), rebuilt from
`results/designs.csv` and validated against the provided-data reference.
Full output in `submission_validation.txt`.

| Run | Exit code | Meaning |
| --- | ---: | --- |
| the ten alone | **0** | every row compliant |
| the ten plus trastuzumab | 1 | 10 of 11 pass; only the control fails, which is correct |

```
>>> EXIT CODE: 0  (0 = every row compliant, 1 = at least one failure)

ERROR  SUB_with_tz.csv: submission is not compliant (10 of 11 row(s) pass)
  - [novelty] tz_reference: CDRH3 SRWGGDGFYAMDY is 0 mutation(s) from
              SRWGGDGFYAMDY in trastuzumab; a novel design needs at least 3
```

This closes the gap the status page flagged: earlier runs validated the
twenty parents, not the ten submitted rows, and the output was never saved.
`tz-opt-15-stab` and `tz-opt-03-stab` are covered here — they were the two
with no evidence file of their own.

**Reproducible.** `scripts/build_final.py` run from a clean copy of
`results/designs.csv` produced a file byte-identical to the committed
`results/SUBMISSION.csv`, with the same six loops in the same order.

## The reference was rebuilt on edit distance, and the verdict held

Both this reference and our exhaustive check originally measured equal-length
loops by Hamming. That is the wrong metric for this rule. Two 10-mers can
differ at three positions and still be one insertion plus one deletion apart —
two edits — and the official checker aligns, so it reads the 2. A design we scored at 3
can therefore be rejected at 2, and a rejected design is never made.

The case is real, not hypothetical: a published design from another entry sits
at Hamming 3 and edit distance 2 from a provided sequence, and the official checker
rejects it. `scripts/verify_novelty_full.py` and `scripts/build_reference.py`
now both use edit distance, pinned by
`tests/test_verify_novelty_full.py::test_equal_length_neighbour_reached_by_a_shift_is_an_edit_distance`.

Rebuilt and re-run end to end:

| Artefact | Result |
| --- | --- |
| `provided_ref_edit.csv` | 6,608 sequences (4,621 within 4 edits of a design, 1,987 controls) |
| exhaustive check, all 525,063 loops | every submitted design still at **exactly 3** — one design's nearest neighbour changed, the verdict did not |
| `official-checker team-validate` on the twenty + control | 20 of 21 pass; only trastuzumab fails |
| `official-checker team-validate` on the ten + control | 10 of 11 pass; only trastuzumab fails |

`validation_edit.json` and `validation_edit10.json` hold the two runs. The
margin is still zero — three mutations against a requirement of three — but it
is now zero under the metric the authority actually uses.

