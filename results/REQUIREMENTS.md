# Entry requirements, per design

`SUBMISSION.csv`, 10 rows, checked against all 525063 provided CDR-H3s.

PASS = checked and holds. FAIL = checked and does not. UNKNOWN = only `abguard` can answer; never reported as a pass.

| Design | format | length_max_250 | no_unpaired_cysteines | cdrh3_identity_under_70pct | mutations_from_provided_data | immunum_numbering_confidence | identity_to_known_antibodies |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `tz-clean-02-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-01-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-03-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-07-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-06-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-04-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-02` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-01` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-03` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-clean-07` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |

## Margins

Every design sits at the minimum distance the rule allows. The margin column is zero by construction, because the generator was built to stay as close to trastuzumab as the rules permit.

| Design | CDR-H3 | Distance to nearest provided | Margin |
| --- | --- | ---: | ---: |
| `tz-clean-02-stab` | RGGEGFYANE | 3 | 0 mutation(s) of margin |
| `tz-clean-01-stab` | YGGDYFYGKD | 3 | 0 mutation(s) of margin |
| `tz-clean-03-stab` | FGGGGNYAKD | 3 | 0 mutation(s) of margin |
| `tz-clean-07-stab` | HGGDRFYKRD | 3 | 0 mutation(s) of margin |
| `tz-clean-06-stab` | YGGAGFHAVD | 3 | 0 mutation(s) of margin |
| `tz-clean-04-stab` | HTGDEFYAID | 3 | 0 mutation(s) of margin |
| `tz-clean-02` | RGGEGFYANE | 3 | 0 mutation(s) of margin |
| `tz-clean-01` | YGGDYFYGKD | 3 | 0 mutation(s) of margin |
| `tz-clean-03` | FGGGGNYAKD | 3 | 0 mutation(s) of margin |
| `tz-clean-07` | HGGDRFYKRD | 3 | 0 mutation(s) of margin |
