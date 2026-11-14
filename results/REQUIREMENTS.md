# Entry requirements, per design

`SUBMISSION.csv`, 10 rows, checked against all 525063 provided CDR-H3s.

PASS = checked and holds. FAIL = checked and does not. UNKNOWN = only `abguard` can answer; never reported as a pass.

| Design | format | length_max_250 | no_unpaired_cysteines | cdrh3_identity_under_70pct | mutations_from_provided_data | immunum_numbering_confidence | identity_to_known_antibodies |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `tz-opt-02-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-18-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-07-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-15-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-01-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-03-stab` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-02` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-18` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-07` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |
| `tz-opt-01` | PASS | PASS | PASS | PASS | PASS | UNKNOWN | UNKNOWN |

## Margins

Every design sits at the minimum distance the rule allows. The margin column is zero by construction, because the generator was built to stay as close to trastuzumab as the rules permit.

| Design | CDR-H3 | Distance to nearest provided | Margin |
| --- | --- | ---: | ---: |
| `tz-opt-02-stab` | MGPEPFYAMD | 3 | 0 mutation(s) of margin |
| `tz-opt-18-stab` | WGGWPAYALD | 3 | 0 mutation(s) of margin |
| `tz-opt-07-stab` | FGGWPWYAMD | 3 | 0 mutation(s) of margin |
| `tz-opt-15-stab` | FGGDMWYWMD | 3 | 0 mutation(s) of margin |
| `tz-opt-01-stab` | APGDGFQIMD | 3 | 0 mutation(s) of margin |
| `tz-opt-03-stab` | LGGDGAYMVD | 3 | 0 mutation(s) of margin |
| `tz-opt-02` | MGPEPFYAMD | 3 | 0 mutation(s) of margin |
| `tz-opt-18` | WGGWPAYALD | 3 | 0 mutation(s) of margin |
| `tz-opt-07` | FGGWPWYAMD | 3 | 0 mutation(s) of margin |
| `tz-opt-01` | APGDGFQIMD | 3 | 0 mutation(s) of margin |
