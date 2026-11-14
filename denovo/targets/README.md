# Targets
* `her2_domIV.cif`: HER2 domain IV as ONE chain `A`, 109 modelled residues (489-580, 591-607 of 1N8Z; 581-590 unresolved).
  `label_seq_id` is contiguous 1..109; `auth_seq_id` keeps HER2 numbering. Map in `numbering_map.json`.
* `her2_domIV.fasta`: the 109 modelled residues. (The full 119 aa domain IV sequence for folding is in `her2_domIV.a3m` query.)
* `her2_domIV.a3m`: existing MSA for domain IV (119 aa, 489-607) used by the earlier refolds.
* `ref_1n8z.cif`: the crystal complex, chains A=VL, B=VH, C+D=HER2.
* `epitope.json`: HER2 residues within 4.5 A of trastuzumab Fv heavy atoms (19 residues, from the crystal, not
  from literature lists) and three hotspot sets: hs1_core (6 most-contacted of the dominant patch), hs2_broad (12 most-contacted),
  hs3_alt_patch (6 of the other patch).
Regenerate with `baseline/.venv/bin/python denovo/targets/make_targets.py`.
