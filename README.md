<h1 align="center">fabsolutely binding</h1>

<p align="center"><em>De novo and optimised antibodies against HER2 domain IV (the trastuzumab epitope)</em><br>
AI × Antibodies Hackathon, Oslo, 14 November 2026</p>

<p align="center">
  <a href="docs/index.html"><b>Landing page</b></a> ·
  <a href="fabsolutely_binding.csv"><b>Submission CSV</b></a> ·
  <a href="#the-submission-30-sequences">Submission</a> ·
  <a href="#optimization-track">Optimization</a> ·
  <a href="#de-novo-track">De novo</a>
</p>

Everything on this page is **computational**. No score here predicts affinity, epitope preservation or melting temperature of a real
molecule; the measurements come after the event.

## The submission (30 sequences)

One file, both tracks, in the organisers' format (`sequence_id,kind,heavy_sequence,light_sequence,design_description,track`):
[`fabsolutely_binding.csv`](fabsolutely_binding.csv) with **15 optimization** and **15 de novo** rows. The letter after the prefix of
each id marks which team member's approach a row comes from (`opt-J..`, `dn-M..`, `dn-K..`).

Checked with the organisers' `abguard team-validate`: **30/30 rows valid, 0 hard errors**. All 15 de novo rows are
classed `de_novo` (identity to the trastuzumab VH 0.63–0.68); the 15 optimization rows are classed `novel`
(identity 0.97–0.98, at least 3 CDR-H3 mutations from every known sequence).

| id | track | kind | heavy | light | description |
|---|---|---|---:|---:|---|
| `opt-J01` | optimization | paired | 120 | 107 | CDRH3 variant WDYQGFYYYD with Vernier L46L>M |
| `opt-J02` | optimization | paired | 120 | 107 | CDRH3 variant WDYQGFYYYD with Vernier L49Y>F |
| `opt-J03` | optimization | paired | 120 | 107 | CDRH3 variant WSSQGFYGRD with Vernier H47W>F |
| `opt-J04` | optimization | paired | 120 | 107 | CDRH3 variant WGSSGFYAWY with Vernier H94R>K |
| `opt-J05` | optimization | paired | 120 | 107 | CDRH3 variant WDYQGFYYYD with Vernier L49Y>H |
| `opt-J06` | optimization | paired | 120 | 107 | CDRH3 variant WGYQGFYGDD with Vernier H93S>T |
| `opt-J07` | optimization | paired | 120 | 107 | CDRH3 variant WGDQGFYYYY with Vernier H94R>Q |
| `opt-J08` | optimization | paired | 120 | 107 | CDRH3 variant WGYQGFYGDD with Vernier L46L>F |
| `opt-J09` | optimization | paired | 120 | 107 | CDRH3 variant WDSQGFYAWR with Vernier H94R>Q |
| `opt-J10` | optimization | paired | 120 | 107 | CDRH3 variant WGSAGFYGWD with Vernier L49Y>H |
| `opt-J11` | optimization | paired | 120 | 107 | CDRH3 variant WASQGFYRWD with Vernier L46L>I |
| `opt-J12` | optimization | paired | 120 | 107 | CDRH3 variant WGSAGFYGWD with Vernier H47W>Y |
| `opt-J13` | optimization | paired | 120 | 107 | CDRH3 variant WASQGFYRWD with Vernier H93S>A |
| `opt-J14` | optimization | paired | 120 | 107 | CDRH3 variant WASQGFYRWD with Vernier L46L>V |
| `opt-J15` | optimization | paired | 120 | 107 | CDRH3 variant WDSYGFYGYD with Vernier L46L>M |
| `dn-M01` | denovo | vhh | 119 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzuma… |
| `dn-M02` | denovo | vhh | 123 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzuma… |
| `dn-M03` | denovo | vhh | 121 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzuma… |
| `dn-M04` | denovo | vhh | 125 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzuma… |
| `dn-M05` | denovo | vhh | 127 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact hNbBCII10-FGLA framework against HER2 domain IV (trastuzumab e… |
| `dn-M06` | denovo | vhh | 126 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact hNbBCII10-FGLA framework against HER2 domain IV (trastuzumab e… |
| `dn-M07` | denovo | vhh | 125 |  | De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzuma… |
| `dn-K01` | denovo | vhh | 120 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AMGY… |
| `dn-K02` | denovo | vhh | 127 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 ARYD… |
| `dn-K03` | denovo | vhh | 128 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AAAL… |
| `dn-K04` | denovo | vhh | 120 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AARK… |
| `dn-K05` | denovo | vhh | 118 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 NASR… |
| `dn-K06` | denovo | vhh | 123 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AARK… |
| `dn-K07` | denovo | vhh | 119 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 NAKR… |
| `dn-K08` | denovo | vhh | 123 |  | de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AAWG… |

<details>
<summary>The full CSV (copy-paste ready)</summary>

```csv
sequence_id,kind,heavy_sequence,light_sequence,design_description,track
opt-J01,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWDYQGFYYYDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKMLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WDYQGFYYYD with Vernier L46L>M,optimization
opt-J02,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWDYQGFYYYDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIFSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WDYQGFYYYD with Vernier L49Y>F,optimization
opt-J03,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEFVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWSSQGFYGRDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WSSQGFYGRD with Vernier H47W>F,optimization
opt-J04,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSKWGSSGFYAWYYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGSSGFYAWY with Vernier H94R>K,optimization
opt-J05,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWDYQGFYYYDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIHSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WDYQGFYYYD with Vernier L49Y>H,optimization
opt-J06,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCTRWGYQGFYGDDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGYQGFYGDD with Vernier H93S>T,optimization
opt-J07,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSQWGDQGFYYYYYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGDQGFYYYY with Vernier H94R>Q,optimization
opt-J08,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGYQGFYGDDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKFLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGYQGFYGDD with Vernier L46L>F,optimization
opt-J09,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSQWDSQGFYAWRYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WDSQGFYAWR with Vernier H94R>Q,optimization
opt-J10,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGSAGFYGWDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIHSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGSAGFYGWD with Vernier L49Y>H,optimization
opt-J11,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWASQGFYRWDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKILIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WASQGFYRWD with Vernier L46L>I,optimization
opt-J12,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEYVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGSAGFYGWDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WGSAGFYGWD with Vernier H47W>Y,optimization
opt-J13,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARWASQGFYRWDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WASQGFYRWD with Vernier H93S>A,optimization
opt-J14,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWASQGFYRWDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKVLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WASQGFYRWD with Vernier L46L>V,optimization
opt-J15,paired,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWDSYGFYGYDYWGQGTLVTVSS,DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKMLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK,CDRH3 variant WDSYGFYGYD with Vernier L46L>M,optimization
dn-M01,vhh,EVQLVESGGGLVQPGGSLRLSCAASGPSFSSYAMGWFRQAPGKGRELVAAISSGGSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAASTHSSYYYDYDYWGQGTQVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzumab epitope, hotspot set hs1_core); CDR3 AASTHSSYYYDYDY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M02,vhh,EVQLVESGGGLVQPGGSLRLSCAASTSYAFSDFNMGWFRQAPGKGRELVASIKTSGSTSYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAYRTDTWNPDTGYLYYWGQGTQVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzumab epitope, hotspot set hs3_alt_patch); CDR3 AAYRTDTWNPDTGYLYY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M03,vhh,EVQLVESGGGLVQPGGSLRLSCAASTSYTLSHMGWFRQAPGKGRELVASISSGGSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAATLYTNSSSTIVFDYWGQGTQVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzumab epitope, hotspot set hs2_broad); CDR3 AAATLYTNSSSTIVFDY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M04,vhh,EVQLVESGGGLVQPGGSLRLSCAASGPTFTDYDMGWFRQAPGKGRELVALIYSTGAGSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAIYWPARTTGSLSEFDYWGQGTQVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzumab epitope, hotspot set hs2_broad); CDR3 AAIYWPARTTGSLSEFDY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M05,vhh,QVQLVESGGGLVQPGGSLRLSCAASGGSLTSMATLGWFRQAPGQGLEAVAAIRVYQGHETTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAALTHPSDWRDPSGYDYWGQGTLVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact hNbBCII10-FGLA framework against HER2 domain IV (trastuzumab epitope, hotspot set hs1_core); CDR3 AAALTHPSDWRDPSGYDY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M06,vhh,QVQLVESGGGLVQPGGSLRLSCAASGVSSSTFSSMMMGWFRQAPGQGLEAVAAMRGSTGTTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAALRDVDWTDSSAYDYWGQGTLVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact hNbBCII10-FGLA framework against HER2 domain IV (trastuzumab epitope, hotspot set hs1_core); CDR3 AAALRDVDWTDSSAYDY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-M07,vhh,EVQLVESGGGLVQPGGSLRLSCAASGPNFTDYWMGWFRQAPGKGRELVAMIRDSSSGRTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAGELSPSGSSHYVVFNYWGQGTQVTVS,,"De novo VHH, RFD4-Proteina generated CDR1-3 on an intact caplacizumab-like framework against HER2 domain IV (trastuzumab epitope, hotspot set hs3_alt_patch); CDR3 AAGELSPSGSSHYVVFNY; passes hard gates, Pareto front 1 of binding, developability, yield and Tm. Computational evidence only",denovo
dn-K01,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRSFSSYVMGWFRQAPGKGRELVASIGGYGSSTNYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAMGYGGSREDGTVYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AMGYGGSREDGTVY",denovo
dn-K02,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRETFAAYGMGWFRQAPGKGRELVAAISGSGLNTNYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCARYDRSYSGGSAHSISAYDYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 ARYDRSYSGGSAHSISAYDY",denovo
dn-K03,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRLNFNSFVMGWFRQAPGKGRELVAAISSSGSYTWYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAALSSNFYNARLHTASKYDYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AAALSSNFYNARLHTASKYDY",denovo
dn-K04,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRTFTDYAMGWFRQAPGKGRELVAWISHDGSSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARKTGTSSTDFDYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AARKTGTSSTDFDY",denovo
dn-K05,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRTFTDYWMGWFRQAPGKGRELVASINSSGSYTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCNASRPHVYAANYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 NASRPHVYAANY",denovo
dn-K06,vhh,EVQLVESGGGLVQPGGSLRLSCAASGRTFTDYGMGWFRQAPGKGRELVAYISSSGSSTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAARKSYDYSSETSGFDYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AARKSYDYSSETSGFDY",denovo
dn-K07,vhh,EVQLVESGGGLVQPGGSLRLSCAASVSGYSFSDYYMGWFRQAPGKGRELVAAINNTGSDTYYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCNAKRVGSTTSYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 NAKRVGSTTSY",denovo
dn-K08,vhh,EVQLVESGGGLVQPGGSLRLSCAASGNSTIFNNYAMGWFRQAPGKGRELVAVISSTGSTSYPDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAAWGSSTSSSSSAFDYWGQGTQVTVS,,"de novo VHH, RFD4-Proteina on the caplacizumab-like framework against HER2 domain IV (trastuzumab epitope); CDR-H3 AAWGSSTSSSSSAFDY",denovo
```

</details>

## Optimization track

*CDR-H3 + Vernier-zone redesign of trastuzumab* (branch [`jakob`](../../tree/jakob), full method in its README). Start from the best
experimentally tested CDR-H3 (`WGSQGFYAYD`, SPR 0.023 nM against 0.114 nM for trastuzumab), change 3–4 of the six positions that do not contact
HER2 with a BLOSUM-45 recipe, add Vernier-zone swaps where a Vernier residue touches a changed position, fold every sequence with JapanFold
OpenDDE antibody-antigen and keep the designs that dock like the crystal.

| step | designs |
|---|---:|
| raw designs (3–4 BLOSUM-45 changes) | 23,760 |
| no Cys, liabilities ≤ WT, ≥ 3 from every known sequence | 5,247 |
| diverse shortlist × Vernier genotypes | 528 |
| pass the organisers' checker | 517 |
| folded (+ 4 controls) | 521 |
| docked (ipTM ≥ 0.9) and within 3 Å of the crystal epitope | 201 |
| keep the crystal contacts | 38 |
| salt-bridge rule | 24 |
| final list | 20 (15 submitted) |

Controls: the base loop folds at ipTM 0.947 with epitope RMSD 1.3 Å; trastuzumab itself folds at 0.687 and sits 29.8 Å away (it does not dock in
this model). One seed per complex, so the top of the list is within noise.

## De novo track

Two independent VHH approaches, converged into one list of 15 sequences (7 + 8):

| approach | branch | idea |
|---|---|---|
| **A**, `dn-M..` (7) | [`moritz`](../../tree/moritz) | RFD4-Proteina on two intact VHH frameworks (1,544 designs), hard gates, OpenFold3 and ESMFold-2 on the VM, OpenDDE on a cropped HER2, TNP / NanoMelt / monoyield profiling, Pareto ranking |
| **B**, `dn-K..` (8) | [`krishagen`](../../tree/krishagen) | RFD4-Proteina VHHs on the caplacizumab framework, selected by its own rule: predicted epitope overlap ≥ 0.2 and a refold within 10 Å of the designed pose |

Also explored, **not submitted**: CDR-H3-only variants on the trastuzumab Fv (branch [`Guled`](../../tree/Guled)). The organisers' checker classes them
`novel` (97 % identical to trastuzumab), and 10 of the 15 shortlisted candidates are rejected for CDR-H3 novelty (only 2 mutations from a known
sequence; at least 3 are required). A later version of that shortlist (with the light-chain changes M46L and R66G) is rejected in the same 10 of 15 cases and stays `novel`. They belong to the optimization style of work, so they are not in the de novo list.

What we know and what we do not:

* The de novo VHHs are 63–68 % identical to the trastuzumab VH and pass every format and novelty rule.
* **No scorer we ran on approach A reproduces a designed VHH pose** (OpenFold3 without an MSA misfolds the VHH; OpenDDE docks the VHHs elsewhere on HER2).
  Scores are weak filters, not evidence of binding. A check on the optimization branch's designs shows the same OpenDDE settings recover the crystal pose for paired
  Fvs (14 of 20 within 3 Å) but not for a lone VH (0 of 20 within 10 Å).
* Two of the eight rows from approach B (`dn-K07`, `dn-K08`) do not meet that approach's own epitope/refold selection rule, and `dn-K06` is 91 %
  identical to `dn-K04`; they were taken to fill the slots with checker-valid de novo designs.

## Branches

| branch | contents |
|---|---|
| [`jakob`](../../tree/jakob) | optimization pipeline: enumerate, panel, fold, contacts, pose RMSD, rank (`run_pipeline.py`) |
| [`moritz`](../../tree/moritz) | de novo VHH pipeline (`denovo/`), Claude skill (`.claude/skills/denovo-vhh-pipeline`), result tables and the top-20 structures (`results_denovo/`) |
| [`krishagen`](../../tree/krishagen) | de novo VHH work and the delivery files |
| [`Guled`](../../tree/Guled) | CDR-H3 design code and shortlist |

```bash
git clone https://github.com/jakyric/antibody_hackathon.git
cd antibody_hackathon
git checkout moritz   # or jakob, krishagen, Guled
```
