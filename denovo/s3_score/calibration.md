# Calibration of the scorers

AUC = P(positive scores better than negative), 90% bootstrap interval in brackets. n = positives/negatives.


## af2ig (af2ig)

| family | metric | n | AUC [5%, 95%] |
|---|---|---|---|
| vhh-real | iptm | 10/34 | 0.44 [0.26, 0.61] |
| vhh-real | ipsae | 10/34 | 0.40 [0.23, 0.59] |
| vhh-real | ipsae15 | 10/34 | 0.40 [0.22, 0.58] |
| vhh-real | ipsae20 | 10/34 | 0.45 [0.27, 0.61] |
| vhh-real | ipae_tm | 10/34 | 0.42 [0.25, 0.60] |
| vhh-real | complex_plddt | 10/34 | 0.56 [0.41, 0.73] |
| vhh-real | binder_plddt | 10/34 | 0.56 [0.41, 0.73] |
| vhh-real | min_ipae | 10/34 | 0.41 [0.24, 0.59] |
| vhh-real | mean_ipae | 10/34 | 0.44 [0.27, 0.60] |
| vhh-real | ipae_contacts | 10/34 | 0.44 [0.27, 0.61] |
| vhh-real | interface_pae | 10/34 | 0.44 [0.27, 0.60] |
| vhh-real | ipsae_chn10 | 10/34 | 0.39 [0.22, 0.56] |
| vhh-real | ipsae_chn15 | 10/34 | 0.39 [0.20, 0.57] |
| vhh-real | pose_rmsd_vhh | 10/34 | 0.53 [0.37, 0.71] |
| vhh-hard | iptm | 10/12 | 0.42 [0.20, 0.65] |
| vhh-hard | ipsae | 10/12 | 0.31 [0.10, 0.54] |
| vhh-hard | ipsae15 | 10/12 | 0.34 [0.12, 0.58] |
| vhh-hard | ipsae20 | 10/12 | 0.40 [0.19, 0.63] |
| vhh-hard | ipae_tm | 10/12 | 0.42 [0.21, 0.65] |
| vhh-hard | complex_plddt | 10/12 | 0.55 [0.34, 0.77] |
| vhh-hard | binder_plddt | 10/12 | 0.55 [0.34, 0.77] |
| vhh-hard | min_ipae | 10/12 | 0.33 [0.11, 0.56] |
| vhh-hard | mean_ipae | 10/12 | 0.43 [0.22, 0.65] |
| vhh-hard | ipae_contacts | 10/12 | 0.42 [0.21, 0.64] |
| vhh-hard | interface_pae | 10/12 | 0.43 [0.22, 0.65] |
| vhh-hard | ipsae_chn10 | 10/12 | 0.28 [0.06, 0.51] |
| vhh-hard | ipsae_chn15 | 10/12 | 0.33 [0.12, 0.57] |
| vhh-hard | pose_rmsd_vhh | 10/12 | 0.39 [0.19, 0.59] |

Primary metric for af2ig: **interface_pae** (chosen a priori; vhh-real AUC 0.44; best of all metrics: complex_plddt 0.56, i.e. differences are within noise).


**Gate for af2ig**: interface_pae <= 16.152 (95th percentile of 34 decoys); passes 20% of the 10 cognate crystal VHH pairs.


## dde (opendde-abag)

| family | metric | n | AUC [5%, 95%] |
|---|---|---|---|
| fv-spr | iptm | 12/12 | 0.67 [0.48, 0.85] |
| fv-spr | ipsae | 12/12 | 0.65 [0.46, 0.83] |
| fv-spr | ipsae15 | 12/12 | 0.69 [0.53, 0.85] |
| fv-spr | ipsae20 | 12/12 | 0.69 [0.53, 0.86] |
| fv-spr | ipae_tm | 12/12 | 0.64 [0.45, 0.83] |
| fv-spr | confidence | 12/12 | 0.66 [0.47, 0.85] |
| fv-spr | complex_plddt | 12/12 | 0.59 [0.38, 0.79] |
| fv-spr | binder_plddt | 12/12 | 0.55 [0.32, 0.76] |
| fv-spr | min_ipae | 12/12 | 0.69 [0.49, 0.87] |
| fv-spr | mean_ipae | 12/12 | 0.64 [0.46, 0.83] |
| fv-spr | ipae_contacts | 12/12 | 0.64 [0.45, 0.83] |
| fv-spr | ipsae_chn10 | 12/12 | 0.65 [0.44, 0.85] |
| fv-spr | ipsae_chn15 | 12/12 | 0.56 [0.36, 0.76] |
| fv-facs | iptm | 16/16 | 0.16 [0.05, 0.28] |
| fv-facs | ipsae | 16/16 | 0.14 [0.04, 0.27] |
| fv-facs | ipsae15 | 16/16 | 0.15 [0.04, 0.29] |
| fv-facs | ipsae20 | 16/16 | 0.15 [0.04, 0.28] |
| fv-facs | ipae_tm | 16/16 | 0.27 [0.11, 0.44] |
| fv-facs | confidence | 16/16 | 0.16 [0.05, 0.28] |
| fv-facs | complex_plddt | 16/16 | 0.23 [0.11, 0.39] |
| fv-facs | binder_plddt | 16/16 | 0.26 [0.12, 0.43] |
| fv-facs | min_ipae | 16/16 | 0.18 [0.07, 0.32] |
| fv-facs | mean_ipae | 16/16 | 0.14 [0.04, 0.27] |
| fv-facs | ipae_contacts | 16/16 | 0.27 [0.11, 0.43] |
| fv-facs | ipsae_chn10 | 16/16 | 0.15 [0.04, 0.29] |
| fv-facs | ipsae_chn15 | 16/16 | 0.14 [0.03, 0.28] |
| vhh-real | iptm | 10/34 | 0.65 [0.53, 0.80] |
| vhh-real | ipsae | 10/34 | 0.65 [0.50, 0.79] |
| vhh-real | ipsae15 | 10/34 | 0.65 [0.53, 0.79] |
| vhh-real | ipsae20 | 10/34 | 0.67 [0.54, 0.80] |
| vhh-real | ipae_tm | 10/34 | 0.67 [0.53, 0.81] |
| vhh-real | confidence | 10/34 | 0.65 [0.53, 0.80] |
| vhh-real | complex_plddt | 10/34 | 0.60 [0.46, 0.75] |
| vhh-real | binder_plddt | 10/34 | 0.66 [0.52, 0.80] |
| vhh-real | min_ipae | 10/34 | 0.72 [0.58, 0.87] |
| vhh-real | mean_ipae | 10/34 | 0.65 [0.53, 0.80] |
| vhh-real | ipae_contacts | 10/34 | 0.67 [0.53, 0.81] |
| vhh-real | ipsae_chn10 | 10/34 | 0.66 [0.53, 0.80] |
| vhh-real | ipsae_chn15 | 10/34 | 0.62 [0.46, 0.78] |
| vhh-hard | iptm | 10/12 | 0.37 [0.17, 0.58] |
| vhh-hard | ipsae | 10/12 | 0.38 [0.17, 0.60] |
| vhh-hard | ipsae15 | 10/12 | 0.39 [0.19, 0.59] |
| vhh-hard | ipsae20 | 10/12 | 0.41 [0.22, 0.62] |
| vhh-hard | ipae_tm | 10/12 | 0.46 [0.25, 0.66] |
| vhh-hard | confidence | 10/12 | 0.37 [0.17, 0.57] |
| vhh-hard | complex_plddt | 10/12 | 0.38 [0.20, 0.58] |
| vhh-hard | binder_plddt | 10/12 | 0.43 [0.21, 0.65] |
| vhh-hard | min_ipae | 10/12 | 0.54 [0.30, 0.75] |
| vhh-hard | mean_ipae | 10/12 | 0.41 [0.22, 0.60] |
| vhh-hard | ipae_contacts | 10/12 | 0.45 [0.25, 0.65] |
| vhh-hard | ipsae_chn10 | 10/12 | 0.40 [0.20, 0.62] |
| vhh-hard | ipsae_chn15 | 10/12 | 0.40 [0.20, 0.61] |

Primary metric for dde: **min_ipae** (chosen a priori; vhh-real AUC 0.72; best of all metrics: min_ipae 0.72, i.e. differences are within noise).


**Gate for dde**: min_ipae <= 0.434 (95th percentile of 46 decoys); passes 40% of the 10 cognate crystal VHH pairs.


## esm2 (esmfold2)

| family | metric | n | AUC [5%, 95%] |
|---|---|---|---|
| fv-spr | iptm | 12/12 | 0.70 [0.49, 0.88] |
| fv-spr | ipsae | 12/12 | 0.72 [0.52, 0.88] |
| fv-spr | ipsae15 | 12/12 | 0.74 [0.56, 0.90] |
| fv-spr | ipsae20 | 12/12 | 0.77 [0.60, 0.92] |
| fv-spr | ipae_tm | 12/12 | 0.72 [0.54, 0.88] |
| fv-spr | binder_plddt | 12/12 | 0.65 [0.44, 0.85] |
| fv-spr | min_ipae | 12/12 | 0.76 [0.59, 0.90] |
| fv-spr | mean_ipae | 12/12 | 0.71 [0.49, 0.88] |
| fv-spr | ipae_contacts | 12/12 | 0.71 [0.53, 0.88] |
| fv-spr | ipsae_chn10 | 12/12 | 0.71 [0.51, 0.87] |
| fv-spr | ipsae_chn15 | 12/12 | 0.76 [0.58, 0.90] |
| fv-facs | iptm | 16/16 | 0.38 [0.22, 0.54] |
| fv-facs | ipsae | 16/16 | 0.41 [0.25, 0.57] |
| fv-facs | ipsae15 | 16/16 | 0.40 [0.24, 0.57] |
| fv-facs | ipsae20 | 16/16 | 0.39 [0.23, 0.55] |
| fv-facs | ipae_tm | 16/16 | 0.37 [0.22, 0.54] |
| fv-facs | binder_plddt | 16/16 | 0.36 [0.20, 0.53] |
| fv-facs | min_ipae | 16/16 | 0.39 [0.23, 0.56] |
| fv-facs | mean_ipae | 16/16 | 0.38 [0.20, 0.54] |
| fv-facs | ipae_contacts | 16/16 | 0.37 [0.21, 0.53] |
| fv-facs | ipsae_chn10 | 16/16 | 0.39 [0.23, 0.57] |
| fv-facs | ipsae_chn15 | 16/16 | 0.38 [0.21, 0.54] |
| vhh-real | iptm | 10/34 | 0.62 [0.46, 0.79] |
| vhh-real | ipsae | 10/34 | 0.62 [0.46, 0.76] |
| vhh-real | ipsae15 | 10/34 | 0.62 [0.45, 0.79] |
| vhh-real | ipsae20 | 10/34 | 0.61 [0.44, 0.79] |
| vhh-real | ipae_tm | 10/34 | 0.64 [0.43, 0.82] |
| vhh-real | binder_plddt | 10/34 | 0.61 [0.44, 0.77] |
| vhh-real | min_ipae | 10/34 | 0.64 [0.45, 0.81] |
| vhh-real | mean_ipae | 10/34 | 0.62 [0.45, 0.78] |
| vhh-real | ipae_contacts | 10/34 | 0.63 [0.42, 0.81] |
| vhh-real | ipsae_chn10 | 10/34 | 0.67 [0.49, 0.83] |
| vhh-real | ipsae_chn15 | 10/34 | 0.64 [0.46, 0.79] |
| vhh-hard | iptm | 10/12 | 0.38 [0.16, 0.60] |
| vhh-hard | ipsae | 10/12 | 0.41 [0.19, 0.63] |
| vhh-hard | ipsae15 | 10/12 | 0.41 [0.17, 0.63] |
| vhh-hard | ipsae20 | 10/12 | 0.40 [0.18, 0.63] |
| vhh-hard | ipae_tm | 10/12 | 0.50 [0.25, 0.72] |
| vhh-hard | binder_plddt | 10/12 | 0.42 [0.22, 0.63] |
| vhh-hard | min_ipae | 10/12 | 0.49 [0.24, 0.71] |
| vhh-hard | mean_ipae | 10/12 | 0.40 [0.17, 0.62] |
| vhh-hard | ipae_contacts | 10/12 | 0.49 [0.25, 0.72] |
| vhh-hard | ipsae_chn10 | 10/12 | 0.52 [0.32, 0.73] |
| vhh-hard | ipsae_chn15 | 10/12 | 0.45 [0.24, 0.66] |

Primary metric for esm2: **ipsae_chn15** (chosen a priori; vhh-real AUC 0.64; best of all metrics: ipsae_chn10 0.67, i.e. differences are within noise).


**Gate for esm2**: ipsae_chn15 >= 0.848 (95th percentile of 46 decoys); passes 20% of the 10 cognate crystal VHH pairs.


## of3 (openfold3)

| family | metric | n | AUC [5%, 95%] |
|---|---|---|---|
| fv-spr | iptm | 12/12 | 0.48 [0.30, 0.67] |
| fv-spr | ipsae | 12/12 | 0.56 [0.42, 0.71] |
| fv-spr | ipsae15 | 12/12 | 0.65 [0.49, 0.85] |
| fv-spr | ipsae20 | 12/12 | 0.56 [0.38, 0.77] |
| fv-spr | ipae_tm | 12/12 | 0.80 [0.65, 0.97] |
| fv-spr | confidence | 12/12 | 0.51 [0.33, 0.72] |
| fv-spr | complex_plddt | 12/12 | 0.43 [0.24, 0.65] |
| fv-spr | binder_plddt | 12/12 | 0.40 [0.20, 0.59] |
| fv-spr | min_ipae | 12/12 | 0.57 [0.39, 0.77] |
| fv-spr | mean_ipae | 12/12 | 0.67 [0.49, 0.87] |
| fv-spr | ipae_contacts | 12/12 | 0.80 [0.65, 0.96] |
| fv-spr | ipsae_chn10 | 12/12 | 0.56 [0.42, 0.71] |
| fv-spr | ipsae_chn15 | 12/12 | 0.65 [0.49, 0.84] |
| fv-facs | iptm | 16/16 | 0.35 [0.18, 0.54] |
| fv-facs | ipsae | 16/16 | 0.52 [0.38, 0.66] |
| fv-facs | ipsae15 | 16/16 | 0.44 [0.29, 0.64] |
| fv-facs | ipsae20 | 16/16 | 0.41 [0.27, 0.59] |
| fv-facs | ipae_tm | 16/16 | 0.33 [0.17, 0.52] |
| fv-facs | confidence | 16/16 | 0.35 [0.18, 0.53] |
| fv-facs | complex_plddt | 16/16 | 0.55 [0.37, 0.72] |
| fv-facs | binder_plddt | 16/16 | 0.45 [0.29, 0.63] |
| fv-facs | min_ipae | 16/16 | 0.43 [0.27, 0.63] |
| fv-facs | mean_ipae | 16/16 | 0.34 [0.18, 0.50] |
| fv-facs | ipae_contacts | 16/16 | 0.34 [0.18, 0.52] |
| fv-facs | ipsae_chn10 | 16/16 | 0.51 [0.37, 0.66] |
| fv-facs | ipsae_chn15 | 16/16 | 0.45 [0.29, 0.64] |
| vhh-real | iptm | 10/34 | 0.54 [0.35, 0.74] |
| vhh-real | ipsae | 10/34 | 0.65 [0.50, 0.80] |
| vhh-real | ipsae15 | 10/34 | 0.55 [0.37, 0.74] |
| vhh-real | ipsae20 | 10/34 | 0.57 [0.40, 0.75] |
| vhh-real | ipae_tm | 10/34 | 0.56 [0.39, 0.75] |
| vhh-real | confidence | 10/34 | 0.55 [0.37, 0.74] |
| vhh-real | complex_plddt | 10/34 | 0.52 [0.33, 0.70] |
| vhh-real | binder_plddt | 10/34 | 0.58 [0.42, 0.75] |
| vhh-real | min_ipae | 10/34 | 0.61 [0.40, 0.79] |
| vhh-real | mean_ipae | 10/34 | 0.53 [0.37, 0.73] |
| vhh-real | ipae_contacts | 10/34 | 0.56 [0.39, 0.75] |
| vhh-real | ipsae_chn10 | 10/34 | 0.66 [0.50, 0.81] |
| vhh-real | ipsae_chn15 | 10/34 | 0.59 [0.39, 0.75] |
| vhh-hard | iptm | 10/12 | 0.45 [0.21, 0.70] |
| vhh-hard | ipsae | 10/12 | 0.55 [0.37, 0.75] |
| vhh-hard | ipsae15 | 10/12 | 0.42 [0.20, 0.65] |
| vhh-hard | ipsae20 | 10/12 | 0.45 [0.22, 0.68] |
| vhh-hard | ipae_tm | 10/12 | 0.45 [0.22, 0.70] |
| vhh-hard | confidence | 10/12 | 0.43 [0.20, 0.66] |
| vhh-hard | complex_plddt | 10/12 | 0.40 [0.18, 0.60] |
| vhh-hard | binder_plddt | 10/12 | 0.47 [0.27, 0.68] |
| vhh-hard | min_ipae | 10/12 | 0.51 [0.28, 0.73] |
| vhh-hard | mean_ipae | 10/12 | 0.38 [0.17, 0.62] |
| vhh-hard | ipae_contacts | 10/12 | 0.45 [0.22, 0.70] |
| vhh-hard | ipsae_chn10 | 10/12 | 0.57 [0.38, 0.77] |
| vhh-hard | ipsae_chn15 | 10/12 | 0.48 [0.26, 0.69] |

Primary metric for of3: **ipsae_chn15** (chosen a priori; vhh-real AUC 0.59; best of all metrics: ipsae_chn10 0.66, i.e. differences are within noise).


**Gate for of3**: ipsae_chn15 >= 0.474 (95th percentile of 77 decoys); passes 30% of the 10 cognate crystal VHH pairs.


## Per generator

| generator | model | n | median primary metric | share past gate (incl. epitope overlap) | median epitope overlap |
|---|---|---|---|---|---|
| control | af2ig | 1 | 25.935 | 0% | nan |
| control | dde | 14 | 7.408 | 0% | 0.00 |
| control | esm2 | 14 | 0.150 | 0% | 0.27 |
| control | of3 | 14 | 0.071 | 0% | 0.23 |
| rfd4jf | dde | 31 | 4.917 | 0% | 0.00 |
| rfd4jf | esm2 | 31 | 0.145 | 0% | 0.29 |
| rfd4jf | of3 | 305 | 0.138 | 0% | 0.21 |

## AF2-IG (re-scores a given pose; deterministic, no MSA)

| group | n | iptm median [max] | interface PAE median [min] | min_ipae median | ipsae_chn15 median | pose RMSD vs input, VHH median [A] |
|---|---|---|---|---|---|---|
| calib:vhhhard | 12 | 0.32 [0.72] | 21.0 [9.0] | 6.5 | 0.26 | 5.8 |
| calib:vhhneg | 22 | 0.22 [0.62] | 23.3 [12.1] | 9.2 | 0.24 | 20.6 |
| calib:vhhpos | 10 | 0.17 [0.74] | 24.5 [10.9] | 10.6 | 0.21 | 14.4 |
| control(OF3 pose) | 1 | 0.11 [0.11] | 25.9 [25.9] | 16.2 | 0.00 | 36.5 |

Caveat: designs are scored on the generator's own pose, controls and calibration decoys on an OpenFold3 pose, so a gap between the groups mixes pose quality with scorer behaviour; the calibration rows (cognate vs shuffled, both on OF3 poses) isolate the scorer.

