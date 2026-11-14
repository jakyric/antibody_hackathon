# Pipeline status, 16:53:47 CET

## Funnel (designs that reached each stage)

| generator | designs | rules | of3 | orth | props | gate_vm | gate_strict | composite |
|---|---|---|---|---|---|---|---|---|
| rfd4jf | 1408 | 1408 | 1153 | 1153 | 1 | 723 | 43 | 737 |
| rfd4jf-loop | 136 | 72 | 62 | 62 | 0 | 0 | 0 | 3 |

## OF3 metrics: controls vs designs (median)

| kind | of3_iptm | of3_ipsae | of3_min_ipae | of3_epi_overlap |
|---|---|---|---|---|
| control | 0.134 | 0.0 | 15.443 | 0.216 |
| design | 0.117 | 0.0 | 14.563 | 0.222 |

### Control rows

| design_id | of3_iptm | of3_ipsae | of3_min_ipae | of3_epi_overlap |
|---|---|---|---|---|
| ctl-bcii10-wt | 0.112 | 0.0 | 16.054 | 0.0 |
| ctl-fwrand-01 | 0.139 | 0.0 | 16.052 | 0.25 |
| ctl-fwrand-02 | 0.121 | 0.0 | 18.734 | 0.0 |
| ctl-fwrand-03 | 0.129 | 0.0 | 12.677 | 0.222 |
| ctl-fwrand-04 | 0.144 | 0.0 | 13.303 | 0.389 |
| ctl-fwrand-05 | 0.112 | 0.0 | 15.945 | 0.2 |
| ctl-fwrand-06 | 0.106 | 0.0 | 16.798 | 0.083 |
| ctl-fwrand-07 | 0.138 | 0.0 | 13.453 | 0.235 |
| ctl-fwrand-08 | 0.169 | 0.0 | 13.49 | 0.308 |
| ctl-fwrand-09 | 0.155 | 0.0 | 15.895 | 0.211 |
| ctl-fwrand-10 | 0.194 | 0.0 | 13.409 | 0.261 |
| ctl-fwrand-11 | 0.167 | 0.0 | 12.899 | 0.273 |
| ctl-fwrand-12 | 0.088 | 0.0 | 16.668 | 0.067 |
| ctl-kri-01 | 0.089 | 0.0 | 15.146 | 0.176 |
| ctl-kri-02 | 0.546 | 0.128 | 3.763 | 0.0 |
| ctl-kri-03 | 0.112 | 0.0 | 17.513 | 0.0 |
| ctl-kri-04 | 0.105 | 0.0 | 15.74 | 0.3 |
| ctl-kri-05 | 0.088 | 0.0 | 16.365 | 0.24 |
| ctl-kri-06 | 0.149 | 0.0 | 13.797 | 0.0 |
| ctl-tz-fv | 0.555 | 0.0 | 10.2 | 0.72 |

## Current top 10 by composite

| design_id | framework | cdr3 | composite | t_bind | iptm_mean | ipsae_mean | of3_iptm | dde_iptm | esm2_iptm | af2ig_iptm | ptx_iptm | of3_ipsae | dde_ipsae | esm2_ipsae | af2ig_ipsae | vmof3_pair_iptm | vmof3_ipsae_chn15 | vmesm2_pair_iptm | vmesm2_ipsae_chn15 | ddex_pdockq2 | orth_dde_own_vhh_rmsd_fold | orth_dde_own_vhh_rmsd_her2frame | tnp_n_red | tnp_n_amber | my_yield | nm_tm | nb9_score | nb9_of | gate_dde | gate_strict | n_flags | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bcii10-s401032-002 | bcii10 | 18 | 0.86 | 0.85 | 0.29 | 0.08 | 0.1 |  |  |  |  | 0.0 |  |  |  | 0.24 | 0.0 |  |  |  |  |  | 0.0 | 0.0 | 0.22 | 70.75 |  |  | 1 | 1 | 2 | below_decoy;clash |
| capla-s408040-006 | capla | 18 | 0.84 | 0.87 | 0.25 | 0.16 | 0.16 | 0.88 |  |  |  | 0.15 | 0.8 |  |  | 0.25 | 0.15 |  |  | 0.48 | 1.48 | 51.53 | 0.0 | 0.0 | -0.56 | 71.6 |  |  | 1 | 1 | 2 | below_decoy;dde_fold_rmsd |
| capla-s407024-000 | capla | 14 | 0.83 | 0.94 | 0.51 | 0.39 | 0.12 | 0.39 |  |  |  | 0.0 | 0.2 |  |  | 0.27 | 0.14 |  |  | 0.02 | 0.47 | 17.73 | 0.0 | 0.0 | -0.79 | 71.55 |  |  | 1 | 1 | 1 | below_decoy |
| capla-s307000-000 | capla | 18 | 0.82 | 0.81 | 0.19 | 0.17 | 0.08 |  |  |  |  | 0.0 |  |  |  | 0.22 | 0.18 |  |  |  |  |  | 0.0 | 0.0 | -0.07 | 70.49 |  |  | 1 | 1 | 2 | below_decoy;cys_in_cdr |
| capla-s100612-007 | capla | 18 | 0.81 | 0.84 | 0.28 | 0.08 | 0.11 | 0.41 | 0.38 | 0.46 |  | 0.16 | 0.27 | 0.18 | 0.52 | 0.16 | 0.0 |  |  | 0.04 | 0.66 | 41.47 | 0.0 | 0.0 | 0.09 | 69.68 |  |  | 1 | 1 | 1 | below_decoy |
| bcii10-s201008-000 | bcii10 | 17 | 0.81 | 0.7 | 0.19 | 0.13 | 0.1 |  |  |  |  | 0.0 |  |  |  | 0.24 | 0.14 |  |  |  |  |  | 0.0 | 1.0 | 0.49 | 72.59 |  |  | 1 | 1 | 2 | below_decoy;clash |
| capla-s101014-007 | capla | 17 | 0.8 | 0.9 | 0.41 | 0.35 | 0.17 | 0.22 |  | 0.53 |  | 0.0 | 0.18 |  | 0.49 | 0.59 | 0.55 |  |  | 0.01 | 0.7 | 49.11 | 0.0 | 0.0 | -0.54 | 69.76 |  |  | 1 | 1 | 1 | below_decoy |
| capla-s411008-004 | capla | 14 | 0.76 | 0.97 | 0.44 | 0.22 | 0.1 |  |  |  |  | 0.0 |  |  |  | 0.42 | 0.21 |  |  |  |  |  | 0.0 | 0.0 | -0.64 | 68.42 |  |  | 1 | 1 | 3 | below_decoy;clash;dde_fold_rmsd |
| capla-s411008-006 | capla | 18 | 0.75 | 0.64 | 0.17 | 0.12 | 0.29 | 0.93 |  | 0.53 |  | 0.22 | 0.81 |  | 0.49 | 0.1 | 0.0 | 0.27 | 0.21 | 0.7 | 2.07 | 48.87 | 0.0 | 0.0 | -0.31 | 71.42 |  |  | 1 | 1 | 2 | below_decoy;dde_fold_rmsd |
| capla-s409024-000 | capla | 17 | 0.75 | 0.64 | 0.16 | 0.14 | 0.2 | 0.75 |  |  |  | 0.14 | 0.57 |  |  | 0.15 | 0.14 | 0.17 | 0.16 | 0.44 | 1.69 | 42.02 | 0.0 | 0.0 | -0.56 | 72.42 |  |  | 1 | 1 | 1 | below_decoy |

## Pareto (last step, `denovo/pareto.py`)

Population 43 gate-passing designs; objectives bind,dev,yield,therm; front sizes {1: 8, 2: 11, 3: 12, 4: 6}; fully measured 43.

JapanFold balance: $76.53
