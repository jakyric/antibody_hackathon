# Pipeline status, 15:00:52 CET

## Funnel (designs that reached each stage)

| generator | designs | rules | of3 | orth | props | gate_vm | gate_strict | composite |
|---|---|---|---|---|---|---|---|---|
| rfd4jf | 449 | 449 | 367 | 367 | 449 | 386 | 283 | 449 |

## OF3 metrics: controls vs designs (median)

| kind | of3_iptm | of3_ipsae | of3_min_ipae | of3_epi_overlap |
|---|---|---|---|---|
| control | 0.138 | 0.0 | 14.692 | 0.229 |
| design | 0.115 | 0.0 | 14.553 | 0.214 |

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
| ctl-tz-fv | 0.555 | 0.0 | 10.2 | 0.72 |

## Current top 10 by composite

| metric | capla-s100600-003 | bcii10-s100500-000 | capla-s100600-004 | bcii10-s2101-007 | capla-s100900-000 | capla-s3100-006 | bcii10-s100200-002 | capla-s101000-002 | capla-s100800-007 | capla-s101003-006 |
|---|---|---|---|---|---|---|---|---|---|---|
| composite | 0.839 | 0.811 | 0.745 | 0.726 | 0.722 | 0.722 | 0.716 | 0.713 | 0.705 | 0.701 |
| coverage | 0.857 | 0.857 | 0.857 | 1 | 0.857 | 1 | 0.857 | 0.857 | 0.857 | 0.857 |
| t_bind | 0.98 | 0.992 | 0.786 | 0.74 | 0.783 | 0.852 | 0.902 | 0.98 | 0.94 | 0.97 |
| t_comp | 0.749 | 0.595 | 0.913 | 0.594 | 0.497 | 0.586 | 0.243 | 0.337 | 0.598 | 0.243 |
| t_therm | 0.959 | 0.9 | 0.803 | 0.965 | 0.926 | 0.45 | 0.924 | 0.679 | 0.624 | 0.894 |
| t_yield |  |  |  | 0.63 |  | 0.593 |  |  |  |  |
| t_dev | 1 | 1 | 0.85 | 0.5 | 1 | 1 | 0.85 | 0.85 | 0.85 | 0.55 |
| t_rules | 0.8 | 0.8 | 0.8 | 0.8 | 0.8 | 1 | 0.8 | 0.8 | 0.4 | 0.667 |
| of3_consensus_rank_score | 0.981 | 0.992 | 0.814 | 0.915 | 0.723 | 0.908 | 0.898 | 0.977 | 0.929 | 0.995 |
| of3_binder_plddt | 84.3 | 89.7 | 84.9 | 87.5 | 85.7 | 87 | 87.4 | 87.5 | 86.1 | 90.9 |
| of3_cdr3_cent_max | 7.92 | 4.27 | 31.1 | 31.6 | 39.9 | 35.5 | 29.4 | 2.86 | 3.34 | 1.19 |
| of3_cdr3_plddt | 51.1 | 76.2 | 61 | 69.2 | 62.2 | 68.3 | 76.8 | 67.3 | 55.3 | 79.4 |
| of3_cdr_fv_frac4 | 0.197 | 0.128 | 0.764 | 0.117 | 0.179 | 0.157 | 0 | 0.0595 | 0.0948 | 0 |
| of3_clash25 | 50 | 23 | 595 | 21 | 134 | 25 | 0 | 8 | 20 | 0 |
| of3_clash30 | 61 | 28 | 647 | 26 | 155 | 31 | 0 | 13 | 23 | 0 |
| of3_complex_plddt | 0.811 | 0.864 | 0.818 | 0.832 | 0.816 | 0.824 | 0.824 | 0.837 | 0.812 | 0.86 |
| of3_confidence | 0.434 | 0.576 | 0.261 | 0.34 | 0.208 | 0.268 | 0.29 | 0.389 | 0.299 | 0.529 |
| of3_cov_hs1_core | 0 | 0.5 | 1 | 0.167 | 0.333 | 0.167 | 0 | 0.333 | 0 | 0 |
| of3_cov_hs2_broad | 0.167 | 0.333 | 0.75 | 0.25 | 0.333 | 0.25 | 0 | 0.167 | 0.167 | 0 |
| of3_cov_hs3_alt_patch | 0.5 | 0.333 | 0.5 | 0.5 | 0 | 0.333 | 0 | 0 | 0.5 | 0 |
| of3_epi_coverage | 0.368 | 0.421 | 0.579 | 0.368 | 0.211 | 0.316 | 0 | 0.105 | 0.263 | 0 |
| of3_epi_overlap | 0.583 | 0.333 | 0.579 | 0.35 | 0.19 | 0.3 | 0 | 0.0909 | 0.385 | 0 |
| of3_epi_overlap_min | 0.389 | 0.304 | 0 | 0 | 0 | 0.263 | 0 | 0.0909 | 0.312 | 0 |
| of3_fv_contact40 | 70 | 40 | 727 | 36 | 203 | 45 | 0 | 21 | 31 | 0 |
| of3_her2_fit_rmsd | 1.22 | 1.29 | 1.5 | 1.21 | 1.43 | 1.32 | 1.63 | 1.49 | 1.13 | 1.17 |
| of3_hewl_ipsae15 |  |  |  | 0.0554 |  | 0.0942 |  |  |  |  |
| of3_hewl_iptm |  |  |  | 0.224 |  | 0.308 |  |  |  |  |
| of3_hewl_min_ipae |  |  |  | 13.3 |  | 10.2 |  |  |  |  |
| of3_hewl_n_contact_pairs |  |  |  | 61 |  | 49 |  |  |  |  |
| of3_ipae_cdr_epitope | 15.8 | 11.1 | 19.7 | 17.5 | 25.3 | 19.3 | 22.3 | 18.9 | 18.4 | 15.5 |
| of3_ipae_contacts | 13.1 | 8.32 | 18.6 | 14.7 | 22.1 | 17.8 | 16.6 | 13.9 | 15.9 | 7.83 |
| of3_ipae_tm | 0.286 | 0.495 | 0.16 | 0.235 | 0.118 | 0.171 | 0.194 | 0.272 | 0.205 | 0.518 |
| of3_ipsae | 0.0161 | 0.16 | 0 | 0 | 0 | 0.0129 | 0 | 0.0335 | 0.0124 | 0.174 |
| of3_ipsae15 | 0.0543 | 0.182 | 0.00579 | 0.0428 | 0.00627 | 0.0126 | 0.01 | 0.112 | 0.0182 | 0.186 |
| of3_ipsae20 | 0.0733 | 0.184 | 0.0403 | 0.0663 | 0.00917 | 0.0345 | 0.0499 | 0.117 | 0.0472 | 0.187 |
| of3_ipsae_chn10 | 0.334 | 0.381 | 0 | 0 | 0 | 0.284 | 0 | 0.305 | 0.28 | 0.389 |
| of3_ipsae_chn15 | 0.234 | 0.329 | 0.154 | 0.171 | 0.163 | 0.191 | 0.171 | 0.238 | 0.199 | 0.357 |
| of3_iptm | 0.387 | 0.548 | 0.19 | 0.264 | 0.126 | 0.197 | 0.215 | 0.344 | 0.233 | 0.491 |
| of3_iptm_sd | 0.0783 | 0.038 | 0.039 | 0.0635 | 0.0175 | 0.0525 | 0.0496 | 0.0474 | 0.00455 | 0.0381 |
| of3_mean_ipae | 18.7 | 13.6 | 22.2 | 20.4 | 25 | 22.9 | 21.6 | 19.5 | 21.6 | 14.3 |
| of3_min_ipae | 5.92 | 4.9 | 12 | 10.7 | 11.7 | 8.67 | 10.6 | 6.92 | 9.26 | 4.82 |
| of3_n_contact_pairs | 33 | 63 | 45 | 51 | 48 | 42 | 39 | 38 | 31 | 32 |
| of3_n_her2_res_contact | 12 | 24 | 19 | 20 | 21 | 20 | 16 | 22 | 13 | 13 |
| of3_n_vhh_res_contact | 18 | 27 | 20 | 25 | 22 | 15 | 22 | 15 | 18 | 18 |
| of3_own_vhh_rmsd_fold | 3.58 | 0.905 | 0.963 | 0.852 | 0.791 | 0.902 | 0.832 | 3.58 | 0.776 | 0.467 |
| of3_own_vhh_rmsd_her2frame | 49 | 44.3 | 11.1 | 41.2 | 34.9 | 51.7 | 47.3 | 39 | 49.5 | 40 |
| of3_pass | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| of3_pct | 0.98 | 0.992 | 0.744 | 0.905 | 0.874 | 0.943 | 0.908 | 0.986 | 0.955 | 0.998 |
| of3_pose_rmsd_max | 6.87 | 9.14 | 47.5 | 66 | 54.7 | 52 | 57 | 8.15 | 4.84 | 5.47 |
| of3_ptm | 0.601 | 0.688 | 0.549 | 0.59 | 0.525 | 0.553 | 0.573 | 0.569 | 0.562 | 0.68 |
| of3_seed_cdr3_dist |  |  |  |  |  |  |  |  |  | 27.4 |
| of3_seed_consistent |  |  |  |  |  |  |  |  |  | 0 |
| of3_seed_iptm_min |  |  |  |  |  |  |  |  |  | 0.18 |
| of3_seed_vhh_rmsd |  |  |  |  |  |  |  |  |  | 39.2 |
| of3_spec_iptm_delta |  |  |  | 0.0408 |  | -0.111 |  |  |  |  |
| orth_parent_id |  |  |  |  |  |  |  |  |  |  |
| orth_is_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| orth_denovo_vm_pass | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| orth_tier | 1 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | 1 | 2 |
| orth_consensus_rank_score | 0.98 | 0.992 | 0.744 | 0.477 | 0.874 | 0.769 | 0.908 | 0.986 | 0.955 | 0.932 |
| orth_n_models_pass | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| orth_n_models_scored | 1 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | 1 | 3 |
| orth_of3_consensus_rank_score | 0.981 | 0.992 | 0.814 | 0.915 | 0.723 | 0.908 | 0.898 | 0.977 | 0.929 | 0.995 |
| orth_af2ig_binder_plddt |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_complex_plddt |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_confidence |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_interface_pae |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipae_cdr_epitope |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipae_contacts |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipae_tm |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipsae |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipsae15 |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipsae20 |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipsae_chn10 |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ipsae_chn15 |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_iptm |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_mean_ipae |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_min_ipae |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_n_contact_pairs |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_n_her_res |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_pass |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_pct |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_pose_rmsd_cdr3 |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_pose_rmsd_vhh |  |  |  |  |  |  |  |  |  |  |
| orth_af2ig_ptm |  |  |  |  |  |  |  |  |  |  |
| orth_dde_binder_plddt |  |  |  | 91.3 |  | 92.3 |  |  |  | 96.5 |
| orth_dde_cdr3_cent_max |  |  |  | 52.7 |  | 21.5 |  |  |  | 2.41 |
| orth_dde_cdr3_plddt |  |  |  | 78.2 |  | 81.2 |  |  |  | 93.3 |
| orth_dde_cdr_fv_frac4 |  |  |  | 0 |  | 0.99 |  |  |  | 0 |
| orth_dde_clash25 |  |  |  | 0 |  | 834 |  |  |  | 0 |
| orth_dde_clash30 |  |  |  | 0 |  | 882 |  |  |  | 0 |
| orth_dde_complex_plddt |  |  |  | 0.869 |  | 0.866 |  |  |  | 0.926 |
| orth_dde_confidence |  |  |  | 0.288 |  | 0.38 |  |  |  | 0.762 |
| orth_dde_cov_hs1_core |  |  |  | 0 |  | 0.833 |  |  |  | 0 |
| orth_dde_cov_hs2_broad |  |  |  | 0 |  | 0.583 |  |  |  | 0 |
| orth_dde_cov_hs3_alt_patch |  |  |  | 0 |  | 0.333 |  |  |  | 0 |
| orth_dde_epi_coverage |  |  |  | 0 |  | 0.474 |  |  |  | 0 |
| orth_dde_epi_overlap |  |  |  | 0 |  | 0.529 |  |  |  | 0 |
| orth_dde_epi_overlap_min |  |  |  | 0 |  | 0.2 |  |  |  | 0 |
| orth_dde_fv_contact40 |  |  |  | 0 |  | 923 |  |  |  | 0 |
| orth_dde_her2_fit_rmsd |  |  |  | 1.38 |  | 1.08 |  |  |  | 1.39 |
| orth_dde_hewl_ipsae15 |  |  |  |  |  |  |  |  |  |  |
| orth_dde_hewl_iptm |  |  |  |  |  |  |  |  |  |  |
| orth_dde_hewl_min_ipae |  |  |  |  |  |  |  |  |  |  |
| orth_dde_hewl_n_contact_pairs |  |  |  |  |  |  |  |  |  |  |
| orth_dde_ipae_cdr_epitope |  |  |  | 25.4 |  | 16.6 |  |  |  | 10.1 |
| orth_dde_ipae_contacts |  |  |  | 19.4 |  | 14.3 |  |  |  | 3.49 |
| orth_dde_ipae_tm |  |  |  | 0.149 |  | 0.277 |  |  |  | 0.837 |
| orth_dde_ipsae |  |  |  | 0 |  | 0.0439 |  |  |  | 0.503 |
| orth_dde_ipsae15 |  |  |  | 0.00782 |  | 0.107 |  |  |  | 0.498 |
| orth_dde_ipsae20 |  |  |  | 0.0549 |  | 0.114 |  |  |  | 0.492 |
| orth_dde_ipsae_chn10 |  |  |  | 0 |  | 0.554 |  |  |  | 0.69 |
| orth_dde_ipsae_chn15 |  |  |  | 0.198 |  | 0.222 |  |  |  | 0.675 |
| orth_dde_iptm |  |  |  | 0.21 |  | 0.315 |  |  |  | 0.744 |
| orth_dde_iptm_sd |  |  |  | 0.00701 |  | 0.0475 |  |  |  | 0.00177 |
| orth_dde_mean_ipae |  |  |  | 24.8 |  | 20 |  |  |  | 9.11 |
| orth_dde_min_ipae |  |  |  | 11.7 |  | 4.01 |  |  |  | 1 |
| orth_dde_n_contact_pairs |  |  |  | 42 |  | 33 |  |  |  | 48 |
| orth_dde_n_her2_res_contact |  |  |  | 20 |  | 17 |  |  |  | 19 |
| orth_dde_n_vhh_res_contact |  |  |  | 20 |  | 18 |  |  |  | 20 |
| orth_dde_own_vhh_rmsd_fold |  |  |  | 2.72 |  | 1.09 |  |  |  | 0.611 |
| orth_dde_own_vhh_rmsd_her2frame |  |  |  | 48.8 |  | 18.7 |  |  |  | 38.2 |
| orth_dde_pass |  |  |  | 0 |  | 0 |  |  |  | 0 |
| orth_dde_pct |  |  |  | 0.0479 |  | 0.596 |  |  |  | 0.952 |
| orth_dde_pose_rmsd_max |  |  |  | 69.3 |  | 47.4 |  |  |  | 3.33 |
| orth_dde_ptm |  |  |  | 0.602 |  | 0.643 |  |  |  | 0.836 |
| orth_esm2_binder_plddt |  |  |  |  |  |  |  |  |  | 86.4 |
| orth_esm2_cdr3_cent_max |  |  |  |  |  |  |  |  |  | 5.26 |
| orth_esm2_cdr3_plddt |  |  |  |  |  |  |  |  |  | 76.5 |
| orth_esm2_cdr_fv_frac4 |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_clash25 |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_clash30 |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_complex_plddt |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_confidence |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_cov_hs1_core |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_cov_hs2_broad |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_cov_hs3_alt_patch |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_epi_coverage |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_epi_overlap |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_epi_overlap_min |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_fv_contact40 |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_her2_fit_rmsd |  |  |  |  |  |  |  |  |  | 1.45 |
| orth_esm2_hewl_ipsae15 |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_hewl_iptm |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_hewl_min_ipae |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_hewl_n_contact_pairs |  |  |  |  |  |  |  |  |  |  |
| orth_esm2_ipae_cdr_epitope |  |  |  |  |  |  |  |  |  | 24.1 |
| orth_esm2_ipae_contacts |  |  |  |  |  |  |  |  |  | 16.1 |
| orth_esm2_ipae_tm |  |  |  |  |  |  |  |  |  | 0.2 |
| orth_esm2_ipsae |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_ipsae15 |  |  |  |  |  |  |  |  |  | 0.00653 |
| orth_esm2_ipsae20 |  |  |  |  |  |  |  |  |  | 0.0418 |
| orth_esm2_ipsae_chn10 |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_ipsae_chn15 |  |  |  |  |  |  |  |  |  | 0.172 |
| orth_esm2_iptm |  |  |  |  |  |  |  |  |  | 0.263 |
| orth_esm2_iptm_sd |  |  |  |  |  |  |  |  |  | 0.00318 |
| orth_esm2_mean_ipae |  |  |  |  |  |  |  |  |  | 22.6 |
| orth_esm2_min_ipae |  |  |  |  |  |  |  |  |  | 12.8 |
| orth_esm2_n_contact_pairs |  |  |  |  |  |  |  |  |  | 39 |
| orth_esm2_n_her2_res_contact |  |  |  |  |  |  |  |  |  | 14 |
| orth_esm2_n_vhh_res_contact |  |  |  |  |  |  |  |  |  | 20 |
| orth_esm2_own_vhh_rmsd_fold |  |  |  |  |  |  |  |  |  | 0.427 |
| orth_esm2_own_vhh_rmsd_her2frame |  |  |  |  |  |  |  |  |  | 40.5 |
| orth_esm2_pass |  |  |  |  |  |  |  |  |  | 0 |
| orth_esm2_pct |  |  |  |  |  |  |  |  |  | 0.845 |
| orth_esm2_pose_rmsd_max |  |  |  |  |  |  |  |  |  | 7.09 |
| orth_esm2_ptm |  |  |  |  |  |  |  |  |  | 0.589 |
| orth_of3_pass | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| orth_of3_pct | 0.98 | 0.992 | 0.744 | 0.905 | 0.874 | 0.943 | 0.908 | 0.986 | 0.955 | 0.998 |
| orth_of3s1_binder_plddt |  |  |  |  |  |  |  |  |  | 70 |
| orth_of3s1_cdr3_cent_max |  |  |  |  |  |  |  |  |  | 2.43 |
| orth_of3s1_cdr3_plddt |  |  |  |  |  |  |  |  |  | 49.7 |
| orth_of3s1_cdr_fv_frac4 |  |  |  |  |  |  |  |  |  | 0 |
| orth_of3s1_clash25 |  |  |  |  |  |  |  |  |  | 0 |
| orth_of3s1_clash30 |  |  |  |  |  |  |  |  |  | 1 |
| orth_of3s1_complex_plddt |  |  |  |  |  |  |  |  |  | 0.626 |
| orth_of3s1_confidence |  |  |  |  |  |  |  |  |  | -99.8 |
| orth_of3s1_cov_hs1_core |  |  |  |  |  |  |  |  |  | 1 |
| orth_of3s1_cov_hs2_broad |  |  |  |  |  |  |  |  |  | 1 |
| orth_of3s1_cov_hs3_alt_patch |  |  |  |  |  |  |  |  |  | 1 |
| orth_of3s1_epi_coverage |  |  |  |  |  |  |  |  |  | 1 |
| orth_of3s1_epi_overlap |  |  |  |  |  |  |  |  |  | 0.161 |
| orth_of3s1_epi_overlap_min |  |  |  |  |  |  |  |  |  | 0.16 |
| orth_of3s1_fv_contact40 |  |  |  |  |  |  |  |  |  | 2 |
| orth_of3s1_her2_fit_rmsd |  |  |  |  |  |  |  |  |  | 17.6 |
| orth_of3s1_ipae_cdr_epitope |  |  |  |  |  |  |  |  |  | 25.9 |
| orth_of3s1_ipae_contacts |  |  |  |  |  |  |  |  |  | 23.9 |
| orth_of3s1_ipae_tm |  |  |  |  |  |  |  |  |  | 0.105 |
| orth_of3s1_ipsae |  |  |  |  |  |  |  |  |  | 0.0107 |
| orth_of3s1_ipsae15 |  |  |  |  |  |  |  |  |  | 0.0149 |
| orth_of3s1_ipsae20 |  |  |  |  |  |  |  |  |  | 0.0466 |
| orth_of3s1_ipsae_chn10 |  |  |  |  |  |  |  |  |  | 0.255 |
| orth_of3s1_ipsae_chn15 |  |  |  |  |  |  |  |  |  | 0.206 |
| orth_of3s1_iptm |  |  |  |  |  |  |  |  |  | 0.18 |
| orth_of3s1_iptm_sd |  |  |  |  |  |  |  |  |  | 0.0102 |
| orth_of3s1_mean_ipae |  |  |  |  |  |  |  |  |  | 23.8 |
| orth_of3s1_min_ipae |  |  |  |  |  |  |  |  |  | 9.99 |
| orth_of3s1_n_contact_pairs |  |  |  |  |  |  |  |  |  | 1.12e+04 |
| orth_of3s1_n_her2_res_contact |  |  |  |  |  |  |  |  |  | 118 |
| orth_of3s1_n_vhh_res_contact |  |  |  |  |  |  |  |  |  | 130 |
| orth_of3s1_own_vhh_rmsd_fold |  |  |  |  |  |  |  |  |  | 13.5 |
| orth_of3s1_own_vhh_rmsd_her2frame |  |  |  |  |  |  |  |  |  | 31.5 |
| orth_of3s1_pose_rmsd_max |  |  |  |  |  |  |  |  |  | 5.97 |
| orth_of3s1_ptm |  |  |  |  |  |  |  |  |  | 0.466 |
| orth_x_dde_vs_of3_cdr3_dist |  |  |  | 37.1 |  | 28.3 |  |  |  | 2.59 |
| orth_x_dde_vs_of3_vhh_rmsd |  |  |  | 66.1 |  | 42.5 |  |  |  | 9.14 |
| orth_x_esm2_vs_of3_cdr3_dist |  |  |  |  |  |  |  |  |  | 2.12 |
| orth_x_esm2_vs_of3_vhh_rmsd |  |  |  |  |  |  |  |  |  | 3.24 |
| orth_xmodel_cdr3_dist_max |  |  |  | 37.1 |  | 28.3 |  |  |  | 2.59 |
| orth_xmodel_n_agree |  |  |  | 0 |  | 0 |  |  |  | 2 |
| orth_xmodel_vhh_rmsd_mean |  |  |  | 66.1 |  | 42.5 |  |  |  | 6.19 |
| nm_tm | 73.8 | 73.1 | 72.3 | 74 | 73.3 | 70.1 | 73.3 | 71.4 | 71 | 73 |
| my_yield |  |  |  | -0.7 |  | -0.739 |  |  |  |  |
| my_yield_sd |  |  |  | 0.136 |  | 0.182 |  |  |  |  |
| my_mdratio |  |  |  | 0.559 |  | -0.21 |  |  |  |  |
| my_mdratio_sd |  |  |  | 0.099 |  | 0.326 |  |  |  |  |
| tnp_cdr_len |  |  |  | 38 |  | 30 |  |  |  |  |
| tnp_cdr3_len |  |  |  | 17 |  | 13 |  |  |  |  |
| tnp_cdr3_compact |  |  |  | 1.35 |  | 0.958 |  |  |  |  |
| tnp_psh |  |  |  | 156 |  | 88.5 |  |  |  |  |
| tnp_ppc |  |  |  | 0 |  | 0 |  |  |  |  |
| tnp_pnc |  |  |  | 0 |  | 0.0361 |  |  |  |  |
| tnp_n_red |  |  |  | 1 |  | 0 |  |  |  |  |
| tnp_n_amber |  |  |  | 1 |  | 0 |  |  |  |  |
| tnp_error |  |  |  |  |  |  |  |  |  |  |
| liab_len | 128 | 129 | 127 | 128 | 125 | 120 | 131 | 128 | 124 | 130 |
| liab_n_cys | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| liab_canon_cys | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| liab_unpaired_cys | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_nglyc | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| liab_nglyc_cdr | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| liab_ng | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_ng_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_ns_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_dg | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_dg_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_dp_cdr | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| liab_met_cdr | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_trp_cdr | 1 | 0 | 0 | 1 | 1 | 1 | 0 | 1 | 2 | 3 |
| liab_pI | 5.86 | 6.34 | 6.49 | 6.34 | 8.01 | 8.01 | 5.17 | 4.96 | 6.85 | 5.04 |
| liab_charge74 | -2.29 | -0.51 | -0.4 | -0.51 | 0.59 | 0.59 | -1.51 | -2.4 | -0.37 | -2.41 |
| liab_cdr_pos | 2 | 1 | 2 | 1 | 3 | 2 | 2 | 0 | 0 | 3 |
| liab_cdr_neg | 5 | 2 | 3 | 2 | 3 | 2 | 4 | 4 | 1 | 6 |
| liab_cdr_gravy | -0.853 | -0.177 | -0.222 | -0.163 | -0.591 | -0.7 | 0.137 | -0.105 | -0.479 | -0.81 |
| liab_cdr3_gravy | -1.13 | -0.49 | -0.441 | -0.376 | -1.11 | -1.22 | -0.338 | 0.09 | -0.271 | -1.03 |
| liab_max_hydro_win5 | 0.52 | 0.82 | 1.02 | 1.08 | 0.7 | 0.62 | 1.56 | 1.96 | 0.26 | 0.1 |
| liab_cdr_hyd_run | 1 | 2 | 2 | 2 | 1 | 1 | 1 | 1 | 2 | 2 |
| liab_f_unpaired_cys | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_f_nglyc | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
| liab_f_ng_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_f_dg_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_f_met_trp_cdr | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| liab_f_hydrophobic | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_f_pI | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 1 |
| liab_f_cdr_pos | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_n_flags | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 1 | 1 | 3 |
| liab_no_canon_ss | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| liab_f_no_canon_ss | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| nb9_score | 4 | 4 | 4 | 4 | 4 | 5 | 4 | 4 | 2 | 6 |
| nb9_applicable | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 5 | 9 |
| nb6_score | 1 | 2 | 2 | 1 | 1 | 2 | 2 | 2 | 2 | 4 |
| nb6_applicable | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 6 |
| nb_has_cys_met | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| parent_id |  |  |  |  |  |  |  |  |  |  |
| gate_strict | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| gate_vm | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| needs_fwdiv | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| vhh_len | 128 | 129 | 127 | 128 | 125 | 120 | 131 | 128 | 124 | 130 |
| len_ok | True | True | True | True | True | True | True | True | True | True |
| invalid_aa | False | False | False | False | False | False | False | False | False | False |
| cys_count | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 |
| noncanon_cys | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| unpaired_cys | False | False | False | False | False | False | False | False | False | False |
| canon_cys_ok | True | True | True | True | True | True | True | True | True | True |
| imm_conf | 0.733 | 0.794 | 0.758 | 0.785 | 0.705 | 0.758 | 0.78 | 0.764 | 0.792 | 0.749 |
| cdrh3_len | 22 | 20 | 22 | 17 | 21 | 13 | 21 | 20 | 17 | 22 |
| cdrh3_dist_trast | 16 | 17 | 16 | 15 | 17 | 10 | 16 | 14 | 15 | 17 |
| cdrh3_dist_trast_imgt | 16 | 17 | 16 | 15 | 17 | 10 | 16 | 14 | 15 | 17 |
| identity_to_trastuzumab | 0.648 | 0.674 | 0.653 | 0.672 | 0.624 | 0.675 | 0.664 | 0.68 | 0.669 | 0.631 |
| cdrh3_min_dist_ref | 11 | 10 | 10 | 7 | 10 | 4 | 10 | 9 | 7 | 10 |
| cdrh3_novel_local | True | True | True | True | True | True | True | True | True | True |
| identity_nearest_known | 0.773 | 0.791 | 0.787 | 0.789 | 0.8 | 0.825 | 0.786 | 0.773 | 0.806 | 0.777 |
| identity_nearest_vhh | 0.773 | 0.791 | 0.787 | 0.789 | 0.8 | 0.825 | 0.786 | 0.773 | 0.806 | 0.777 |
| ident_lt70 | False | False | False | False | False | False | False | False | False | False |
| ident_lt68 | False | False | False | False | False | False | False | False | False | False |
| abguard_pass | True | True | True | True | True | True | True | True | True | True |
| abguard_identity_trast | 0.648 | 0.674 | 0.653 | 0.672 | 0.624 | 0.675 | 0.664 | 0.68 | 0.669 | 0.631 |
| abguard_cdrh3_min_dist | 11 | 10 | 10 | 7 | 10 | 4 | 10 | 9 | 7 | 10 |
| abguard_n_high_liab | 1 | 1 | 2 | 1 | 1 | 1 | 2 | 1 | 1 | 2 |
| basic_pass | True | True | True | True | True | True | True | True | True | True |
| denovo_vm_pass | True | True | True | True | True | True | True | True | True | True |
| denovo_strict_pass | False | False | False | False | False | False | False | False | False | False |
| gap_vm_ok_ident_ge70 | True | True | True | True | True | True | True | True | True | True |
| tier | 1 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | 1 | 2 |
| consensus_rank_score | 0.98 | 0.992 | 0.744 | 0.477 | 0.874 | 0.769 | 0.908 | 0.986 | 0.955 | 0.932 |
| n_models_pass | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| n_models_scored | 1 | 1 | 1 | 2 | 1 | 2 | 1 | 1 | 1 | 3 |
| panel_ident_lt70 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Pareto (last step, `denovo/pareto.py`)

Population 283 gate-passing designs; objectives bind,comp,therm,yield,dev; front sizes {1: 45, 2: 62, 3: 60, 4: 36}; fully measured 187.

JapanFold balance: $83.86
