# Failure modes met during the hackathon, and what to do

| Symptom | Cause | Handling |
|---|---|---|
| JapanFold jobs stay `queued` (`account_backlog`, fleet 264/264 chips busy) | shared fleet saturated | it is queueing, not a fault; keep at most ~12 design jobs and ~9 scoring jobs in flight so the stages do not starve each other; never resubmit |
| OpenFold3 interface scores at decoy level (median ipTM about 0.13), `of3_ipsae` exactly 0 | our recipe is harder than upstream's: MSA on, no template (JapanFold and Nebius OpenFold3 take no template); ipSAE with the 10 A PAE cutoff is 0 when best interface PAE > 10 A | treat OpenFold3 as a weak ordering signal; use ipSAE at 15/20 A (`ipsae15`, `ipsae20`) or ipae_tm; do not read it as "designs are bad" or "designs bind"; a templated no-MSA run needs the Apheris OpenFold3 image on a GPU (not done) |
| NanoMelt/TNP `CUDA out of memory` on an 8 GB GPU | ESM models plus other processes | run NanoMelt on CPU (`CUDA_VISIBLE_DEVICES=""`, about 1 s per sequence); TNP can use `--tnp-cpu` |
| `props.csv` shows rows but empty Tm | a tool failed for that batch | columns stay empty (never 0); the `tools_done` column says what ran; rerun with `--retry-failed` |
| abguard class flips between `de_novo` and `novel` for similar designs | identity to trastuzumab VH = 1 - Levenshtein / max(Fv length, 120): the 70% line sits at Fv length about 125 for the BCII10 framework | steer to longer CDR3; read `gate-analysis.md`; the caplacizumab-like framework passes at any length |
| zero-Cys designs pass abguard | zero is an even Cys count, so "no unpaired Cys" passes trivially | canonical-Cys hard gate in `rank.py`/`rules.py` (`canon_cys_ok`) |
| leaderboard stale, `rank.py` crashed with "sort array is read-only" | pandas copy-on-write arrays are read-only | fixed (`np.sort(...)` on a copy); if you see it again check `denovo/rank.py` `pct()` |
| designs measured on one term only topped the composite | weights were renormalised over available terms | missing evidence counts as a 0.3 prior |
| developability objective dropped from Pareto | it required TNP, which covers only the top third | now computed from liabilities for all designs (`liab_n_flags`), TNP is an extra penalty |
| `abguard` VM runs fail with connection resets | the shared VM's sshd throttles bursts | the VM scripts reuse one ssh connection and retry; VM is optional and off by default |
| cancelled generator still shows in tables | its files matched `designs_*.csv` | move them to `parked_*.csv`; `rank.py` and `status.py` also filter known parked generators |
| pkill in a shell script kills your own shell | `pkill -f` matches the calling command line | stop by pid files (`denovo/.run/*.pid`) or `touch denovo/STOP` |
| Two modules named `common` | `optimization/baseline/common.py` vs the pipeline's helpers | the pipeline imports `dn_common` (the file `common.py` is a shim) |
