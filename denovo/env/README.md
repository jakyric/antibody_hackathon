# Optional micromamba environments for the property stage (S4)

`tnp.yml` and `nanomelt.yml` are `micromamba env export --no-builds` of the environments the pipeline was run with (the two tools
themselves were editable installs and are removed from the export; `denovo/setup.sh --tools` clones them at the commits in
`tools.txt` and installs them into the environments). Both run on CPU; NanoMelt predicts about 1 sequence per second on 8 threads.

```bash
bash denovo/setup.sh --tools     # creates envs `tnp` and `nanomelt`, clones tools/TNP and tools/nanomelt, pip install -e
micromamba run -n nanomelt nanomelt predict -i seqs.fasta -align -o tm.csv      # force CPU: CUDA_VISIBLE_DEVICES=""
micromamba run -n tnp TNP --output out_dir --file seqs.fasta --ncores 4
```

Licences: NanoMelt is CC BY-NC-SA 4.0 (non-commercial use only); check TNP's own licence before commercial use. NanoMelt downloads
about 3 GB of ESM weights to the torch hub cache on first use. Without these environments the pipeline still runs: S4 then skips
Tm and TNP (columns stay empty, never 0) and keeps sequence liabilities and the nine-rule score.
