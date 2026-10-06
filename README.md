# ie_finder_site

A reduced copy of the MGE_finder discovery workflow. It runs on its own: the integrase HMMs, the search scripts, and the conda environment are in this repository. `UPSTREAM` records which commit of [MGE_finder](https://github.com/rljech13/MGE_finder) the copy was taken from. A checkout of that repository is not required.

```bash
./run.sh strain.fasta outdir
./run.sh genomes_dir outdir
```

`outdir/` receives three files per assembly:

- `strain.ie.gff3`
- `strain.ie.gbk`
- `strain.ie.report.txt` — one report: how many integrases and tRNA pairs were found, why a candidate was rejected, the published coordinates, and any step log that was not empty

Coordinates are 1-based on the contigs of the input FASTA. Each element has a span, attL, attR, and the integrase CDS. Intermediate tables and the cut-out islands are written to a temporary directory and then removed. `KEEP_WORK=1` keeps that directory. `ANNOTATE_ALL=1` also writes candidates that have an attL but failed a later filter. If the run fails, the report is still written from whatever finished.

Thresholds are in `ie_finder_config.yaml`. On a fresh machine, `USE_CONDA=1 ./run.sh ...` builds `envs/IE_finder.yaml`.
