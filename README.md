# ie_finder_site

A thin front end for putting integrative elements on a genome browser. The search is not copied here. It runs in [MGE_finder](https://github.com/rljech13/MGE_finder), at the commit recorded in `UPSTREAM`. This repository only places attL and attR on the input assembly.

Put a checkout of the main repository next to this directory, or pass the path:

```bash
./run.sh strain.fasta outdir
MGE_FINDER=/path/to/MGE_finder ./run.sh genomes_dir outdir
```

`outdir/` receives three files per assembly:

- `strain.ie.gff3`
- `strain.ie.gbk`
- `strain.ie.report.txt` — one report: how many integrases and tRNA pairs were found, why a candidate was rejected, the published coordinates, and any step log that was not empty

Coordinates are 1-based on the contigs of the input FASTA. Each element has a span, attL, attR, and the integrase CDS. Intermediate tables and the cut-out islands are written to a temporary directory and then removed. `KEEP_WORK=1` keeps that directory. `ANNOTATE_ALL=1` also writes candidates that have an attL but failed a later filter. If the run fails, the report is still written from whatever finished.

Thresholds come from `MGE_finder/finder_pipeline/ie_finder_config.yaml`. The environment file is `MGE_finder/envs/IE_finder.yaml`. On a fresh machine, run `USE_CONDA=1 ./run.sh ...`.
