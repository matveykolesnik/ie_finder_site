#!/bin/bash
# FASTA in, genome-coordinate GFF3 and GenBank out.
# The search itself is MGE_finder/finder_pipeline. This repo only exports attL/attR
# onto the assembly and drops the working files.
#
#   ./run.sh strain.fasta outdir
#   ./run.sh genomes_dir outdir
#   MGE_FINDER=/path/to/MGE_finder ./run.sh strain.fasta outdir
#   ANNOTATE_ALL=1 ./run.sh strain.fasta outdir
#   KEEP_WORK=1 ./run.sh strain.fasta outdir
#
# Extra arguments after the outdir are passed to Snakemake.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ $# -lt 1 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "Usage: ./run.sh genome.fasta|genomes_dir [outdir] [snakemake args...]" >&2
  echo "       MGE_FINDER=/path/to/MGE_finder ./run.sh genome.fasta outdir" >&2
  exit 2
fi

MGE_FINDER="$(realpath "${MGE_FINDER:-$SCRIPT_DIR/../MGE_finder}")"
if [[ ! -f "$MGE_FINDER/finder_pipeline/scripts/predict_orfs.py" ]]; then
  echo "MGE_finder checkout not found at $MGE_FINDER" >&2
  echo "Clone https://github.com/rljech13/MGE_finder.git and set MGE_FINDER." >&2
  exit 1
fi

INPUT="$(realpath "$1")"
shift
if [[ $# -ge 1 && "$1" != -* ]]; then
  OUTDIR="$(realpath -m "$1")"
  shift
else
  OUTDIR="$(realpath -m "$SCRIPT_DIR/outdir")"
fi
mkdir -p "$OUTDIR"

if [[ ! -e "$INPUT" ]]; then
  echo "Input not found: $INPUT" >&2
  exit 1
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ie_finder_site.XXXXXX")"
UPSTREAM_COMMIT=""
if git -C "$MGE_FINDER" rev-parse --short HEAD >/dev/null 2>&1; then
  UPSTREAM_COMMIT="$(git -C "$MGE_FINDER" rev-parse --short HEAD)"
fi

write_reports() {
  shopt -s nullglob
  for sample_dir in "$RESULTS"/*/; do
    sample="$(basename "$sample_dir")"
    if [[ ! -f "$sample_dir/integrase_hits_summary.tsv" && ! -f "$sample_dir/ie_filter_audit.tsv" ]]; then
      continue
    fi
    python3 "$SCRIPT_DIR/scripts/write_report.py" \
      --sample-dir "$sample_dir" \
      --sample "$sample" \
      --mge-finder "$MGE_FINDER" \
      --upstream-commit "$UPSTREAM_COMMIT" \
      --out "$OUTDIR/${sample}.ie.report.txt" \
      || echo "Could not write a report for $sample" >&2
  done
  shopt -u nullglob
}

cleanup() {
  local status=$?
  if [[ $status -ne 0 && -d "${RESULTS:-}" ]]; then
    write_reports || true
  fi
  if [[ $status -eq 0 && "${KEEP_WORK:-0}" != "1" ]]; then
    rm -rf "$WORK"
  elif [[ -d "$WORK" ]]; then
    echo "Work left at $WORK" >&2
  fi
  exit "$status"
}
trap cleanup EXIT

GENOMES="$WORK/genomes"
RESULTS="$WORK/results"
mkdir -p "$GENOMES" "$RESULTS"

link_fasta() {
  local src="$1"
  local base sample dest
  base="$(basename "$src")"
  sample="${base%.*}"
  if [[ -z "$sample" || "$sample" == .* ]]; then
    echo "Could not derive a sample name from $src" >&2
    exit 1
  fi
  dest="$GENOMES/${sample}.fna"
  if [[ -e "$dest" ]]; then
    echo "Two inputs share the sample name $sample" >&2
    exit 1
  fi
  ln -s "$src" "$dest"
}

shopt -s nullglob
if [[ -d "$INPUT" ]]; then
  files=("$INPUT"/*.fna "$INPUT"/*.fa "$INPUT"/*.fasta "$INPUT"/*.fsa)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "No FASTA files in $INPUT" >&2
    exit 1
  fi
  for src in "${files[@]}"; do
    link_fasta "$src"
  done
else
  if [[ ! -s "$INPUT" ]]; then
    echo "FASTA is empty: $INPUT" >&2
    exit 1
  fi
  link_fasta "$INPUT"
fi
shopt -u nullglob

CONFIG="$WORK/config.yaml"
python3 - "$MGE_FINDER/finder_pipeline/ie_finder_config.yaml" "$CONFIG" "$MGE_FINDER" "$GENOMES" "$RESULTS" "${ANNOTATE_ALL:-0}" << 'PY'
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required")

src, dst, mge, genomes, results, annotate_all = sys.argv[1:7]
cfg = yaml.safe_load(Path(src).read_text()) or {}
mge = Path(mge)
cfg["mge_finder"] = str(mge)
cfg.setdefault("paths", {})
cfg["paths"]["config_file"] = dst
cfg["paths"]["genomes_dir"] = genomes
cfg["paths"]["results_dir"] = results
cfg.setdefault("execution", {})
cfg["execution"]["conda_env"] = str(mge / "envs" / "IE_finder.yaml")
cfg["annotate"] = {"all_candidates": annotate_all == "1"}
cfg["pfam_profiles"] = [
    str(mge / "pfam" / "PF00589.hmm"),
    str(mge / "pfam" / "PF22022.hmm"),
]
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

SM="${SNAKEMAKE:-snakemake}"
if ! command -v "$SM" >/dev/null 2>&1; then
  if [[ -x "$HOME/miniforge3/envs/snakemake/bin/snakemake" ]]; then
    SM="$HOME/miniforge3/envs/snakemake/bin/snakemake"
  else
    echo "snakemake is not on PATH" >&2
    exit 1
  fi
fi

CONDA_ARGS=()
if [[ "${USE_CONDA:-0}" == "1" ]]; then
  CONDAP="${SNAKEMAKE_CONDA_PREFIX:-${MINIFORGE_ENVS:-$HOME/miniforge3/envs}}"
  CONDA_ARGS=(--use-conda --conda-prefix "$CONDAP")
fi

"$SM" \
  --snakefile "$SCRIPT_DIR/Snakefile" \
  --configfile "$CONFIG" \
  --cores "$(nproc)" \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
  "${CONDA_ARGS[@]}" \
  "$@"

shopt -s nullglob
found=0
for gff in "$RESULTS"/*/*.ie.gff3; do
  sample="$(basename "$(dirname "$gff")")"
  cp -f "$gff" "$OUTDIR/${sample}.ie.gff3"
  cp -f "$RESULTS/$sample/${sample}.ie.gbk" "$OUTDIR/${sample}.ie.gbk"
  echo "GFF3:    $OUTDIR/${sample}.ie.gff3"
  echo "GenBank: $OUTDIR/${sample}.ie.gbk"
  echo "Report:  $OUTDIR/${sample}.ie.report.txt"
  found=1
done
shopt -u nullglob
write_reports
if [[ $found -eq 0 ]]; then
  echo "Snakemake finished without a GFF3" >&2
  exit 1
fi
