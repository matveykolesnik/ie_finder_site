#!/bin/bash
# FASTA in, genome-coordinate GFF3 and GenBank out.
# This directory is a reduced copy of the MGE_finder discovery steps and runs on its own.
#
#   ./run.sh strain.fasta outdir
#   ./run.sh genomes_dir outdir
#   ANNOTATE_ALL=1 ./run.sh strain.fasta outdir
#   KEEP_WORK=1 ./run.sh strain.fasta outdir
#   SEARCH_PARAMS=strict.yaml ./run.sh strain.fasta outdir
#
# The first run creates the conda environment named in envs/IE_finder.yaml
# and puts it on PATH. Later runs reuse it. Search thresholds are in
# search_params.yaml. SEARCH_PARAMS names a file whose keys override them.
#
# Extra arguments after the outdir are passed to Snakemake.

set -euo pipefail

CALL_DIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Relative paths are taken from the directory where the script was started,
# not from this repository after the cd below.
abspath() {
  case "$1" in
    /*) realpath -m "$1" ;;
    *) realpath -m "$CALL_DIR/$1" ;;
  esac
}

cd "$SCRIPT_DIR"

if [[ $# -lt 1 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "Usage: ./run.sh genome.fasta|genomes_dir [outdir] [snakemake args...]" >&2
  exit 2
fi

INPUT="$(abspath "$1")"
shift
if [[ $# -ge 1 && "$1" != -* ]]; then
  OUTDIR="$(abspath "$1")"
  shift
else
  OUTDIR="$(realpath -m "$SCRIPT_DIR/outdir")"
fi
mkdir -p "$OUTDIR"

if [[ ! -e "$INPUT" ]]; then
  echo "Input not found: $INPUT" >&2
  exit 1
fi

BASE_PARAMS="$SCRIPT_DIR/search_params.yaml"
USER_PARAMS=""
if [[ -n "${SEARCH_PARAMS:-}" ]]; then
  USER_PARAMS="$(abspath "$SEARCH_PARAMS")"
  if [[ ! -f "$USER_PARAMS" ]]; then
    echo "Search parameters not found: $USER_PARAMS" >&2
    exit 1
  fi
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ie_finder_site.XXXXXX")"
PARAMS_USED="$WORK/search_params.yaml"
UPSTREAM_COMMIT="$(awk -F': ' '/^commit:/{print $2; exit}' "$SCRIPT_DIR/UPSTREAM")"

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
      --mge-finder "$SCRIPT_DIR" \
      --upstream-commit "$UPSTREAM_COMMIT" \
      --params "${PARAMS_USED:-}" \
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

find_conda() {
  local candidate
  if [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
    echo "$CONDA_EXE"
    return 0
  fi
  candidate="$(command -v conda 2>/dev/null || true)"
  if [[ -n "$candidate" && -x "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi
  for candidate in \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/mambaforge/bin/conda" \
    "$HOME/miniconda3/bin/conda" \
    "$HOME/anaconda3/bin/conda" \
    "/opt/miniforge3/bin/conda" \
    "/opt/conda/bin/conda"
  do
    if [[ -x "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

install_miniforge() {
  local dest arch url installer
  dest="${MINIFORGE_ROOT:-$HOME/miniforge3}"
  arch="$(uname -m)"
  case "$arch" in
    x86_64|aarch64) ;;
    *)
      echo "conda is not installed, and Miniforge has no build for $arch" >&2
      exit 1
      ;;
  esac
  if [[ -x "$dest/bin/conda" ]]; then
    echo "$dest/bin/conda"
    return 0
  fi
  url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname -s)-${arch}.sh"
  installer="$(mktemp "${TMPDIR:-/tmp}/Miniforge3.XXXXXX.sh")"
  echo "Installing Miniforge to $dest" >&2
  curl -fL --retry 3 -o "$installer" "$url"
  bash "$installer" -b -p "$dest"
  rm -f "$installer"
  echo "$dest/bin/conda"
}

ensure_runtime() {
  local conda_bin base name prefix env_file stamp want
  if ! conda_bin="$(find_conda)"; then
    conda_bin="$(install_miniforge)"
  fi
  env_file="$SCRIPT_DIR/envs/IE_finder.yaml"
  name="$(awk '/^name:/{print $2; exit}' "$env_file")"
  base="$("$conda_bin" info --base)"
  prefix="$base/envs/$name"
  # The environment records the checksum of the file it was built from, so
  # an edit to envs/IE_finder.yaml (a new pin, say) triggers an update.
  stamp="$prefix/.ie_finder_site_env.sha256"
  want="$(sha256sum "$env_file" | cut -d' ' -f1)"
  if [[ ! -d "$prefix" ]]; then
    echo "Building conda environment $name from envs/IE_finder.yaml" >&2
    "$conda_bin" env create -f "$env_file"
    echo "$want" > "$stamp"
  elif [[ "$(cat "$stamp" 2>/dev/null)" != "$want" || ! -x "$prefix/bin/snakemake" || ! -x "$prefix/bin/prodigal" || ! -x "$prefix/bin/hmmsearch" || ! -x "$prefix/bin/aragorn" || ! -x "$prefix/bin/blastn" ]]; then
    echo "Updating conda environment $name to match envs/IE_finder.yaml" >&2
    "$conda_bin" env update -p "$prefix" -f "$env_file" --prune
    echo "$want" > "$stamp"
  fi
  if [[ ! -x "$prefix/bin/snakemake" || ! -x "$prefix/bin/python3" ]]; then
    echo "Environment $name has no snakemake after creation: $prefix" >&2
    exit 1
  fi
  export PATH="$prefix/bin:${PATH:-}"
  hash -r
}

ensure_runtime

# The values this run uses: search_params.yaml with SEARCH_PARAMS on top.
merge_args=(--base "$BASE_PARAMS" --out "$PARAMS_USED")
if [[ -n "$USER_PARAMS" ]]; then
  merge_args+=(--override "$USER_PARAMS")
fi
python3 "$SCRIPT_DIR/scripts/merge_search_params.py" "${merge_args[@]}"

SM="${SNAKEMAKE:-snakemake}"
if ! command -v "$SM" >/dev/null 2>&1; then
  echo "snakemake is not on PATH after building the environment" >&2
  exit 1
fi

PYTHONWARNINGS="${PYTHONWARNINGS:-ignore::FutureWarning}" "$SM" \
  --snakefile "$SCRIPT_DIR/Snakefile" \
  --directory "$WORK" \
  --config \
    genomes_dir="$GENOMES" \
    results_dir="$RESULTS" \
    search_params="$PARAMS_USED" \
    all_candidates="${ANNOTATE_ALL:-0}" \
  --cores "$(nproc)" \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
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
