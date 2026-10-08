#!/bin/bash
# FASTA in, genome-coordinate GFF3 and GenBank out.
# This directory is a reduced copy of the MGE_finder discovery steps and runs on its own.
#
#   ./run.sh strain.fasta outdir
#   ./run.sh genomes_dir outdir --cores 8
#   ANNOTATE_ALL=1 ./run.sh strain.fasta outdir
#   KEEP_WORK=1 ./run.sh strain.fasta outdir
#   SEARCH_PARAMS=strict.yaml ./run.sh strain.fasta outdir
#
# The first run creates the conda environment named in envs/IE_finder.yaml
# and puts it on PATH; later runs reuse it. Then scripts/ie_finder.py does the
# work. Arguments are passed on to it (see scripts/ie_finder.py --help).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 ]]; then
  echo "Usage: ./run.sh genome.fasta|genomes_dir [outdir] [--cores N]" >&2
  exit 2
fi
# Checked here too, so a typo does not wait for the environment to build.
if [[ "$1" != -* && ! -e "$1" ]]; then
  echo "Input not found: $1" >&2
  exit 1
fi
if [[ -n "${SEARCH_PARAMS:-}" && ! -f "$SEARCH_PARAMS" ]]; then
  echo "Search parameters not found: $SEARCH_PARAMS" >&2
  exit 1
fi

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
  elif [[ "$(cat "$stamp" 2>/dev/null)" != "$want" || ! -x "$prefix/bin/python3" || ! -x "$prefix/bin/prodigal" || ! -x "$prefix/bin/hmmsearch" || ! -x "$prefix/bin/aragorn" || ! -x "$prefix/bin/blastn" ]]; then
    echo "Updating conda environment $name to match envs/IE_finder.yaml" >&2
    "$conda_bin" env update -p "$prefix" -f "$env_file" --prune
    echo "$want" > "$stamp"
  fi
  if [[ ! -x "$prefix/bin/python3" ]]; then
    echo "Environment $name has no python3 after creation: $prefix" >&2
    exit 1
  fi
  export PATH="$prefix/bin:${PATH:-}"
  hash -r
}

ensure_runtime

flags=()
if [[ -n "${SEARCH_PARAMS:-}" ]]; then
  flags+=(--search-params "$SEARCH_PARAMS")
fi
if [[ "${ANNOTATE_ALL:-0}" == "1" ]]; then
  flags+=(--all-candidates)
fi
if [[ "${KEEP_WORK:-0}" == "1" ]]; then
  flags+=(--keep-work)
fi
exec python3 "$SCRIPT_DIR/scripts/ie_finder.py" "$@" ${flags[@]+"${flags[@]}"}
