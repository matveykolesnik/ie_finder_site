#!/bin/bash
# Run the whole pipeline on real assemblies and compare each GFF3 with the
# copy in tests/expected/. The assemblies are not in the repository; pass
# their paths. The sample name is the file name without its last suffix.
#
#   tests/regression.sh T_oshimai.fasta T_thermophilus_HB27.fasta
#   tests/regression.sh --update T_oshimai.fasta   # accept the new output
#
# Extra environment variables (CONDA_EXE, SEARCH_PARAMS, ...) reach run.sh.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPECTED="$ROOT/tests/expected"

update=0
if [[ "${1:-}" == "--update" ]]; then
  update=1
  shift
fi
if [[ $# -eq 0 ]]; then
  echo "Usage: tests/regression.sh [--update] genome.fasta..." >&2
  exit 2
fi

OUT="$(mktemp -d "${TMPDIR:-/tmp}/ie_finder_regression.XXXXXX")"
trap 'rm -rf "$OUT"' EXIT

# One input directory, so a single run covers every genome.
mkdir -p "$OUT/genomes"
for fasta in "$@"; do
  if [[ ! -f "$fasta" ]]; then
    echo "Not found: $fasta" >&2
    exit 1
  fi
  ln -s "$(realpath "$fasta")" "$OUT/genomes/$(basename "$fasta")"
done

if ! "$ROOT/run.sh" "$OUT/genomes" "$OUT/results" > "$OUT/run.log" 2>&1; then
  cat "$OUT/run.log" >&2
  echo "run.sh failed" >&2
  exit 1
fi

failed=0
mkdir -p "$EXPECTED"
for fasta in "$@"; do
  base="$(basename "$fasta")"
  sample="${base%.*}"
  got="$OUT/results/$sample.ie.gff3"
  want="$EXPECTED/$sample.ie.gff3"
  if [[ $update -eq 1 ]]; then
    cp "$got" "$want"
    echo "updated  $sample"
  elif [[ ! -f "$want" ]]; then
    echo "missing  $sample: no $want (run with --update to create it)"
    failed=1
  elif diff -u "$want" "$got"; then
    echo "ok       $sample ($(grep -c 'mobile_genetic_element' "$got") elements)"
  else
    echo "CHANGED  $sample"
    failed=1
  fi
done
exit $failed
