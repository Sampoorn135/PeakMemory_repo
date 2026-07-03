#!/usr/bin/env bash
# =============================================================================
# Annotate the final StrictMemory peaks (both datasets, all 4 stringent
# conditions) to mm10 genes using HOMER's annotatePeaks.pl.
#
# RUN THIS ON YOUR MAC (HOMER is installed there, not in Claude's sandbox):
#     bash "/Users/sampoorn/Claude/Projects/Memory_peak/gene_annotation/run_homer_annotation.sh"
#
# It reads the 8 HOMER-format peak files in gene_annotation/peaks/ and writes
# annotated tables to gene_annotation/homer_out/. Once it finishes, tell Claude
# and it will parse the output and compute the gene-level overlap.
# =============================================================================
set -uo pipefail

BASE="/Users/sampoorn/Claude/Projects/Memory_peak/gene_annotation"
PEAKS="$BASE/peaks"
OUT="$BASE/homer_out"
HOMER_BIN="/Users/sampoorn/homer/bin"
GENOME="mm10"

mkdir -p "$OUT"
export PATH="$HOMER_BIN:$PATH"

# Locate annotatePeaks.pl (use the path you confirmed, fall back to PATH).
ANNOT="$HOMER_BIN/annotatePeaks.pl"
if [ ! -x "$ANNOT" ]; then
  ANNOT="$(command -v annotatePeaks.pl || true)"
fi
if [ -z "${ANNOT:-}" ] || [ ! -x "$ANNOT" ]; then
  echo "ERROR: annotatePeaks.pl not found. Edit HOMER_BIN at the top of this script." >&2
  exit 1
fi
echo "Using: $ANNOT"
echo "Genome: $GENOME"
echo "Output: $OUT"
echo

n=0
for pk in "$PEAKS"/*.homer.txt; do
  name="$(basename "$pk" .homer.txt)"
  echo ">>> Annotating $name ..."
  "$ANNOT" "$pk" "$GENOME" \
      -annStats "$OUT/${name}.annStats.txt" \
      > "$OUT/${name}.annotated.txt" \
      2> "$OUT/${name}.log"
  if [ $? -eq 0 ] && [ -s "$OUT/${name}.annotated.txt" ]; then
      rows=$(( $(wc -l < "$OUT/${name}.annotated.txt") - 1 ))
      echo "    OK -> ${name}.annotated.txt  (${rows} annotated peaks)"
      n=$((n+1))
  else
      echo "    FAILED for $name -- see ${name}.log" >&2
  fi
done

echo
echo "Done: $n / 8 files annotated. Outputs in: $OUT"
echo "If genome mm10 is not configured, run once:  perl $HOMER_BIN/../configureHomer.pl -install mm10"
