#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <file.bim>" >&2
  exit 2
fi

bim="$1"
if [ ! -f "$bim" ]; then
  echo "ERROR: BIM file not found: $bim" >&2
  exit 1
fi

awk '
BEGIN {
  total = 0
  chr23 = chr24 = chr25 = chr26 = 0
  autosomal = nonzero_pos = 0
}
{
  total++
  chr[$1]++
  if ($1 == 23 || $1 == "X" || $1 == "x") chr23++
  if ($1 == 24 || $1 == "Y" || $1 == "y") chr24++
  if ($1 == 25 || $1 == "XY" || $1 == "xy") chr25++
  if ($1 == 26 || $1 == "MT" || $1 == "Mt" || $1 == "mt" || $1 == "M") chr26++
  if ($1 ~ /^([1-9]|1[0-9]|2[0-2])$/) autosomal++
  if ($4 != 0) nonzero_pos++
}
END {
  printf("BIM file: %s\n", FILENAME)
  printf("Variants checked: %d\n", total)
  printf("Variants with non-zero positions: %d\n", nonzero_pos)
  printf("Autosomal variants: %d\n", autosomal)
  printf("Sex/MT variants: chr23/X=%d chr24/Y=%d chr25/XY=%d chr26/MT=%d\n", chr23, chr24, chr25, chr26)
  print("Reference build inference: not performed by this self-contained pipeline helper.")
  print("Use the manifest build summary above, or set CHECK_REFBUILD_SH to a project-local reference-build checker if stricter inference is required.")
}
' "$bim"
