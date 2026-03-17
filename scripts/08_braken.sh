#!/bin/bash

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

KRAKEN_DB=$BASE/kraken2
KRAKEN_OUT=$BASE/results/kraken2_outputs
OUT=$BASE/results/bracken_outputs

READ_LEN=150

mkdir -p "$OUT"

# bracken

for report in "$KRAKEN_OUT"/*_kraken2.report
do
  base=$(basename "$report" _kraken2.report)

  apptainer exec "$CONTAINERS/bracken.sif" bracken \
    -d "$KRAKEN_DB" \
    -i "$report" \
    -o "$OUT/${base}.bracken" \
    -r "$READ_LEN" \
    -l S \
    -t 10
done

