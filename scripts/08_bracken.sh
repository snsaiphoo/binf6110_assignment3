#!/bin/bash

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

KRAKEN_DB=$BASE/kraken2_16db
KRAKEN_OUT=$BASE/results/kraken2_16gb_0.10_outputs
OUT=$BASE/results/bracken_0.10_outputs

READ_LEN=150

mkdir -p "$OUT"

# bracken

for report in "$KRAKEN_OUT"/*_0.10_kraken2.report
do
  base=$(basename "$report" _0.10_kraken2.report)

  apptainer exec "$CONTAINERS/bracken.sif" bracken \
    -d "$KRAKEN_DB" \
    -i "$report" \
    -o "$OUT/${base}_0.10.bracken" \
    -r "$READ_LEN" \
    -l S \
    -t 10
done

