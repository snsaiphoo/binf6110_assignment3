#!/bin/bash

module load apptainer

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

# Standard database used on 10 samples of both vegan and omnivore

KRAKEN_DB=$BASE/kraken2_full

# standard DB, 0.15 confidence
KRAKEN_OUT=$BASE/results/kraken2_full_0.15_2outputs

# output folder
OUT=$BASE/results/bracken_0.15_full_2outputs

READ_LEN=150

mkdir -p "$OUT"

# Bracken on both omnivore and vegan samples
for report in "$KRAKEN_OUT"/*.report
do
  base=$(basename "$report" .report)

  apptainer exec "$CONTAINERS/bracken.sif" bracken \
    -d "$KRAKEN_DB" \
    -i "$report" \
    -o "$OUT/${base}_0.15_full.bracken" \
    -r "$READ_LEN" \
    -l S \
    -t 10
done
