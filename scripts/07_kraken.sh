#!/bin/bash

module load apptainer

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

KRAKEN_DB=$BASE/kraken2_16db
OUT=$BASE/results/kraken2_16gb_outputs

THREADS=2

mkdir -p $OUT

# kraken2 for both vegan and omnivore

for GROUP in omnivore vegan
do
  for f in $BASE/raw_data/$GROUP/*_1.trimmed.fastq.gz
  do
    base=$(basename $f _1.trimmed.fastq.gz)

    apptainer exec $CONTAINERS/kraken2.sif kraken2 \
      --db $KRAKEN_DB \
      --threads $THREADS \
      --confidence 0.05 \
      --paired \
      --gzip-compressed \
      --report $OUT/${base}_${GROUP}_16gb_kraken2.report \
      --output $OUT/${base}_${GROUP}_16gb_kraken2.output \
      $BASE/raw_data/$GROUP/${base}_1.trimmed.fastq.gz \
      $BASE/raw_data/$GROUP/${base}_2.trimmed.fastq.gz
      sleep 20
  done
done

