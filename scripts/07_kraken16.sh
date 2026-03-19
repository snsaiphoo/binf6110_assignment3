#!/bin/bash
#SBATCH --job-name=kraken16
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --output=kraken_%j.out
#SBATCH --error=kraken_%j.err

module load apptainer

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

KRAKEN_DB=$BASE/kraken2_16db

# changed for the 0.10 confidence threshold
OUT=$BASE/results/kraken2_16gb_0.05_outputs

THREADS=2

mkdir -p $OUT

cd $BASE/raw_data

for GROUP in omnivore vegan
do
  cd $GROUP

  for R1 in *_1.trimmed.fastq.gz
  do
    BASE_NAME=${R1%%_1.trimmed.fastq.gz}

    apptainer exec $CONTAINERS/kraken2.sif kraken2 \
      --db $KRAKEN_DB \
      --threads $THREADS \
      --confidence 0.10 \
      --paired \
      --gzip-compressed \
      --report $OUT/${BASE_NAME}_${GROUP}_16gb_kraken2_0.05.report \
      --output $OUT/${BASE_NAME}_${GROUP}_16gb_kraken2_0.05.output \
      $BASE/raw_data/$GROUP/${BASE_NAME}_1.trimmed.fastq.gz \
      $BASE/raw_data/$GROUP/${BASE_NAME}_2.trimmed.fastq.gz

  done
  cd ..

done
