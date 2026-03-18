#!/bin/bash
#SBATCH --job-name=kraken16_10_all
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/ssaiphoo/metagenomics/logs/kraken_015_all_%j.out
#SBATCH --error=/scratch/ssaiphoo/metagenomics/logs/kraken_015_all_%j.err

module load apptainer

BASE=/scratch/ssaiphoo/metagenomics
CONTAINERS=/home/ssaiphoo/work/metagenomics/containers

KRAKEN_DB=$BASE/kraken2_full
OUT=$BASE/results/kraken2_full_0.15_outputs
LOGDIR=$BASE/logs

THREADS=2

mkdir -p $OUT
mkdir -p $LOGDIR

# Loop through groups and samples
for GROUP in omnivore vegan
do
  for f in $(ls $BASE/raw_data/$GROUP/*_1.trimmed.fastq.gz | head -n 1)
  do
    base=$(basename $f _1.trimmed.fastq.gz)

    echo "Starting $base ($GROUP) at $(date)"

    apptainer exec $CONTAINERS/kraken2.sif kraken2 \
      --db $KRAKEN_DB \
      --threads $THREADS \
      --confidence 0.15 \
      --memory-mapping \
      --paired \
      --gzip-compressed \
      --report $OUT/${base}_${GROUP}_full_0.15_kraken2.report \
      --output $OUT/${base}_${GROUP}_full_0.15_kraken2.output \
      $BASE/raw_data/$GROUP/${base}_1.trimmed.fastq.gz \
      $BASE/raw_data/$GROUP/${base}_2.trimmed.fastq.gz

    echo "Finished $base ($GROUP) at $(date)"

    sleep 60

  done
done

