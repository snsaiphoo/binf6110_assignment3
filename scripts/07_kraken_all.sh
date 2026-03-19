#!/bin/bash

# This script is meant to be run inside an interactive job:
# In the scratch environment, run the following line first and then the script
# srun --cpus-per-task=16 --mem=120G --time=2:00:00 --pty bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken2.sif
BASE=/home/ssaiphoo/scratch/metagenomics
DB=$BASE/kraken2_full
OUT=$BASE/results/kraken2_full_0.15_outputs

mkdir -p $OUT

# vegan samples - ran on the full standard database with a 0.15 confidence threshold
cd $BASE/raw_data/vegan

for SAMPLE in SRR8146944 SRR8146951 SRR8146954
do
  apptainer exec $CONTAINER kraken2 \
    --db $DB \
    --threads 16 \
    --confidence 0.15 \
    --quick \
    --paired \
    --gzip-compressed \
    ${SAMPLE}_1.trimmed.fastq.gz ${SAMPLE}_2.trimmed.fastq.gz \
    --report $OUT/${SAMPLE}_vegan.report \
    --output $OUT/${SAMPLE}_vegan.kraken
done

# omnivore samples - ran on the full standard database with a 0.15 confidence threshold
cd $BASE/raw_data/omnivore

for SAMPLE in SRR8146935 SRR8146936 SRR8146938
do
  apptainer exec $CONTAINER kraken2 \
    --db $DB \
    --threads 16 \
    --confidence 0.15 \
    --quick \
    --paired \
    --gzip-compressed \
    ${SAMPLE}_1.trimmed.fastq.gz ${SAMPLE}_2.trimmed.fastq.gz \
    --report $OUT/${SAMPLE}_omnivore.report \
    --output $OUT/${SAMPLE}_omnivore.kraken
done

