#!/bin/bash

module load apptainer

BASE_DIR=$SCRATCH/metagenomics
CONTAINERS=$HOME/work/metagenomics/containers

THREADS=8

# vegan trimming samples

cd $BASE_DIR/raw_data/vegan

for srr in SRR8146944 SRR8146951 SRR8146954
do
    apptainer exec $CONTAINERS/fastp.sif fastp \
        -i ${srr}_1.fastq.gz \
        -I ${srr}_2.fastq.gz \
        -o ${srr}_1.trimmed.fastq.gz \
        -O ${srr}_2.trimmed.fastq.gz \
        -q 20 \
        -l 50 \
        -w $THREADS \
        -h ${srr}_fastp.html \
        -j ${srr}_fastp.json
done

# omnivore trimming samples

cd $BASE_DIR/raw_data/omnivore

for srr in SRR8146935 SRR8146936 SRR8146938
do
    apptainer exec $CONTAINERS/fastp.sif fastp \
        -i ${srr}_1.fastq.gz \
        -I ${srr}_2.fastq.gz \
        -o ${srr}_1.trimmed.fastq.gz \
        -O ${srr}_2.trimmed.fastq.gz \
        -q 20 \
        -l 50 \
        -w $THREADS \
        -h ${srr}_fastp.html \
        -j ${srr}_fastp.json
done
