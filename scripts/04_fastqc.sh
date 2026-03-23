#!/bin/bash

module load apptainer

# Define base directories
BASE_DIR=$SCRATCH/metagenomics
CONTAINERS=$HOME/work/metagenomics/containers
REPORTS=$HOME/work/metagenomics

# create directories if not made for the fastqc reports
# one for the vegans and one for the omnivore

mkdir -p $REPORTS/fastqc_reports/vegan
mkdir -p $REPORTS/fastqc_reports/omnivore

# for the vegan samples, 8 threads
apptainer exec $CONTAINERS/fastqc.sif fastqc \
    $BASE_DIR/raw_data/vegan/*.fastq.gz \
    -o $REPORTS/fastqc_reports/vegan/ \
    -t 8

# for the omnivore samples, 8 threads
apptainer exec $CONTAINERS/fastqc.sif fastqc \
    $BASE_DIR/raw_data/omnivore/*.fastq.gz \
    -o $REPORTS/fastqc_reports/omnivore/ \
    -t 8

