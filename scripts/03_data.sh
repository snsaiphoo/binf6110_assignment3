#!/bin/bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif
PIGZ_CONTAINER=/home/ssaiphoo/work/metagenomics/containers/pigz.sif

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data
THREADS=8

# converting SRA to FASTQ
# zipping the FASTQ vegan samples

cd $SCRATCH_DIR/vegan

for srr in SRR8146944 SRR8146951 SRR8146954
do
    # Convert SRA to FASTQ
    apptainer exec $CONTAINER fasterq-dump \
        --split-files \
        --threads $THREADS \
        --temp $SCRATCH \
        $srr

    # Using pigz (parallel gzip) for faster compression.
    # See 00_buildcontainers.sh for pigz container setup
    apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}_*.fastq

done

# omnivore samples repeating steps above

cd $SCRATCH_DIR/omnivore

for srr in SRR8146935 SRR8146936 SRR8146938
do
    # Convert SRA → FASTQ
    apptainer exec $CONTAINER fasterq-dump \
        --split-files \
        --threads $THREADS \
        --temp $SCRATCH \
        $srr

    # Compress FASTQ files
    apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}_*.fastq

done
