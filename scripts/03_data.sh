#!/bin/bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif
PIGZ_CONTAINER=/home/ssaiphoo/work/metagenomics/containers/pigz.sif

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data
THREADS=8

vegan_samples=(SRR8146944 SRR8146951 SRR8146954 SRR8146952 SRR8146955 SRR8146959 SRR8146960 SRR8146961 SRR8146963 SRR8146965)
omnivore_samples=(SRR8146935 SRR8146936 SRR8146938 SRR8146956 SRR8146957 SRR8146969 SRR8146970 SRR8146971 SRR8146972 SRR8146975)

cd $SCRATCH_DIR/vegan
for srr in "${vegan_samples[@]}"
do
    apptainer exec $CONTAINER fasterq-dump --split-files --threads $THREADS --temp $SCRATCH $srr &&
    apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}*.fastq
done

cd $SCRATCH_DIR/omnivore
for srr in "${omnivore_samples[@]}"
do
    apptainer exec $CONTAINER fasterq-dump --split-files --threads $THREADS --temp $SCRATCH $srr &&
    apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}*.fastq
done
