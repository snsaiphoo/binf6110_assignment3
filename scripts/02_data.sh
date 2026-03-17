#!/bin/bash

module load apptainer

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data

# vegan samples validation

cd $SCRATCH_DIR/vegan

for srr in SRR8146944 SRR8146951 SRR8146954
do
    vdb-validate $srr
done

# omnivore sample validation

cd $SCRATCH_DIR/omnivore

for srr in SRR8146935 SRR8146936 SRR8146938
do
    vdb-validate $srr
done

# If validation fails, run:
    # rm -rf $srr
    # rm -rf ~/ncbi/public/sra/${srr}*
    # prefetch $srr --max-size 100G


