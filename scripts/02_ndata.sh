#!/bin/bash

module load apptainer

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data
CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif

failed_samples=()

# validate the vegan samples were prefetched properly

vegan_samples=(
    SRR8146952
    SRR8146955
    SRR8146959
    SRR8146960
    SRR8146961
    SRR8146963
    SRR8146965
)

cd $SCRATCH_DIR/vegan

for srr in "${vegan_samples[@]}"
do
    echo "running $srr - vegan"
    if ! apptainer exec $CONTAINER vdb-validate $srr; then
        echo "FAILED $srr - vegan"
        failed_samples+=("$srr (vegan)")
    else
        echo "PASSED $srr - vegan"
    fi
done

# validating the omnivore samples
omnivore_samples=(
    SRR8146956
    SRR8146957
    SRR8146969
    SRR8146970
    SRR8146971
    SRR8146972
    SRR8146975
)

cd $SCRATCH_DIR/omnivore

for srr in "${omnivore_samples[@]}"
do
    echo "running $srr - omnivore"
    
    if ! apptainer exec $CONTAINER vdb-validate $srr; then
        echo "FAILED  $srr - omnivore"
        failed_samples+=("$srr (omnivore)")
    else
        echo "PASSED $srr"
    fi
done

echo "-----------------------------"

if [ ${#failed_samples[@]} -eq 0 ]; then
    echo "All samples passed"
else
    echo "failed samples"
    for sample in "${failed_samples[@]}"
    do
        echo " - $sample"
    done
fi
