#!/bin/bash

module load apptainer

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data
CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif

failed_samples=()

# validate the vegan samples were prefetched properly

vegan_samples=(
    SRR8146944
    SRR8146951
    SRR8146954
    SRR8146952
    SRR8146955
    SRR8146959
    SRR8146960
    SRR8146961
    SRR8146963
    SRR8146965
)

cd $SCRATCH_DIR/vegan

# for loop iterating through the different samples to validate their download
# same is done for the omnivore samples

for srr in "${vegan_samples[@]}"
do
    if ! apptainer exec $CONTAINER vdb-validate $srr; then
        echo "failed $srr - vegan"
        failed_samples+=("$srr (vegan)")
    else
        echo "passed $srr - vegan"
    fi
done

# validating the omnivore samples
omnivore_samples=(
    SRR8146935
    SRR8146936
    SRR8146938
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
    if ! apptainer exec $CONTAINER vdb-validate $srr; then
        echo "failed  $srr - omnivore"
        failed_samples+=("$srr (omnivore)")
    else
        echo "passed $srr"
    fi
done

if [ ${#failed_samples[@]} -eq 0 ]; then
    echo "samples passed"
else
    echo "failed samples"
    for sample in "${failed_samples[@]}"
    do
        echo "$sample"
    done
fi
