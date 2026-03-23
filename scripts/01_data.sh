#!/bin/bash
# Using the SRA Toolkit to prefetch the SRR samples
# 10 vegan and 10 omnivores

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif
SCRATCH_DIR=$SCRATCH/metagenomics/raw_data

mkdir -p $SCRATCH_DIR/vegan
mkdir -p $SCRATCH_DIR/omnivore

# vegan samples

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

for sample in "${vegan_samples[@]}"
do
    apptainer exec $CONTAINER prefetch $sample
done

# omnivore samples

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

for sample in "${omnivore_samples[@]}"
do
    apptainer exec $CONTAINER prefetch $sample
done
