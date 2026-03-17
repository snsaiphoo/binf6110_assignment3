#!/bin/bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif
SCRATCH_DIR=$SCRATCH/metagenomics/raw_data

mkdir -p $SCRATCH_DIR/vegan
mkdir -p $SCRATCH_DIR/omnivore

# vegan samples
cd $SCRATCH_DIR/vegan

# subject ID VOV54
apptainer exec $CONTAINER prefetch SRR8146944

# subject ID 02BA
apptainer exec $CONTAINER prefetch SRR8146951

# subject ID 18PR
apptainer exec $CONTAINER prefetch SRR8146954

# omnivore samples
cd $SCRATCH_DIR/omnivore

# subject ID VOV08
apptainer exec $CONTAINER prefetch SRR8146935

# subject ID VOV104
apptainer exec $CONTAINER prefetch SRR8146936

# subject ID VOV03
apptainer exec $CONTAINER prefetch SRR8146938

