#!/bin/bash

# move to the raw data folder

cd ..
mkdir -p raw_data
cd raw_data


module load apptainer

mkdir -p vegan
cd vegan

# vegan samples
# subject ID VOV54
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146944

# subject ID 02BA
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146951

# subject ID 18PR
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146954

cd ..
mkdir -p omnivore
cd omnivore

# subject ID VOV08
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146935

# subject ID VOV104
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146936

# subject ID VOV03
apptainer exec ../../containers/sra-toolkit.sif prefetch SRR8146938

