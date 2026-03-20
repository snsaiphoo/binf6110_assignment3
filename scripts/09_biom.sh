#!/bin/bash

BASE=/scratch/ssaiphoo/metagenomics
IN=$BASE/results/kraken2_full_0.15_2outputs
CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken-biom.sif

cd $IN

module load apptainer

apptainer exec $CONTAINER kraken-biom \
  *_bracken_species.report \
  --fmt json \
  -o $IN/table_full_0.15.biom

# Changed format to json because HDF5 was having comptability issues in R
# To create BIOM tables for other runs (16GB of different confidence levels),
# just change the input folder
