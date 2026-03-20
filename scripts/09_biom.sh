#!/bin/bash

BASE=/scratch/ssaiphoo/metagenomics
IN=$BASE/results/kraken2_16gb_0.05_outputs
CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken-biom.sif

cd $IN

module load apptainer

apptainer exec $CONTAINER kraken-biom \
  *_bracken_species.report \
  --fmt json \
  -o $IN/bracken_0.05.biom

# To create BIOM tables for other runs (16GB of different confidence levels),
# just change: the input folder and the output file name
