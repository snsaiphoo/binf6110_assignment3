#!/bin/bash

BASE=/scratch/ssaiphoo/metagenomics
IN=$BASE/results/kraken2_full_0.15_outputs
OUT=$BASE/results
CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken_biom.sif

cd $IN

module load apptainer

apptainer exec $CONTAINER kraken-biom \
  *_bracken_species.report \
  --fmt json \
  -o $OUT/bracken_0.15_full.biom

# To create BIOM tables for other runs (16GB of different confidence levels),
# just change: the input folder and the output file name
