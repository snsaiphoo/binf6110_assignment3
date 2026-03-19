#!/bin/bash

# Both the standard and the standard 16gb databases will be tested
# at different confidence thresholds, so both can be downloaded through this script
# future updates can be found here: https://benlangmead.github.io/aws-indexes/k2

BASE=/scratch/ssaiphoo/metagenomics

# Create folders
mkdir -p $BASE/kraken2_16gb
mkdir -p $BASE/kraken2_full

echo "Starting downloads at $(date)"

# 16GB database
cd $BASE/kraken2_16gb

wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260226.tar.gz
tar -xvzf k2_standard_16_GB_20260226.tar.gz

# standard database
cd $BASE/kraken2_full

wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260226.tar.gz
tar -xvzf k2_standard_20260226.tar.gz
