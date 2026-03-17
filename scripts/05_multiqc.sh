#!/bin/bash

module load apptainer

BASE_DIR=$HOME/work/metagenomics
CONTAINERS=$BASE_DIR/containers

# create output directory
mkdir -p $BASE_DIR/multiqc

# run MultiQC on both vegan and omnivore FastQC results
apptainer exec $CONTAINERS/multiqc.sif \
    multiqc $BASE_DIR/fastqc_reports/ \
    -o $BASE_DIR/multiqc \
    -n microbiome_multiqc_report.html \
    --title "Vegan vs Omnivore Microbiome QC" \
    -f
