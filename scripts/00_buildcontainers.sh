#!/bin/bash

cd ..
mkdir -p containers
cd containers

module load apptainer

# SRA Toolkit version 3.2.1
singularity pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
mv sra-tools_3.2.1--h4304569_1.sif sra-toolkit.sif

# PIGZ (Parallel compressing)
apptainer pull pigz.sif docker://genevera/docker-pigz

# FASTQC version 0.11.0_cv8
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8
mv fastqc_v0.11.9_cv8.sif fastqc.sif

# MultiQC version 1.33 with PDF support
apptainer pull multiqc.sif docker://multiqc/multiqc:pdf-v1.33
