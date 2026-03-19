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

# FASTP version 0.23.4
apptainer pull fastp.sif docker://quay.io/biocontainers/fastp:0.23.4--h5f740d0_0

# Apptainer Kraken2 version 2.17.1
apptainer pull kraken2_2.17.1.sif docker://staphb/kraken2:2.17.1
mv kraken2_2.17.1.sif kraken2.sif

# Bracken version 3.1
apptainer pull bracken_3.1.sif docker://staphb/bracken:3.1
mv bracken_3.1.sif bracken.sif

# Kraken-Biom version 1.0.1
apptainer pull kraken-biom.sif docker://shaunchuah/kraken_biom:v1.0.1
