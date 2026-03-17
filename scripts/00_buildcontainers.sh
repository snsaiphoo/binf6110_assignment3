#!/bin/bash

cd ..
mkdir -p containers
cd containers

module load apptainer

# SRA Toolkit version 3.2.1
singularity pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
mv sra-tools_3.2.1--h4304569_1.sif sra-toolkit.sif


