#!/bin/bash

module load apptainer

# vegan samples
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551665
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551664
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551663

# omnivore samples
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551662
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551661
apptainer exec ../containers/sra-toolkit.sif fasterq-dump SRR10551660

gzip *.fastq
