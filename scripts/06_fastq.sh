#!/bin/bash

module load apptainer

BASE_DIR=$SCRATCH/metagenomics
CONTAINERS=$HOME/work/metagenomics/containers

THREADS=8

failed_samples=()

# vegan samples

vegan_samples=(
    SRR8146944
    SRR8146951
    SRR8146954
    SRR8146952
    SRR8146955
    SRR8146959
    SRR8146960
    SRR8146961
    SRR8146963
    SRR8146965
)

cd $BASE_DIR/raw_data/vegan

for srr in "${vegan_samples[@]}"
do
    echo "Trimming $srr - vegan"

    # Check input files exist
    if [[ -f "${srr}_1.fastq.gz" && -f "${srr}_2.fastq.gz" ]]; then

        if apptainer exec $CONTAINERS/fastp.sif fastp \
            -i ${srr}_1.fastq.gz \
            -I ${srr}_2.fastq.gz \
            -o ${srr}_1.trimmed.fastq.gz \
            -O ${srr}_2.trimmed.fastq.gz \
            -q 20 \
            -l 50 \
            -w $THREADS \
            -h ${srr}_fastp.html \
            -j ${srr}_fastp.json
        then
            echo "Completed $srr"
        else
            echo "fastp FAILED for $srr"
            failed_samples+=("$srr (vegan)")
        fi

    else
        echo "Missing FASTQ files for $srr"
        failed_samples+=("$srr (vegan - missing input)")
    fi
done

# omnivore samples

omnivore_samples=(
    SRR8146935
    SRR8146936
    SRR8146938
    SRR8146956
    SRR8146957
    SRR8146969
    SRR8146970
    SRR8146971
    SRR8146972
    SRR8146975
)

cd $BASE_DIR/raw_data/omnivore

for srr in "${omnivore_samples[@]}"
do
    echo "Trimming $srr - omnivore"

    if [[ -f "${srr}_1.fastq.gz" && -f "${srr}_2.fastq.gz" ]]; then

        if apptainer exec $CONTAINERS/fastp.sif fastp \
            -i ${srr}_1.fastq.gz \
            -I ${srr}_2.fastq.gz \
            -o ${srr}_1.trimmed.fastq.gz \
            -O ${srr}_2.trimmed.fastq.gz \
            -q 20 \
            -l 50 \
            -w $THREADS \
            -h ${srr}_fastp.html \
            -j ${srr}_fastp.json
        then
            echo "Completed $srr"
        else
            echo "fastp FAILED for $srr"
            failed_samples+=("$srr (omnivore)")
        fi

    else
        echo "Missing FASTQ files for $srr"
        failed_samples+=("$srr (omnivore - missing input)")
    fi
done

if [ ${#failed_samples[@]} -eq 0 ]; then
    echo "All samples trimmed successfully"
else
    echo "The following samples had issues:"
    for sample in "${failed_samples[@]}"
    do
        echo " - $sample"
    done
fi
