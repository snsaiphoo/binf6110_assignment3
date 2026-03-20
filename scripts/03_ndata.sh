#!/bin/bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/sra-toolkit.sif
PIGZ_CONTAINER=/home/ssaiphoo/work/metagenomics/containers/pigz.sif

SCRATCH_DIR=$SCRATCH/metagenomics/raw_data
THREADS=8

failed_samples=()

############################
# VEGAN SAMPLES
############################

vegan_samples=(
    SRR8146952
    SRR8146955
    SRR8146959
    SRR8146960
    SRR8146961
    SRR8146963
    SRR8146965
)

cd $SCRATCH_DIR/vegan

for srr in "${vegan_samples[@]}"
do
    echo "=============================="
    echo "Processing $srr (vegan)"
    echo "=============================="

    TMP_DIR=$SCRATCH/tmp_$srr
    mkdir -p $TMP_DIR

    # Step 1: Convert SRA → FASTQ
    if apptainer exec $CONTAINER fasterq-dump \
        --split-files \
        --threads $THREADS \
        --temp $TMP_DIR \
        $srr
    then
        echo "Conversion successful for $srr"

        # Step 2: Check if FASTQ files exist
        files=(${srr}*.fastq)

        if [ -e "${files[0]}" ]; then
            echo "Compressing FASTQ files for $srr"

            apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}*.fastq

            echo "Compression complete for $srr"
        else
            echo "❌ No FASTQ files found for $srr"
            failed_samples+=("$srr (no FASTQ output)")
        fi

    else
        echo "❌ Conversion FAILED for $srr"
        failed_samples+=("$srr (conversion failed)")
    fi

    # Optional: clean temp
    rm -rf $TMP_DIR

done

############################
# OMNIVORE SAMPLES
############################

omnivore_samples=(
    SRR8146956
    SRR8146957
    SRR8146969
    SRR8146970
    SRR8146971
    SRR8146972
    SRR8146975
)

cd $SCRATCH_DIR/omnivore

for srr in "${omnivore_samples[@]}"
do
    echo "=============================="
    echo "Processing $srr (omnivore)"
    echo "=============================="

    TMP_DIR=$SCRATCH/tmp_$srr
    mkdir -p $TMP_DIR

    if apptainer exec $CONTAINER fasterq-dump \
        --split-files \
        --threads $THREADS \
        --temp $TMP_DIR \
        $srr
    then
        echo "Conversion successful for $srr"

        files=(${srr}*.fastq)

        if [ -e "${files[0]}" ]; then
            echo "Compressing FASTQ files for $srr"

            apptainer exec $PIGZ_CONTAINER pigz -p $THREADS ${srr}*.fastq

            echo "Compression complete for $srr"
        else
            echo "❌ No FASTQ files found for $srr"
            failed_samples+=("$srr (no FASTQ output)")
        fi

    else
        echo "❌ Conversion FAILED for $srr"
        failed_samples+=("$srr (conversion failed)")
    fi

    rm -rf $TMP_DIR

done

############################
# SUMMARY
############################

echo ""
echo "=============================="
echo "PIPELINE SUMMARY"
echo "=============================="

if [ ${#failed_samples[@]} -eq 0 ]; then
    echo "All samples processed successfully ✅"
else
    echo "The following samples had issues:"
    for sample in "${failed_samples[@]}"
    do
        echo " - $sample"
    done
fi
