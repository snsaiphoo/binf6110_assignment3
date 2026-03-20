#!/bin/bash

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken2.sif
BASE=/home/ssaiphoo/scratch/metagenomics
DB=$BASE/kraken2_full
OUT=$BASE/results/kraken2_full_0.15_outputs

mkdir -p $OUT

failed_samples=()

############################
# VEGAN SAMPLES
############################

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

cd $BASE/raw_data/vegan

for SAMPLE in "${vegan_samples[@]}"
do
    echo "=============================="
    echo "Running Kraken2: $SAMPLE (vegan)"
    echo "=============================="

    if [[ -f "${SAMPLE}_1.trimmed.fastq.gz" && -f "${SAMPLE}_2.trimmed.fastq.gz" ]]; then

        if apptainer exec $CONTAINER kraken2 \
            --db $DB \
            --threads 16 \
            --confidence 0.15 \
            --quick \
            --paired \
            --gzip-compressed \
            ${SAMPLE}_1.trimmed.fastq.gz ${SAMPLE}_2.trimmed.fastq.gz \
            --report $OUT/${SAMPLE}_vegan.report \
            --output $OUT/${SAMPLE}_vegan.kraken
        then
            echo "✅ Completed $SAMPLE"
        else
            echo "❌ Kraken2 FAILED for $SAMPLE"
            failed_samples+=("$SAMPLE (vegan)")
        fi

    else
        echo "❌ Missing trimmed files for $SAMPLE"
        failed_samples+=("$SAMPLE (vegan - missing input)")
    fi
done

############################
# OMNIVORE SAMPLES
############################

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

cd $BASE/raw_data/omnivore

for SAMPLE in "${omnivore_samples[@]}"
do
    echo "=============================="
    echo "Running Kraken2: $SAMPLE (omnivore)"
    echo "=============================="

    if [[ -f "${SAMPLE}_1.trimmed.fastq.gz" && -f "${SAMPLE}_2.trimmed.fastq.gz" ]]; then

        if apptainer exec $CONTAINER kraken2 \
            --db $DB \
            --threads 16 \
            --confidence 0.15 \
            --quick \
            --paired \
            --gzip-compressed \
            ${SAMPLE}_1.trimmed.fastq.gz ${SAMPLE}_2.trimmed.fastq.gz \
            --report $OUT/${SAMPLE}_omnivore.report \
            --output $OUT/${SAMPLE}_omnivore.kraken
        then
            echo "✅ Completed $SAMPLE"
        else
            echo "❌ Kraken2 FAILED for $SAMPLE"
            failed_samples+=("$SAMPLE (omnivore)")
        fi

    else
        echo "❌ Missing trimmed files for $SAMPLE"
        failed_samples+=("$SAMPLE (omnivore - missing input)")
    fi
done

############################
# SUMMARY
############################

echo ""
echo "=============================="
echo "KRAKEN2 SUMMARY"
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
