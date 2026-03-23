#!/bin/bash
#SBATCH --job-name=kraken2_full_0.15
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=/home/ssaiphoo/scratch/metagenomics/logs/kraken2_full_%j.out
#SBATCH --error=/home/ssaiphoo/scratch/metagenomics/logs/kraken2_full_%j.err

module load apptainer

CONTAINER=/home/ssaiphoo/work/metagenomics/containers/kraken2.sif
BASE=/home/ssaiphoo/scratch/metagenomics
DB=$BASE/kraken2_full
OUT=$BASE/results/kraken2_full_0.15_2outputs

mkdir -p $OUT

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

# this for loop ensures the paired-end reads exist for Kraken2 to run 
# reads are processed in paired, compressed format with a confidence threshold of 0.15 and quick mode enabled.
# fail/success is tracked per run 
# process repeated for the omnivore samples after the vegans 

cd $BASE/raw_data/vegan

for SAMPLE in "${vegan_samples[@]}"
do
    echo "$SAMPLE - vegan"

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
            echo "completed $SAMPLE"
        else
            echo "failed $SAMPLE"
            failed_samples+=("$SAMPLE vegan")
        fi

    else
        failed_samples+=("$SAMPLE vegan missing input")
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

cd $BASE/raw_data/omnivore

for SAMPLE in "${omnivore_samples[@]}"
do
    echo "$SAMPLE (omnivore)"
    
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
            echo "completed $SAMPLE"
        else
            echo "failed $SAMPLE"
            failed_samples+=("$SAMPLE omnivore")
        fi

    else
        echo "Missing input: $SAMPLE"
        failed_samples+=("$SAMPLE omnivore - missing input")
    fi
done

# if samples failed print them out 

if [ ${#failed_samples[@]} -eq 0 ]; then
    echo "successful"
else
    echo "samples with issues:"
    for sample in "${failed_samples[@]}"
    do
        echo "$sample"
    done
fi
