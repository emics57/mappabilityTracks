#!/bin/bash
bamArray=("1000" "5000" "10000" "15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
for readSize in ${bamArray[@]}; do
    if [ ! -d "chrMap${readSize}bp" ]; then
        mkdir chrMap${readSize}bp
    fi
done

source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/emxu/.conda/envs/nucfreq
# create mappability tracks for BAM files of different read sizes
while IFS=$'\t' read -r chr start end; do
    bamArray=("1000" "5000" "10000") 
    #"15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
    for val in ${bamArray[@]}; do
        python3 /private/groups/migalab/emxu/tools/mappabilityTracks.py -b readAlignments/${chr}/bams/${val}.bam -o chrMap${val}bp/${chr}.mappability.bed
    done
done < '/private/groups/migalab/emxu/hg002/hg002v1.0.1.activehor.testrun.bed'
conda deactivate

# merge chr.mappability.bed all together into one for each read size
for readSize in ${bamArray[@]}; do
    cat chrMap${readSize}bp/*.bed > mappability.${readSize}bp.bed
done

echo mappability tracks are done.