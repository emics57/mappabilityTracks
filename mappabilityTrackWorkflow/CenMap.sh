#!/bin/bash
#SBATCH --job-name=CenMap
#SBATCH --partition=long
#SBATCH --mail-user=emxu@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=cenMap.%j.log
#SBATCH --time=168:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=emxu@ucsc.edu

# --------------------------------------------------------------------
# cenMap2.sh 
# |__wgsimJobArray.sh --> Job Array: 46 jobs total (runs 15 at a time)
#   |__ 11 (for each read size) sequential wgsim simulations per job 
# |__alignmentJobArray.sh --> Job Array: 46 jobs total (runs 10 at a time)
#   |__ 11 (for each read size) sequential alignments per job
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# OUTPUTS:
# --readAlignment/
#   |__chr1_MATERNAL/
#       |__fastqs/
#       |__alignmentBams/
#   |__ ...
#   |__chrY_PATERNAL/
#       |__fastqs/
#       |__alignmentBams/
#
# --chr_fastas/   # REMOVED
#   |__chr1_MATERNAL.subset.fa
#   |__ ...
#   |__chrY_PATERNAL.subset.fa
#
# --mappabilityTracks/
#   |__chrMap300000bp/
#   |__chrMap200000bp/
#   |__chrMap100000bp/
#   |__chrMap80000bp/
#   |__chrMap60000bp/
#   |__chrMap40000bp/
#   |__chrMap20000bp/
#   |__chrMap15000bp/
#   |__chrMap10000bp/
#   |__chrMap5000bp/
#   |__chrMap1000bp/
#
# --mappability.1000bp.bed  # MAPPABILITY TRACKS
# --mappability.5000bp.bed
# --mappability.10000bp.bed
# --mappability.15000bp.bed
# --mappability.20000bp.bed
# --mappability.40000bp.bed
# --mappability.60000bp.bed
# --mappability.80000bp.bed
# --mappability.100000bp.bed
# --mappability.200000bp.bed
# --mappability.300000bp.bed
# --------------------------------------------------------------------

# run this script in a folder you want readAlignment/, chr_fastas/, and mappabilityTracks to be in

# User Inputs:
genomePath=''
coordinates=''
coverage=''
aligner=''

# sbatch CenMap3.sh -g <fasta path> -b <BED coordinates>

while getopts 'g:b:n:a:' flag; do
  case "${flag}" in
    g) genomePath="${OPTARG}" ;; # [required] full path to genome
    b) coordinates="${OPTARG}" ;; # [optional] default to generating tracks genomewide 
    n) coverage="${OPTARG}" ;; # [optional] default to generating 30x coverage
    a) aligner="${OPTARG}" ;; # [optional] default to minimap2. to use winnowmap, specify '-a winnow'
  esac
done

# if coordinate param (BED file of coordinates) is not given, generate coordinates for entire genome
if [ -z ${coordinates} ]; then
    samtools faidx ${genomePath} > ${genomePath}.fai
    coordinates=${genomePath}.fai.bed
    awk '{print $1"\t0\t"$2}' ${genomePath}.fai > ${coordinates}
fi

# if coverage param is not given, then generate 30x coverage for all read sizes
if [ -z ${coverage} ]; then
    coverage=30 # default 30x
fi

# if -a param is not given (aligner type), then default to minimap2 
if [ -z ${aligner} ]; then
    aligner=mini # default minimap2
fi

# create contig file of chromosomes and region size for the job array
chrArray=chrList.txt
awk 'BEGIN {OFS="\t"; print "ArrayTaskID", "jobName", "regionSize", "start", "end"} $1 != "" {regionSize = $3 - $2; print NR, $1, regionSize, $2, $3}' ${coordinates} > ${chrArray}

# create chr_fastas folder
if [ ! -d "chr_fastas" ]; then
  mkdir chr_fastas
fi

# subset centromere fasta files and create chromosome folders for fastqs and bams
source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/mcechova/.conda/envs/methylation
while IFS=$'\t' read -r chr start end; do 
    samtools faidx ${genomePath} "$chr:$start-$end" > chr_fastas/${chr}.subset.fasta
    echo Creating folder for ${chr}
    if [ ! -d "readAlignments/${chr}" ]; then
        mkdir readAlignments/${chr}
    fi
    if [ ! -d "readAlignments/${chr}/fastqs" ]; then
        mkdir readAlignments/${chr}/fastqs
    fi
    if [ ! -d "readAlignments/${chr}/bams" ]; then
        mkdir readAlignments/${chr}/bams
    fi
done < ${coordinates}
conda deactivate

# run jobArrayScript for wgsim
sbatch --wait wgsim30xCovJobArray.sh ${chrArray} ${coverage}
echo wgsim jobs have completed.

# remove intermediate fasta files
rm -r chr_fastas

# run jobArrayScript using minimap2 to align simulated reads to reference genome
sbatch --wait alignmentJobArray.sh ${genomePath} ${chrArray} ${aligner}
echo alignment jobs have completed.

# create mappability track folders for different read sizes
bamArray=("1000" "5000" "10000" "15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
for readSize in ${bamArray[@]}; do
    if [ ! -d "chrMap${readSize}bp" ]; then
        mkdir chrMap${readSize}bp
        echo created mappability folder for ${readSize}bp
    fi
done

# create mappability tracks for BAM files of different read sizes
source /opt/miniconda/etc/profile.d/conda.sh
conda activate /private/home/emxu/.conda/envs/nucfreq
while IFS=$'\t' read -r chr start end; do
    bamArray=("1000" "5000" "10000" "15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
    for val in ${bamArray[@]}; do
        python3 mappabilityTracks.py -b readAlignments/${chr}/bams/${val}.bam -o chrMap${val}bp/${chr}.mappability.bed
    done
done < ${coordinates}
conda deactivate

# remove intermediate readAlignment folder
# rm -r readAlignment

# merge chr.mappability.bed all together into one BED file for each read size
for readSize in ${bamArray[@]}; do
    cat chrMap${readSize}bp/*.bed > mappability.${readSize}bp.bed
done

echo mappability tracks are done.