#!/bin/bash
#SBATCH --job-name=wgsim30x
#SBATCH --partition=medium
#SBATCH --mail-user=emxu@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=5gb
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --output=wgsim30x.%j.log
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=emxu@ucsc.edu
#SBATCH --array=[1-46]%15

inputList=$1
coverage=$2

# obtain array task id and chr name (1st and 2nd column)
chr=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' ${inputList})
regionSize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' ${inputList})
startCoord=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' ${inputList})
endCoord=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' ${inputList})

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID and the chromosome name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the chromosome name is ${chr}, the region size is ${regionSize}"

# set up list of read sizes, fasta path, fastq path
bamArray=("1000" "5000" "10000" "15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
faPath="chr_fastas/${chr}.subset.fasta" 
fqPath="readAlignments/${chr}/fastqs"

# activate cenmap env 
source /opt/miniconda/etc/profile.d/conda.sh
conda activate cenmap

# simulate fastq files for all read sizes
for val in ${bamArray[@]}; do
    echo "	Running wgsim for ${val} read length"
    numReads=$(echo "${regionSize}*${coverage}/${val}" | bc) # calculate number of reads required to get 30x coverage for the specific read length
    echo "  generating ${numReads} reads"
    wgsim -N ${numReads} -1 ${val} -d 0 -S 11 -e 0 -r 0 ${faPath} ${fqPath}/${val}.fq /dev/null
    # samtools faidx ${genomePath} "$chr:$startCoord-$endCoord" | wgsim -N ${numReads} -1 ${val} -d 0 -S 11 -e 0 -r 0 - ${fqPath}/${val}.fq /dev/null
done
conda deactivate

echo "Done."