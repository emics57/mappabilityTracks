#!/bin/bash
#SBATCH --job-name=minimap30x
#SBATCH --partition=long
#SBATCH --mail-user=emxu@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --output=minimap30x.%j.log
#SBATCH --time=168:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=emxu@ucsc.edu
#SBATCH --array=[1-46]%10

# inputs
genome=$1
inputList=$2
aligner=$3

# obtain array task id and chr name (1st and 2nd column)
chr=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' ${inputList})

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID and the chromosome name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the chromosome name is ${chr}"

source /opt/miniconda/etc/profile.d/conda.sh
conda activate cenmap
# if using winnowmap then create kmer file using meryl
if [ ${aligner} == 'winnow' ]; then
    meryl count k=15 output merylDB ${genome}
    meryl print greater-than distinct=0.9998 merylDB > readAlignments/repetitive_k15.txt
fi

# alignment
bamArray=("1000" "5000" "10000" "15000" "20000" "40000" "60000" "80000" "100000" "200000" "300000")
for val in ${bamArray[@]}; do
    fastq=readAlignments/${chr}/fastqs/${val}.fq # fastqPath
    output=readAlignments/${chr}/bams/${val}.bam # outputBamPath

    if [ ${aligner} = 'winnow' ]; then # winnowmap
        echo running winnowmap for ${chr} of size ${val}bp
        winnowmap --split-prefix ${chr}.${val}.tempFile -W readAlignments/repetitive_k15.txt -y -k 15 -t 32 -ax map-ont ${genome} ${fastq} | samtools view -bh - | samtools sort - > ${output}
    else # minimap2
        echo running minimap2 for ${chr} of size ${val}bp
        minimap2 -a -t 32 -x map-ont -I8g ${genome} ${fastq} | samtools view -bh - | samtools sort - > ${output}
    fi
    echo " "
done
conda deactivate

echo "Done."