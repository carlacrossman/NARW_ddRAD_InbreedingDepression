#!/bin/bash
#SBATCH --job-name=new_headers_pt1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=3,30,36,78,85,114,119
#SBATCH --cpus-per-task=6
#SBATCH --mem=25G
#SBATCH --time=1:00:00

SAMPLE_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" reheader_file)

FILE=$(echo $SAMPLE_FILE | awk '{print $2}')
SAMPLE=$(echo $SAMPLE_FILE | awk '{print $1}')

READGROUP="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:RW\tCN:MCGILL\tPL:ILLUMINA"

samtools addreplacerg  -r $READGROUP -o ${SAMPLE}_bwa_rg.bam $FILE

samtools index ${SAMPLE}_bwa_rg.bam
