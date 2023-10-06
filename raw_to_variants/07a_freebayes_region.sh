#!/bin/bash
#SBATCH --job-name=freebayes_regions
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=8,10,13,15,21,23,24,27,32,33,34,38,39,40,43,47,49,53,57,59,62,69,73,74,77,90,91,93,104,106,113-373%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=4-00:00:00

REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/ddRAD/vcf/400_regions)
REFERENCE=~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta
LOG_DIR=~/projects/def-frasiert/ddRAD/QC/std_out

freebayes \
        --bam-list bam_rg_list \
        --fasta-reference $REFERENCE \
        --region $REGION \
        --vcf ~/scratch/ddRAD/vcf/freebayes_${REGION}.vcf \
        --report-monomorphic \
        --min-alternate-count 1 2> ${LOG_DIR}/freebayes_${SLURM_ARRAY_TASK_ID}.out

