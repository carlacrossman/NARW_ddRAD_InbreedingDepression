#!/bin/bash
#SBATCH --job-name=freebayes_regions_complex
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=39,90,91,123,130,144,163,200,281,291,292,311,348
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00

REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/ddRAD/vcf/400_regions)
REFERENCE=~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta
LOG_DIR=~/projects/def-frasiert/ddRAD/QC/std_out

freebayes \
        --bam-list bam_rg_list \
        --fasta-reference $REFERENCE \
        --region $REGION \
        --vcf ~/scratch/ddRAD/vcf/freebayes_${REGION}.vcf \
        --report-monomorphic \
        --min-alternate-count 1 \
	--skip-coverage 10000 \
	--use-best-n-alleles 4 2> ${LOG_DIR}/freebayes_${SLURM_ARRAY_TASK_ID}.out

