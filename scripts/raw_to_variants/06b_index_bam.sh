#!/bin/bash
#SBATCH --job-name=rw_index_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=2-120%20
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=00:30:00

BAMFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/ddRAD/temp_mapping/bamfile_list)

samtools index $BAMFILE
