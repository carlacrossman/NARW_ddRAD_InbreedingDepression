#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

bcftools concat -f NARW_filtered_file_list -Oz --output NARW_freebayes_filtered.vcf.gz

