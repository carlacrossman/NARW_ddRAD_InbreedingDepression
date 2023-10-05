#!/bin/bash
#SBATCH --job-name=concat_dp5
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

bcftools concat -f NARW_filtered_file_list_dp5 -Oz --output NARW_freebayes_filtered_dp5.vcf.gz
