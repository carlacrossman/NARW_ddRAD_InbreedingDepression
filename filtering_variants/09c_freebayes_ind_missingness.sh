#!/bin/bash
#SBATCH --job-name=indv_missing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00


vcftools --gzvcf NARW_freebayes_filtered_dp5.vcf.gz --out freebayes_dp5 --missing-indv

vcftools --gzvcf NARW_freebayes_filtered.vcf.gz --out freebayes --missing-indv
