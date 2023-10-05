#!/bin/bash
#SBATCH --job-name=stacks_missing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

module load vcftools/0.1.16

vcftools --gzvcf gt0.01/NARW_stacks0.01_missing.recode.vcf.gz --out stacks_0.01 --missing-indv

vcftools --gzvcf gt0.001/NARW_stacks0.001_missing.recode.vcf.gz --out stacks_0.001 --missing-indv
