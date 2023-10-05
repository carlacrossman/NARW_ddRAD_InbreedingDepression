#!/bin/bash
#SBATCH --job-name=stacks_filter_missing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00


vcftools --gzvcf gt0.01/NARW_stacks0.01_initial.vcf.gz --max-missing 0.6 --out gt0.01/NARW_stacks0.01_missing --recode --recode-INFO-all

gzip gt0.01/NARW_stacks0.01_missing.recode.vcf

vcftools --gzvcf gt0.001/NARW_stacks0.001_initial.vcf.gz --max-missing 0.6 --out gt0.001/NARW_stacks0.001_missing --recode --recode-INFO-all 

gzip gt0.001/NARW_stacks0.001_missing.recode.vcf
