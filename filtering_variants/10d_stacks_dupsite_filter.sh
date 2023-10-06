#!/bin/bash
#SBATCH --job-name=stacks_dupsite_removal
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1:00:00

module load StdEnv/2020 gcc/9.3.0 bcftools/1.16

bcftools norm -d all NARW_stacks0.001_biallelic_snp_unique_female_affiltered.recode.vcf.gz > NARW_stacks0.001_fully_filtered.vcf
bgzip NARW_stacks0.001_fully_filtered.vcf


bcftools norm -d all NARW_stacks0.01_biallelic_snp_unique_female_affiltered.recode.vcf.gz > NARW_stacks0.01_fully_filtered.vcf
bgzip NARW_stacks0.01_fully_filtered.vcf
