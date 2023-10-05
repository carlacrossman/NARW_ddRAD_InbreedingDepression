#!/bin/bash
#SBATCH --job-name=filter_freebayes_pt1
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-373%40
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --time=3:00:00

module load vcftools/0.1.16

PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/ddRAD/vcf/freebayes_regions_list)

## FILTER NARW ON DEPTH AND REMOVE SITES WITH MISSING GENOTYPES IN >20% OF SAMPLES

vcftools --gzvcf freebayes_${PREFIX}.vcf.gz --keep NARW_list_filters  --minDP 5 --recode --recode-INFO-all --out ${PREFIX}_pt1_filtered_dp5

vcftools --vcf ${PREFIX}_pt1_filtered_dp5.recode.vcf --max-missing 0.6 --out ${PREFIX}_filtered_dp5 --recode --recode-INFO-all

gzip ${PREFIX}_filtered_dp5.recode.vcf
rm ${PREFIX}_pt1_filtered_dp5.recode.vcf
mv ${PREFIX}_filtered_dp5.recode.vcf.gz filtered_dp5_jul24/

