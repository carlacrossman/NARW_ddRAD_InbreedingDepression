#!/bin/bash
#SBATCH --job-name=filter_stacks
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=2-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=18:00:00

THRESHOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/scratch/ddRAD/stacks_thresholds)

cp gt${THRESHOLD}/NARW_stacks${THRESHOLD}_missing.recode.vcf.gz ${SLURM_TMPDIR}
cp stacks${THRESHOLD}_missing_indv ${SLURM_TMPDIR}
cp ~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz ${SLURM_TMPDIR}
cp NARW_good_nodup_stacks${THRESHOLD}_list ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

bcftools view -S ^stacks${THRESHOLD}_missing_indv NARW_stacks${THRESHOLD}_missing.recode.vcf.gz | bcftools filter -e 'F_MISSING>0.2' -Oz -o NARW_stacks${THRESHOLD}_missingness_filtered.recode.vcf.gz

bcftools sort -Oz --output NARW_stacks${THRESHOLD}_filtered_sorted.vcf.gz NARW_stacks${THRESHOLD}_missingness_filtered.recode.vcf.gz

bedtools intersect -v -a NARW_stacks${THRESHOLD}_filtered_sorted.vcf.gz -b Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > NARW_stacks${THRESHOLD}_filtered_sorted_norepeats.vcf.gz

bedtools merge -d 1 -i NARW_stacks${THRESHOLD}_filtered_sorted_norepeats.vcf.gz > NARW_stacks${THRESHOLD}_allsites_called.bed

bcftools view -m2 -M2 -v snps -Oz --output NARW_stacks${THRESHOLD}_biallelic_snps.vcf.gz NARW_stacks${THRESHOLD}_filtered_sorted_norepeats.vcf.gz

bcftools view --samples-file NARW_good_nodup_stacks${THRESHOLD}_list -Oz --output NARW_stacks${THRESHOLD}_biallelic_snps_females_nodups.vcf.gz NARW_stacks${THRESHOLD}_biallelic_snps.vcf.gz

bcftools +fill-tags -Oz --output NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_newaf.vcf.gz NARW_stacks${THRESHOLD}_biallelic_snps_females_nodups.vcf.gz -- -t all

vcftools --gzvcf NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_newaf.vcf.gz  --mac 3 --maf 0.01 --recode --recode-INFO-all --out NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_affiltered

bgzip NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_affiltered.recode.vcf

cp NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_affiltered.recode.vcf.gz ${SLURM_SUBMIT_DIR}
cp NARW_stacks${THRESHOLD}_biallelic_snp_unique_female_newaf.vcf.gz ${SLURM_SUBMIT_DIR}
cp NARW_stacks${THRESHOLD}_allsites_called.bed ${SLURM_SUBMIT_DIR}
cp NARW_stacks${THRESHOLD}_filtered_sorted_norepeats.vcf.gz ${SLURM_SUBMIT_DIR}
cp NARW_stacks${THRESHOLD}_biallelic_snps.vcf.gz ${SLURM_SUBMIT_DIR}
