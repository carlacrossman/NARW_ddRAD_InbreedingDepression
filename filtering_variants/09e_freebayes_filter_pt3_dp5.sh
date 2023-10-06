#!/bin/bash
#SBATCH --job-name=filter_freebayes_pt3_dp5
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=8:00:00

module load StdEnv/2020 gcc/9.3.0 bcftools/1.16 bedtools/2.30

#cp NARW_freebayes_max2_mapq_masked_dp5.vcf.gz ${SLURM_TMPDIR}
cp NARW_freebayes_biallelic_snps_dp5.vcf.gz ${SLURM_TMPDIR}
cp ~/scratch/ddRAD/NARW_good_nodup_sample_list_dp5_jul28 ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

#bcftools filter -e 'INFO/DP>8732' -Oz -o NARW_freebayes_max2_mapq_masked_maxdepth_dp5.vcf.gz  NARW_freebayes_max2_mapq_masked_dp5.vcf.gz

#bcftools +fill-tags -Oz --output NARW_freebayes_filtered_allsites_dp5.vcf.gz NARW_freebayes_max2_mapq_masked_maxdepth_dp5.vcf.gz -- -t all

#bcftools view -m2 --min-ac 1 -v snps -Oz --output NARW_freebayes_biallelic_snps_dp5.vcf.gz NARW_freebayes_filtered_allsites_dp5.vcf.gz

bcftools view --samples-file NARW_good_nodup_sample_list_dp5_jul28 -Oz --output NARW_freebayes_biallelic_snps_females_nodups_dp5.vcf.gz NARW_freebayes_biallelic_snps_dp5.vcf.gz

bcftools +fill-tags -Oz --output NARW_freebayes_biallelic_snp_unique_female_newaf_dp5.vcf.gz NARW_freebayes_biallelic_snps_females_nodups_dp5.vcf.gz -- -t all

module load vcftools/0.1.16

vcftools --gzvcf NARW_freebayes_biallelic_snp_unique_female_newaf_dp5.vcf.gz --mac 3 --maf 0.01 --recode --recode-INFO-all --out NARW_freebayes_biallelic_snps_unique_female_af_filtered_dp5

bgzip NARW_freebayes_biallelic_snps_unique_female_af_filtered_dp5.recode.vcf


#cp NARW_freebayes_filtered_allsites_dp5.vcf.gz ${SLURM_SUBMIT_DIR}
#cp NARW_freebayes_biallelic_snps_dp5.vcf.gz ${SLURM_SUBMIT_DIR}
cp NARW_freebayes_biallelic_snps_unique_female_af_filtered_dp5.recode.vcf.gz ${SLURM_SUBMIT_DIR}
