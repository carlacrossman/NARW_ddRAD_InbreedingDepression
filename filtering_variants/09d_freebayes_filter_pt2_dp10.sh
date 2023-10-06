#!/bin/bash
#SBATCH --job-name=freebayes_filter_pt2
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00


cp NARW_freebayes_filtered.vcf.gz ${SLURM_TMPDIR}
cp fb_missing_indv ${SLURM_TMPDIR}
cp ~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}

module load StdEnv/2020 gcc/9.3.0 bcftools/1.16 bedtools/2.30

bcftools view -S ^fb_missing_indv NARW_freebayes_filtered.vcf.gz | bcftools filter -e 'F_MISSING>0.2' -Oz -o NARW_freebayes_missingness_filtered.vcf.gz 

bcftools view -M2 --exclude 'MQM<30 || MQMR<30' -Oz --output NARW_freebayes_max2_mapq.vcf.gz NARW_freebayes_missingness_filtered.vcf.gz 

bedtools intersect -v -a NARW_freebayes_max2_mapq.vcf.gz -b ~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC.repeatmasker.trf.windowmasker.gff.gz -wa -header > NARW_freebayes_max2_mapq_masked.vcf

bcftools query -f '%DP\n' NARW_freebayes_max2_mapq_masked.vcf > NARW_DP.tsv

bgzip NARW_freebayes_max2_mapq_masked.vcf

cp NARW_DP.tsv ${SLURM_SUBMIT_DIR}
cp NARW_freebayes_max2_mapq_masked.vcf.gz ${SLURM_SUBMIT_DIR}
