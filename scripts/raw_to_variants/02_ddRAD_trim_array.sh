#!/bin/bash
#SBATCH --account=def-frasiert
#SBATCH --job-name=trim_ddrad
#SBATCH --array=2,48
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00


PROJECTS=~/projects/def-frasiert/ddRAD
SCRATCH=~/scratch/ddRAD

SAMPLE_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${PROJECTS}/raw_data/data_files)

FILE=$(echo $SAMPLE_FILE | awk '{print $1}')
SAMPLE=$(echo $SAMPLE_FILE | awk '{print $2}')

java -Xmx8G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE \
    -threads 12 -phred33 \
    ${PROJECTS}/raw_data/Q017961_ddRAD/${FILE}L001_R1_001.fastq.gz \
    ${PROJECTS}/raw_data/Q017961_ddRAD/${FILE}L001_R2_001.fastq.gz \
    ${SCRATCH}/trimmed_reads/${SAMPLE}-FP.fq.gz \
    ${SCRATCH}/temp/${SAMPLE}-FUP.fq.gz \
    ${SCRATCH}/trimmed_reads/${SAMPLE}-RP.fq.gz \
    ${SCRATCH}/temp/${SAMPLE}-RUP.fq.gz \
    ILLUMINACLIP:${PROJECTS}/misc_files/NexteraPE-PE.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
    AVGQUAL:30 MINLEN:36 2> ${PROJECTS}/QC/trimming/${SAMPLE}.trim.out

fastqc -o ${PROJECTS}/QC/fastQC/ ${SCRATCH}/trimmed_reads/${SAMPLE}-FP.fq.gz
fastqc -o ${PROJECTS}/QC/fastQC/ ${SCRATCH}/trimmed_reads/${SAMPLE}-RP.fq.gz

rm ${SCRATCH}/temp/${SAMPLE}-*.fq.gz
