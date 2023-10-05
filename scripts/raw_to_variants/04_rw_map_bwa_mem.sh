#!/bin/bash
#SBATCH --job-name=rw_map_bwa_mem
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=2,32,48,88
#SBATCH --cpus-per-task=20
#SBATCH --mem=15G
#SBATCH --time=04:00:00

INPUT_DIR=~/scratch/ddRAD/trimmed_reads
REFERENCE_GENOME=~/projects/def-frasiert/ddRAD/reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta

LOG_DIR=~/projects/def-frasiert/ddRAD/QC/std_out

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/projects/def-frasiert/ddRAD/sample_list)

cp ${INPUT_DIR}/$SAMPLE-FP.fq.gz ${SLURM_TMPDIR}
cp ${INPUT_DIR}/$SAMPLE-RP.fq.gz ${SLURM_TMPDIR}

cd ${SLURM_TMPDIR}


bwa mem -M -t 20 \
  ${REFERENCE_GENOME} \
  $SAMPLE-FP.fq.gz \
  $SAMPLE-RP.fq.gz > \
  $SAMPLE-aln_21scaffs.sam \
  2> ${LOG_DIR}/$SAMPLE-bwa.err


samtools sort -o ${SAMPLE}_bwa.bam $SAMPLE-aln_21scaffs.sam

cp ${SAMPLE}_bwa.bam ${SLURM_SUBMIT_DIR}
