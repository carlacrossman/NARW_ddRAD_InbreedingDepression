#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-frasiert
#SBATCH --job-name=ref_alignment
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8

bwa index reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta

samtools faidx reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta

samtools dict reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.fasta -o reference/Eubalaena_glacialis_HiC_Mar2023_21scaffs.dict

