#!/bin/bash
#SBATCH --account=def-frasiert
#SBATCH --job-name=prelim_multiqc
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=5
#SBATCH --time=02:00:00

multiqc .
