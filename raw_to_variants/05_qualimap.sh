#!/bin/bash
#SBATCH --job-name=qualimap_testing
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=12:00:00

OUT_DIR=~/projects/def-frasiert/ddRAD/QC/mapping/all/

qualimap multi-bamqc -d qualimap_all_input \
        -outdir ${OUT_DIR} -nr 10000 \
        --java-mem-size=8G -r
