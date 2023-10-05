#!/bin/bash
#SBATCH --job-name=stacks_ref
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=24:00:00

ref_map.pl --samples ~/scratch/ddRAD/temp_mapping/ --popmap ~/scratch/ddRAD/stacks_ref/stacks_ref_pop_map --out-path ~/scratch/ddRAD/stacks_ref/gt0.01 -T 8 -X "populations: -r 0.8 --vcf-all" -X "gstacks: --gt-alpha 0.01 --var-alpha 0.01 --details"

mv populations.all.vcf NARW_stacks0.01_initial.vcf

gzip NARW_stacks0.01_initial.vcf

