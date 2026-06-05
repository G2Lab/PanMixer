#!/bin/bash
#SBATCH --job-name=remove_chm13
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=remove_chm13.out
#SBATCH --error=remove_chm13.err

INPUT="pangenome_no_X.vcf.gz"
OUTPUT="pangenome_no_X_no_chm13.vcf.gz"

bcftools view -s ^chm13 "$INPUT" -Oz -o "$OUTPUT"
bcftools index "$OUTPUT"
