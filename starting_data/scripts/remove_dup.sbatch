#!/bin/bash
#SBATCH --job-name=remove_dup
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=remove_dup_%a.out
#SBATCH --error=remove_dup_%a.err

# Load modules if needed (uncomment if you use modules)

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

module load bcftools

INPUT="${CHR_DIR}/1000g_phased.vcf.gz"
OUTPUT="${CHR_DIR}/1000g_phased_clean.vcf.gz"

bcftools norm -d all -Oz -o $OUTPUT $INPUT
bcftools index $OUTPUT

INPUT="${CHR_DIR}/1000g.vcf.gz"
OUTPUT="${CHR_DIR}/1000g_clean.vcf.gz"

bcftools norm -d all -Oz -o $OUTPUT $INPUT
bcftools index $OUTPUT