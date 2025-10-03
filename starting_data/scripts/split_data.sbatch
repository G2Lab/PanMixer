#!/bin/bash
#SBATCH --job-name=split_chr
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=split_chr_%a.out
#SBATCH --error=split_chr_%a.err

# Load modules if needed (uncomment if you use modules)
module load bcftools

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR=chr${CHR}

# Create the output directory for this chromosome
mkdir -p ${CHR_DIR}

# Input files
VCF1="pangenome_no_X.vcf.gz"
VCF2="1000g.vcf.gz"

# Output files
OUT1="${CHR_DIR}/pangenome.vcf.gz"
OUT2="${CHR_DIR}/1000g.vcf.gz"

# Split using bcftools
bcftools view --regions grch38#chr${CHR} ${VCF1} -Oz -o ${OUT1}
bcftools index ${OUT1}

bcftools view --regions chr${CHR} ${VCF2} -Oz -o ${OUT2}
bcftools index ${OUT2}