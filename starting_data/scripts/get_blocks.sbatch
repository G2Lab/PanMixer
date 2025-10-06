#!/bin/bash
#SBATCH --job-name=get_blocks
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_blocks%a.out
#SBATCH --error=get_blocks%a.err

module load plink/1.9b_6.21-x86_64

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

PHASED_VCF="${CHR_DIR}/1000g_phased.vcf.gz"
$plink_path --vcf ${PHASED_VCF} --make-bed --out ${CHR_DIR}/1000g_phased --allow-extra-chr --vcf-half-call missing

$plink_path --blocks no-pheno-req --bfile ${CHR_DIR}/1000g_phased --out ${CHR_DIR}/blocks --allow-extra-chr