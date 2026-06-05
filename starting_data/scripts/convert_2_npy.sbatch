#!/bin/bash
#SBATCH --job-name=convert_2_npy
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=convert_2_npy%a.out
#SBATCH --error=convert_2_npy%a.err

PYTHON="${PYTHON:-python3}"

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

VCF1="${CHR_DIR}/pangenome.vcf.gz"
VCF2="${CHR_DIR}/PG.vcf.gz"
VCF3="${CHR_DIR}/1000g_phased.vcf.gz"

$PYTHON ../tools/VCFtoNP.py ${VCF1} true
$PYTHON ../tools/VCFtoNP.py ${VCF2} false
$PYTHON ../tools/VCFtoNP.py ${VCF3} true
