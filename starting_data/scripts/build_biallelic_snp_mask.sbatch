#!/bin/bash
#SBATCH --job-name=build_biallelic_snp_mask
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=build_biallelic_snp_mask%a.out
#SBATCH --error=build_biallelic_snp_mask%a.err

PYTHON="${PYTHON:-python3}"

CHR=$SLURM_ARRAY_TASK_ID

$PYTHON ./scripts/src/build_biallelic_snp_mask.py $CHR
