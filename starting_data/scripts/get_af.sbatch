#!/bin/bash
#SBATCH --job-name=get_af
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=140G
#SBATCH --time=08:00:00
#SBATCH --output=get_af%a.out
#SBATCH --error=get_af%a.err

PYTHON="${PYTHON:-python3}"

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

$PYTHON ./scripts/src/get_af.py $CHR
