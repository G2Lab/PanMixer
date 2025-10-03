#!/bin/bash
#SBATCH --job-name=apply_mask
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --output=apply_mask%a.out
#SBATCH --error=apply_mask%a.err

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

python3 ./scripts/src/apply_mask.py $CHR