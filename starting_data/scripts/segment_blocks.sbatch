#!/bin/bash
#SBATCH --job-name=seg_blocks
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=seg_blocks%a.out
#SBATCH --error=seg_blocks%a.err

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

python3 ./scripts/src/get_blocks_simple.py $CHR