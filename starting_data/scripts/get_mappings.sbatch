#!/bin/bash
#SBATCH --job-name=get_mappings
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_mappings%a.out
#SBATCH --error=get_mappings%a.err

plink_path="/gpfs/commons/home/jblindenbach/tools/plink"

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

python3 ./scripts/src/get_mappings.py $CHR