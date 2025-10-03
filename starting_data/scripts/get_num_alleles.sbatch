#!/bin/bash
#SBATCH --job-name=get_num_alleles
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_num_alleles%a.out
#SBATCH --error=get_num_alleles%a.err

CHR=$SLURM_ARRAY_TASK_ID

python3 ./scripts/src/get_possible_alleles.py $CHR