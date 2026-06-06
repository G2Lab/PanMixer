#!/bin/bash
#SBATCH --job-name=get_total_util
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=get_total_util.out
#SBATCH --error=get_total_util.err

PYTHON="${PYTHON:-python3}"

$PYTHON ./scripts/src/get_total_utility_loss.py
