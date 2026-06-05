#!/bin/bash
#SBATCH --job-name=get_pmi_weights
#SBATCH --array=1-968
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=get_pmi_weights_%a.out
#SBATCH --error=get_pmi_weights_s%a.err

PYTHON="${PYTHON:-python3}"

# Get the subject ID from the SLURM_ARRAY_TASK_ID
subject_id=$(( ($SLURM_ARRAY_TASK_ID - 1) % 44))

# Get the chromosome number from the SLURM_ARRAY_TASK_ID
CHR=$(( ($SLURM_ARRAY_TASK_ID - 1) / 44 + 1 ))
SEED="${SEED:-}"
TASK_SEED=""
if [[ -n "$SEED" ]]; then
  TASK_SEED=$((SEED + SLURM_ARRAY_TASK_ID))
fi

CHR_DIR="./chr${CHR}"

if [[ -n "$TASK_SEED" ]]; then
  $PYTHON ./scripts/src/get_support_and_pmi.py $CHR $subject_id --seed "$TASK_SEED"
else
  $PYTHON ./scripts/src/get_support_and_pmi.py $CHR $subject_id
fi
