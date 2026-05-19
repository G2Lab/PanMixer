#!/bin/bash
#SBATCH --job-name=obfuscate_v1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=obfuscate_v1_chr%a.out
#SBATCH --error=obfuscate_v1_chr%a.err

PYTHON=/users/jb4816/.conda/envs/genomics2/bin/python3
BASE_DIR=/groups/gg2845_gp/jb4816/PanMixer5

DATA_DIR="$BASE_DIR/starting_data"
SUBJECT="HG00438"
CAPACITY="0.5"
OUTPUT_DIR="$BASE_DIR/obfuscation_output/v1_HG00438_c0.5"
SEED="${SEED:-}"

CHR=$SLURM_ARRAY_TASK_ID
TASK_SEED=""
if [[ -n "$SEED" ]]; then
  TASK_SEED=$((SEED + CHR))
fi

cd "$BASE_DIR"
if [[ -n "$TASK_SEED" ]]; then
  $PYTHON tools/obfuscate.py "$DATA_DIR" "$SUBJECT" "$CAPACITY" "$CHR" "$OUTPUT_DIR" --seed "$TASK_SEED"
else
  $PYTHON tools/obfuscate.py "$DATA_DIR" "$SUBJECT" "$CAPACITY" "$CHR" "$OUTPUT_DIR"
fi
