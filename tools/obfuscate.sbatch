#!/bin/bash
#SBATCH --job-name=obfuscate
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=obfuscate_chr%a.out
#SBATCH --error=obfuscate_chr%a.err

# Usage:
#   sbatch --export=DATA_DIR=/path/to/starting_data,SUBJECT=HG00438,CAPACITY=0.1,OUTPUT_DIR=/path/to/output,SEED=123 obfuscate.sbatch
#
# Or set defaults below and just run: sbatch obfuscate.sbatch

PYTHON=/users/jb4816/.conda/envs/genomics2/bin/python3
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Defaults (override via --export)
DATA_DIR="${DATA_DIR:-$BASE_DIR/starting_data}"
SUBJECT="${SUBJECT:-HG00438}"
CAPACITY="${CAPACITY:-0.1}"
OUTPUT_DIR="${OUTPUT_DIR:-$BASE_DIR/obfuscation_output}"
SEED="${SEED:-}"

CHR=$SLURM_ARRAY_TASK_ID
TASK_SEED=""
if [[ -n "$SEED" ]]; then
  TASK_SEED=$((SEED + CHR))
fi

echo "=== obfuscate.sbatch ==="
echo "DATA_DIR=$DATA_DIR"
echo "SUBJECT=$SUBJECT"
echo "CAPACITY=$CAPACITY"
echo "OUTPUT_DIR=$OUTPUT_DIR"
echo "CHR=$CHR"
echo "SEED=${TASK_SEED:-unseeded}"
echo "========================"

cd "$BASE_DIR"
export PYTHONPATH="$BASE_DIR:$PYTHONPATH"

if [[ -n "$TASK_SEED" ]]; then
  $PYTHON tools/obfuscate.py "$DATA_DIR" "$SUBJECT" "$CAPACITY" "$CHR" "$OUTPUT_DIR" --seed "$TASK_SEED"
else
  $PYTHON tools/obfuscate.py "$DATA_DIR" "$SUBJECT" "$CAPACITY" "$CHR" "$OUTPUT_DIR"
fi
