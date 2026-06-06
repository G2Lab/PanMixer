#!/bin/bash
#SBATCH --job-name=get_blocks_hg38
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_blocks_hg38%a.out
#SBATCH --error=get_blocks_hg38%a.err

set -euo pipefail

conda activate genomics2
module load plink/1.9 2>/dev/null || true

# Submit this job from starting_data/. Slurm executes a spool
# copy of this file, so resolve the helper from the submission directory.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}/scripts"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

bash "$SCRIPT_DIR/src/get_blocks_grch38.sh" "$SLURM_ARRAY_TASK_ID"
