#!/bin/bash
#SBATCH --job-name=julia_ref_eval
#SBATCH --output=logs/julia_ref_eval_%A_%a.out
#SBATCH --error=logs/julia_ref_eval_%A_%a.err
#SBATCH --array=1-22
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load Julia (adjust if needed for your cluster)
module load julia

# Define chromosome list (1–22, X, Y)
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

# Pick chromosome for this array job
chr=${chromosomes[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "Running chromosome $chr"

# Run the Julia script
julia ../tools/evaluations_julia/evaluations.jl ref $chr