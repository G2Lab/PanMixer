#!/bin/bash
#SBATCH --job-name=get_pmi_weights
#SBATCH --array=1-968
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_pmi_weights_%a.out
#SBATCH --error=get_pmi_weights_s%a.err

plink_path="/gpfs/commons/home/jblindenbach/tools/plink"

# Get the subject ID from the SLURM_ARRAY_TASK_ID
subject_id=$(( ($SLURM_ARRAY_TASK_ID - 1) % 44))

# Get the chromosome number from the SLURM_ARRAY_TASK_ID
CHR=$(( ($SLURM_ARRAY_TASK_ID - 1) / 44 + 1 ))

CHR_DIR="./chr${CHR}"

python3 ./scripts/src/verify_sampled.py $CHR $subject_id
