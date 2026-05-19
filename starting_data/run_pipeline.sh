#!/bin/bash
# Pipeline script: submits all sbatch jobs sequentially using SLURM dependencies.
# Each step waits for the previous one to complete successfully before starting.

set -euo pipefail

cd "$(dirname "$0")"

submit_dependent() {
    local script="$1"
    local dep_jobid="${2:-}"

    if [ -n "$dep_jobid" ]; then
        jobid=$(sbatch --dependency=afterok:${dep_jobid} --parsable "$script")
    else
        jobid=$(sbatch --parsable "$script")
    fi

    echo "Submitted $script -> Job ID: $jobid" >&2
    echo "$jobid"
}

echo "=== Starting pipeline ==="

# 1. Remove X chromosome
JOB1=$(submit_dependent scripts/remove_X.sbatch)

# 2. Remove chm13 sample
JOB2=$(submit_dependent scripts/remove_chm13.sbatch "$JOB1")

# 3. Split data by chromosome
JOB3=$(submit_dependent scripts/split_data.sbatch "$JOB2")

# 4. Identify unique alleles and variants
JOB4=$(submit_dependent scripts/get_num_alleles.sbatch "$JOB3")

# 5. Compute LD blocks
JOB5=$(submit_dependent scripts/get_blocks.sbatch "$JOB4")

# 6. Convert VCF files to Numpy
JOB6=$(submit_dependent scripts/convert_2_npy.sbatch "$JOB5")

# 7. Compute variant mappings
JOB7=$(submit_dependent scripts/get_mappings.sbatch "$JOB6")

# 8. Refine segmented blocks
JOB8=$(submit_dependent scripts/segment_blocks.sbatch "$JOB7")

# 9. Compute allele frequencies
JOB9=$(submit_dependent scripts/get_af.sbatch "$JOB8")

# 10. Compute PMI and utility loss
JOB10=$(submit_dependent scripts/get_pmi_utility.sbatch "$JOB9")

# 11. Compute total utility loss JSON
JOB11=$(submit_dependent scripts/get_total_utility.sbatch "$JOB10")

echo "=== All jobs submitted ==="
echo "Final job ID: $JOB11"
