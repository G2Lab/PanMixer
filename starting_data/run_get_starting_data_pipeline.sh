#!/bin/bash
# Download and preprocessing pipeline for PanMixer starting data.
# Local download/index steps run first, then SLURM jobs are submitted with dependencies.

set -euo pipefail

cd "$(dirname "$0")"

resolve_panmixer_python() {
    if [ -n "${PYTHON:-}" ]; then
        return
    fi

    local env_prefix=""
    if command -v conda >/dev/null 2>&1; then
        env_prefix=$(conda info --envs | awk '$1 == "panmixer" {print $NF; exit}')
    fi

    if [ -n "$env_prefix" ] && [ -x "$env_prefix/bin/python" ]; then
        PYTHON="$env_prefix/bin/python"
        PYTHONNOUSERSITE="${PYTHONNOUSERSITE:-1}"
    else
        PYTHON="python3"
    fi

    export PYTHON
    export PYTHONNOUSERSITE
}

resolve_panmixer_python

submit_dependent() {
    local script="$1"
    local dep_jobid="${2:-}"

    if [ -n "$dep_jobid" ]; then
        jobid=$(sbatch --export=ALL --dependency=afterok:${dep_jobid} --parsable "$script")
    else
        jobid=$(sbatch --export=ALL --parsable "$script")
    fi

    echo "Submitted $script -> Job ID: $jobid" >&2
    echo "$jobid"
}

echo "=== Starting pipeline ==="

# 1. Download input datasets
echo "=== Downloading pangenome VCF ==="
./scripts/get_pangenomes.sh

echo "=== Downloading Pangenie alignments ==="
./scripts/get_pangenie_alignments.sh

echo "=== Downloading 1000 Genomes phased panel ==="
./scripts/get_1000g_phased.sh

echo "=== Indexing pangenome VCF ==="
bcftools index -f pangenome.vcf.gz

# 2. Remove X chromosome
JOB1=$(submit_dependent scripts/remove_X.sbatch)

# 3. Remove chm13 sample
JOB2=$(submit_dependent scripts/remove_chm13.sbatch "$JOB1")

# 4. Split data by chromosome
JOB3=$(submit_dependent scripts/split_data.sbatch "$JOB2")

# 5. Identify unique alleles and variants
JOB4=$(submit_dependent scripts/get_num_alleles.sbatch "$JOB3")

# 6. Compute LD blocks
JOB5=$(submit_dependent scripts/get_blocks.sbatch "$JOB4")

# 7. Convert VCF files to Numpy
JOB6=$(submit_dependent scripts/convert_2_npy.sbatch "$JOB5")

# 8. Compute variant mappings
JOB7=$(submit_dependent scripts/get_mappings.sbatch "$JOB6")

# 9. Refine segmented blocks
JOB8=$(submit_dependent scripts/segment_blocks.sbatch "$JOB7")

# 10. Compute allele frequencies
JOB9=$(submit_dependent scripts/get_af.sbatch "$JOB8")

# 11. Compute PMI and utility loss
JOB10=$(submit_dependent scripts/get_pmi_utility.sbatch "$JOB9")

# 12. Compute total utility loss JSON
JOB11=$(submit_dependent scripts/get_total_utility.sbatch "$JOB10")

echo "=== All jobs submitted ==="
echo "Final job ID: $JOB11"
