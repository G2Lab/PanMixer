import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_with_custom_command_by_rowid
from tools.utils import load_data_multichromosome

#, "left_inference", "center_inference", "right_inference"]

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
    BASE_PATH,
    PYTHON_ENV
)

def get_command(experiment_number):
    return f"""
if [[ $(hostname) == ne1* ]]; then\
    module load Python/3.10.8-GCCcore-12.2.0
    source {PYTHON_ENV}/bin/activate
    module load bcftools
else
    module load jupyter3
    module load tabix
    module load bcftools
    module load gcc/11.2.0
    module load java
fi

set -e

data_prefix={BASE_PATH}/experiments/exp_{experiment_number}/data/

out=$data_prefix/all/$ROW_ID/merged.vcf.gz

bcftools concat -Oz -o $out \
    $data_prefix/chr1/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr2/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr3/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr4/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr5/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr6/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr7/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr8/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr9/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr10/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr11/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr12/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr13/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr14/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr15/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr16/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr17/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr18/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr19/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr20/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr21/$ROW_ID/new_haplotypes.vcf.gz \
    $data_prefix/chr22/$ROW_ID/new_haplotypes.vcf.gz

bcftools index $out
"""

def combine_vcfs(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    id_counter = 0
    command = get_command(experiment_number)
    launch_job_with_custom_command_by_rowid("combine_vcfs", id_counter, command, memory="32g", cpus="1", num_tasks=str(num_tasks))
