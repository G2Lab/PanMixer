import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_with_custom_command_multichromosome
from tools.utils import load_data_multichromosome

#, "left_inference", "center_inference", "right_inference"]

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
    BASE_PATH
)

def get_command(experiment_number):
    return f"""
if [[ $(hostname) == ne1* ]]; then\
    module load Python/3.10.8-GCCcore-12.2.0
    source /gpfs/commons/home/jblindenbach/envs/env-numpy-1-21/bin/activate
    module load all/tabixpp/1.1.2-GCC-12.3.0
    module load Java/11-GCC-12.2.0
else
    module load jupyter3
    module load tabix
    module load bcftools
    module load gcc/11.2.0
    module load java
fi

set -e

beagle_path={BASE_PATH}/downloaded_tools/beagle.27Feb25.75f.jar
vcf_path={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr${{CHR}}/${{ROW_ID}}/new_haplotypes.vcf.gz

master_vcf_no_pangenomes={BASE_PATH}/starting_data/chr${{CHR}}/1000g_snps_no_missing_renamed_norm.vcf.gz

output_path={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr${{CHR}}/${{ROW_ID}}/beagle_refined.vcf.gz

java -Xmx64g -jar $beagle_path gt=$vcf_path ref=$master_vcf_no_pangenomes out=${{output_path}} nthreads=16
"""

def beagle_refinement(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    id_counter = 0
    command = get_command(experiment_number)
    launch_job_with_custom_command_multichromosome("beagle", id_counter, command, memory="64g", cpus="16", num_tasks=str(num_tasks))
