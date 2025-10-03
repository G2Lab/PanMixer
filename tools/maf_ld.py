import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_with_custom_command_multichromosome, launch_job_with_custom_command_by_chromosome
from tools.utils import load_data_multichromosome, get_chromosome_path
import os

from constants import BASE_PATH

def get_command(experiment_number):
    return f"""
module load julia

cd {BASE_PATH}/tools/evaluations_julia
julia evaluations.jl {experiment_number} $CHR $ROW_ID
"""

def maf_ld_by_chrom(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    launch_job_with_custom_command_by_chromosome(
        "maf_ld",
        0,
        get_command(experiment_number),
        memory="32g",
        cpus=1,
        num_tasks=num_tasks,
    )

def maf_ld(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])
    

    launch_job_with_custom_command_multichromosome(
        "maf_ld",
        0,
        get_command(experiment_number),
        memory="32g",
        cpus=1,
        num_tasks=num_tasks,
    )
