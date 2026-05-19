import numpy as np
import os
import gzip
import sys

from tools.VCFtoNP import get_numpy_matrices
from constants import EXPERIMENT_PATH
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data_multichromosome

def VCFtoNP_parallel(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    print(num_tasks)

    args = [experiment_number]

    launch_job_multichromosome("VCFtoNP_parallel", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

if __name__ == "__main__":
    experiment_number = int(sys.argv[1])
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])


    vcf_file_path = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.vcf.gz"
    phased = True
    genotypes, subjects, positions = get_numpy_matrices(vcf_file_path, phased)

    # Define output file names, handling both `.vcf` and `.vcf.gz`
    base_name = vcf_file_path.replace('.vcf.gz', '').replace('.vcf', '')

    np.save(f"{base_name}.npy", genotypes)
    np.save(f"{base_name}_subjects.npy", subjects)
    np.save(f"{base_name}_positions.npy", positions)

    print("Done")
