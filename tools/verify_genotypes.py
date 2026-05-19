import sys
import numpy as np

from tools.utils import load_data_multichromosome 
from tools.slurm_helper import launch_job_multichromosome

from constants import (
    EXPERIMENT_PATH,
    STARTING_DATA_PATH,
)

def verify_all(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    print(num_tasks)

    args = [experiment_number]

    launch_job_multichromosome("verify_genotypes", args, memory="16g", cpus="1", num_tasks=str(num_tasks))

def verify_sample(experiment_number, chromosome, row_id):
    subject_haplotypes = np.load(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")
    possible_alleles = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/num_alleles.npy")

    label = f"{experiment_number}_{chromosome}_{row_id}"

    assert subject_haplotypes.shape[1] == 2, f"{label}: Subject does not have 2 haplotypes"
    assert subject_haplotypes.shape[0] == possible_alleles.shape[0], f"{label}: Subject does not have the same number of haplotypes as possible alleles"

    assert np.all(subject_haplotypes[:, 0] < possible_alleles), f"{label}: Subject has haplotypes that are not in the possible alleles {possible_alleles}"
    assert np.all(subject_haplotypes[:, 1] < possible_alleles), f"{label}: Subject has haplotypes that are not in the possible alleles {possible_alleles}"

if __name__ == "__main__":
    experiment_number = int(sys.argv[1])
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    verify_sample(experiment_number, chromosome, row_id)
