import sys
import json
import numpy as np
import sys


import yaml
import os
def load_config():
    # Load default template
    with open("../config.yaml") as f:
        config = yaml.safe_load(f)

    # If user has a local config, override defaults
    if os.path.exists("../config.local.yaml"):
        with open("../config.local.yaml") as f:
            local_config = yaml.safe_load(f)
        config.update(local_config)

    return config

CONFIG = load_config()
BASE_PATH = CONFIG["base_path"]



def verify_sample(chromosome, subject_id):
    subject_haplotypes = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy")

    possible_alleles = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/num_alleles.npy")

    assert subject_haplotypes.shape[1] == 2, f"Subject {subject_id} does not have 2 haplotypes"
    assert subject_haplotypes.shape[0] == possible_alleles.shape[0], f"Subject {subject_id} does not have the same number of haplotypes as possible alleles"

    assert np.all(subject_haplotypes[:, 0] < possible_alleles), f"Subject {subject_id} has haplotypes that are not in the possible alleles {possible_alleles}"
    assert np.all(subject_haplotypes[:, 1] < possible_alleles), f"Subject {subject_id} has haplotypes that are not in the possible alleles {possible_alleles}"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python verify_sampled.py <chromosome> <subject_id>")
        sys.exit(1)

    chromosome = int(sys.argv[1])
    subject_id = int(sys.argv[2])
    verify_sample(chromosome, subject_id)
    print(f"Subject {subject_id} has valid haplotypes")