import sys
import json
import numpy as np
import pickle
import os

from hmm import HaplotypeHMM

import yaml
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

base_path = f"{BASE_PATH}/starting_data/"
INFINITY = 1e200

def get_support(chromosome, subject_id):
    #1 compute utility vector

    pangenome = np.load(base_path + f"chr{chromosome}/pangenome.npy")
    not_missing = np.sum(pangenome != -1, axis=(0,2))

    valid_positions = np.where(not_missing > 0)[0]

    utility_loss = np.zeros(len(not_missing))
    utility_loss[valid_positions] = 2 / not_missing[valid_positions]

    maximum_utility_loss = 2 * np.sum(utility_loss)
    print(f"Maximum utility loss: {maximum_utility_loss}")

    #3 find a resampling for each block
    original_haplotypes = pangenome[subject_id]
    new_haplotypes = np.zeros_like(pangenome[subject_id])
    new_haplotypes[:] = -1

    blocks = json.load(open(base_path + f"chr{chromosome}/blocks_dict.json"))
    population_af = np.load(base_path + f"chr{chromosome}/allele_frequencies.npy")

    #save new haplotypes
    new_haplotypes = np.load(f"chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy")

    #4 compute support
    support_blocks = np.zeros((len(blocks), 2), dtype=float)
    for j, key in enumerate(blocks.keys()):
        print(f"Completed {j} out of {len(blocks)} blocks", end="\r")
        if len(blocks[key][1]) == 0:
            continue
        if len(blocks[key][1]) == 1:
            support_blocks[j] = utility_loss[blocks[key][1]]
        else:
            new_haplotypes_block = new_haplotypes[blocks[key][1]]
            original_haplotypes_block = original_haplotypes[blocks[key][1]]
            support_block = utility_loss[blocks[key][1]]

            support_block_stacked = np.column_stack((support_block, support_block))
            
            support_blocks[j] = np.sum(support_block_stacked * (new_haplotypes_block != original_haplotypes_block), axis=0)
    
    #save support blocks

    diff = np.mean(new_haplotypes != original_haplotypes)

    print(np.sum(support_blocks), diff * maximum_utility_loss)
    print(diff)

    #np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/support_blocks.npy", support_blocks)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_mappings.py <chromosome> <subject_id>")
        sys.exit(1)
    
    chromosome = sys.argv[1]
    subject_id = int(sys.argv[2])
    get_support(chromosome, subject_id)