import numpy as np
from tools.utils import load_data, load_data_multichromosome
from tools.slurm_helper import launch_job_multichromosome
import pickle
import sys
import json

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
)

STACKER_STRATEGIES = [
    "to_best",
    "to_empty",
    "to_random",
    "to_unedited",
]

def stacker(strategy, experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    print(num_tasks)

    args = [strategy, experiment_number]

    launch_job_multichromosome("stacker", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

def main():
    strategy = sys.argv[1]
    experiment_number = sys.argv[2]
    chromosome = int(sys.argv[3])
    row_id = int(sys.argv[4])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]
    capacity = data["capacity"]

    pangenome_subject = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/pangenome_subjects.npy")
    subject_id = np.where(pangenome_subject == subject_name)[0][0]

    blocks = json.load(open(f"{STARTING_DATA_PATH}/chr{chromosome}/blocks_dict.json", "r"))
    keys = list(blocks.keys())

    pangenome_hapotypes = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/pangenome.npy")
    original_haplotypes = pangenome_hapotypes[subject_id]

    new_variants = original_haplotypes.copy()

    if strategy == "to_best":
        xsol = np.load(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/xsol.npy")
        obfuscated_variants = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy")

        print(xsol.shape, obfuscated_variants.shape, len(blocks))

        true_indices = np.argwhere(xsol == 1)
        singletons = 0

        alleles_touched = 0

        for i,j in true_indices:
            block_indices = blocks[keys[i]][1]
            for block_index in block_indices:
                new_variants[block_index][j] = obfuscated_variants[block_index][j]
                alleles_touched += 1


        diff = np.mean(new_variants != original_haplotypes)
        mask = original_haplotypes != -1

        print(singletons)

        print("Mask entries:", np.sum(mask))
        print("Xsol entries:", np.sum(xsol))

        changed_alleles = np.sum(new_variants[mask] != original_haplotypes[mask])
        changed_aleles_percentage = changed_alleles / np.sum(mask)
        print(f"Mean percentage of changed alleles: {changed_aleles_percentage}")
        print(f"Mean difference between new and original haplotypes: {diff}")

        alleles_touched_percentage = alleles_touched / obfuscated_variants.size

        results = {
            "stacker_strategy": strategy,
            "changed_alleles_percentage": 100 * float(changed_aleles_percentage),
            "mean_difference": float(diff),
            "number_of_changed_alleles": int(changed_alleles),
            "number_of_alleles_touched": int(alleles_touched),
            "number_of_alleles_touched_percentage": 100 * float(alleles_touched_percentage),
        }
        print(results)
        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy", new_variants)
        print("Saved new haplotypes")
        with open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/stacked.json", "w") as f:
            json.dump(results, f, indent=4)
        print("Wrote results to json")
    elif strategy == "to_empty":
        empty_haplotypes = np.zeros_like(original_haplotypes)
        empty_haplotypes.fill(-1)

        new_variants = empty_haplotypes.copy()

        diff = np.mean(new_variants != original_haplotypes)
        mask = original_haplotypes != -1

        changed_alleles = np.sum(new_variants[mask] != original_haplotypes[mask])
        changed_aleles_percentage = changed_alleles / np.sum(mask)
        print(f"Mean percentage of changed alleles: {changed_aleles_percentage}")
        print(f"Mean difference between new and original haplotypes: {diff}")

        results = {
            "stacker_strategy": strategy,
            "changed_alleles_percentage": 100 * float(changed_aleles_percentage),
            "mean_difference": float(diff),
            "number_of_changed_alleles": int(changed_alleles),
        }
        print(results)
        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy", new_variants)
        print("Saved new haplotypes")
        with open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/stacked.json", "w") as f:
            json.dump(results, f, indent=4)
        print("Wrote results to json")
    elif strategy == "to_unedited":
        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy", original_haplotypes)
        results = {
            "stacker_strategy": strategy,
            "changed_alleles_percentage": 0,
            "mean_difference": 0,
            "number_of_changed_alleles": 0,
        }
        print(results)
        with open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/stacked.json", "w") as f:
            json.dump(results, f, indent=4)
        print("Wrote results to json")
    else:
        raise NotImplementedError(f"Strategy {strategy} not implemented")

if __name__ == "__main__":
    main()
