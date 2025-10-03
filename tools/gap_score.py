import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data, add_command, get_chromosome_path, load_data_multichromosome, profile
import pandas
import json

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
    AFS_FREQ_NPY,
    ONEK_MASKED_NPY,
    BASE_PATH,
)

def gap_score_computer(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [experiment_number]
    launch_job_multichromosome("gap_score", args, memory="80g", cpus="1", num_tasks=str(num_tasks))

@profile
def score(individual_haplotype, other_haplotypes, population_weights):
    haplotypes_equal = individual_haplotype == other_haplotypes
    genotypes_equal = np.logical_and(haplotypes_equal[:, :, 0], haplotypes_equal[:, :, 1])
    main_missing_indices = individual_haplotype == -1 
    other_missing_indices = other_haplotypes == -1
    missing_indices = np.logical_or(main_missing_indices, other_missing_indices)
    missing_indices = np.logical_and(missing_indices[:, :, 0], missing_indices[:, :, 1])

    interested_indices = np.logical_and(genotypes_equal, np.logical_not(missing_indices))

    has_both_haplotypes = np.logical_and(individual_haplotype[:, 0] != -1, individual_haplotype[:, 1] != -1)
    non_missing_individual_haplotype = individual_haplotype[has_both_haplotypes]
    non_missing_individual_haplotype = non_missing_individual_haplotype.astype(int)
    # one hot encode the haplotypes
    max_num_alleles = population_weights.shape[1]
    one_hot_individual_haplotype = np.zeros((non_missing_individual_haplotype.shape[0], max_num_alleles, 2))

    one_hot_individual_haplotype[np.arange(non_missing_individual_haplotype.shape[0]), non_missing_individual_haplotype[:, 0], 0] = 1
    one_hot_individual_haplotype[np.arange(non_missing_individual_haplotype.shape[0]), non_missing_individual_haplotype[:, 1], 1] = 1

    score_individual = one_hot_individual_haplotype[:,:,0] * population_weights[has_both_haplotypes] + one_hot_individual_haplotype[:,:,1] * population_weights[has_both_haplotypes]
    score_individual = np.sum(score_individual, axis=1)

    score_individual_uncompressed = np.zeros(individual_haplotype.shape[0])
    score_individual_uncompressed[has_both_haplotypes] = score_individual

    percent_shared_variants = np.mean(interested_indices, axis=1)

    score_all = np.dot(interested_indices, score_individual_uncompressed)
    return score_all, percent_shared_variants

def main():
    experiment_number = sys.argv[1]
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]

    AF_file = STARTING_DATA_PATH + f"/chr{chromosome}/{AFS_FREQ_NPY}"

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"
    
    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    haplotypes = np.load(chromosome_path + "pangenome.npy")

    attack_db_masked = np.load(chromosome_path + ONEK_MASKED_NPY)

    frequencies = np.load(AF_file)

    subject_index_pangenome = np.where(subjects == subject_name)[0][0]

    interested_haplotypes = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")

    pangenome_mask = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy")

    log_frequencies = -1 * np.log(frequencies + 1e-10)
    log_frequencies = np.nan_to_num(log_frequencies, nan=0.0)

    scores, percent_shared_snps = score(interested_haplotypes[pangenome_mask], attack_db_masked, log_frequencies[pangenome_mask])

    original_haplotypes = haplotypes[subject_index_pangenome]
    original_haplotypes = np.expand_dims(original_haplotypes, axis=0)

    score_against_self, percent_shared_snps_against_self = score(interested_haplotypes[pangenome_mask], original_haplotypes[:, pangenome_mask], log_frequencies[pangenome_mask])

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score.npy", scores)

    score_against_self = np.array(score_against_self)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_against_self.npy", score_against_self)

    highest_score = np.max(scores)

    g_to_gstar = score_against_self[0] - highest_score

    results = {
        "g_to_gstar": g_to_gstar,
        "highest_score": highest_score,
        "in_data": True,
        "percent_shared_snps": percent_shared_snps_against_self[0],
    }

    print(results)

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score.json", "w") as f:
        json.dump(results, f, indent=4)
    
    
if __name__ == '__main__':
    main()