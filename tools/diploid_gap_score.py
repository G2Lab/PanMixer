import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data, add_command, get_chromosome_path, load_data_multichromosome, profile
import pandas
import json
import os


from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
    AFS_FREQ_NPY,
    ONEK_MASKED_NPY,
    BASE_PATH,
)

def diploid_gap_score_computer(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [experiment_number]
    launch_job_multichromosome("diploid_gap_score", args, memory="120g", cpus="1", num_tasks=str(num_tasks))

@profile
def score_haplotypes(individual_haplotypes, other_haplotypes, alt_af):
    ref_af = alt_af

    alt_af_broadcasted = alt_af[np.newaxis, :, np.newaxis]
    ref_af_broadcasted = ref_af[np.newaxis, :, np.newaxis]

    other_probs = other_haplotypes * alt_af_broadcasted + (1 - other_haplotypes) * ref_af_broadcasted
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)

    equal_indices = individual_haplotypes == other_haplotypes

    score = log_probs * equal_indices

    scores = np.sum(score, axis=(1,2))
    haplotypes_same = np.mean(equal_indices, axis=(1,2))

    return scores, haplotypes_same

@profile
def score_genotypes(individual_haplotypes, other_haplotypes, alt_af):
    individual_genotypes = np.sum(individual_haplotypes, axis = 1)
    other_genotypes = np.sum(other_haplotypes, axis = 2)

    p = alt_af           # (n_loci,)
    q = 1.0 - p
    p00 = q*q            # P(geno=0)
    p01 = 2*p*q          # P(geno=1)
    p11 = p*p 

    p00, p01, p11 = p00[None, :], p01[None, :], p11[None, :]

    other_probs = np.where(other_genotypes==0, p00, np.where(other_genotypes==1, p01, p11))
    
    print("number of singletons", np.sum(other_probs == 0))
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)
    
    equal_indices = individual_genotypes == other_genotypes

    score = log_probs * equal_indices

    scores = np.sum(score, axis = 1)
    snps_equal = np.mean(equal_indices, axis = 1)

    return scores, snps_equal

@profile
def score_no_weight_haplotype(individual_haplotypes, other_haplotypes, _):
    equal_indices = individual_haplotypes == other_haplotypes
    return np.sum(equal_indices, axis = (1,2)), None

@profile
def score_both_haplotypes_equal(individual_haplotypes, other_haplotypes, alt_af):
    ref_af = alt_af

    alt_af_broadcasted = alt_af[np.newaxis, :, np.newaxis]
    ref_af_broadcasted = ref_af[np.newaxis, :, np.newaxis]

    other_probs = other_haplotypes * alt_af_broadcasted + (1 - other_haplotypes) * ref_af_broadcasted
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)

    equal_indices = individual_haplotypes == other_haplotypes
    equal_indices_both = np.logical_and(equal_indices[:, :, 0], equal_indices[:, :, 1])

    log_probs_both = np.sum(log_probs, axis = 2)

    score = log_probs_both * equal_indices_both

    scores = np.sum(score, axis=1)
    haplotypes_both = np.mean(equal_indices_both, axis=1)

    return scores, haplotypes_both
    
def main():
    experiment_number = sys.argv[1]
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]

    if os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json"):
        return

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

    #are these only snps?
    zero_one_frequencies = np.sum(frequencies[pangenome_mask][:, :2], axis=1)
    snp_indices = zero_one_frequencies == 1

    attack_db_haplotypes = attack_db_masked[:, snp_indices]

    assert np.all(attack_db_haplotypes < 3)

    alt_af = frequencies[pangenome_mask][snp_indices, 1]

    interested_haplotypes = interested_haplotypes[pangenome_mask][snp_indices]

    print(np.sum(interested_haplotypes >= 2))
    assert np.all(interested_haplotypes < 2), f"Print found {np.sum(interested_haplotypes >= 2)} bad alleles"
    assert interested_haplotypes.shape[0] == attack_db_haplotypes.shape[1]
    assert interested_haplotypes.shape[0] == alt_af.shape[0]

    original_haplotypes = haplotypes[subject_index_pangenome]
    original_haplotypes = np.expand_dims(original_haplotypes, axis=0)
    original_haplotypes = original_haplotypes[:, pangenome_mask][:, snp_indices]


    ##### memory management
    attack_db_haplotypes = attack_db_haplotypes.astype(np.float32)
    interested_haplotypes = interested_haplotypes.astype(np.float32)
    original_haplotypes = original_haplotypes.astype(np.float32)

    import gc; gc.collect()

    #do chunking by 10 to reduce memory footprint:

    N_ATTACK = attack_db_haplotypes.shape[0]
    NUM_CHUNKS = 10
    chunks = np.array_split(np.arange(N_ATTACK), NUM_CHUNKS)
    
    ##########################
    # genotypes comparison:
    ##########################
    # if obfuscated is 0/1 and database is 0/1 = true
    # ignores phasing...


    scores_genotypes = np.zeros(N_ATTACK)
    for i in range(NUM_CHUNKS):
        scores_genotypes_chunk, _ = score_genotypes(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
        scores_genotypes[chunks[i]] = scores_genotypes_chunk
    
    score_genotypes_against_self, percent_shared_genotypes_against_self = score_genotypes(interested_haplotypes, original_haplotypes, alt_af)
    highest_score_genotypes = np.max(scores_genotypes)
    genotypes_g_to_gstar = score_genotypes_against_self[0] - highest_score_genotypes

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores.npy", scores_genotypes)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores_self.npy", score_genotypes_against_self)

    ########################
    # haplotype comparison
    ########################

    scores_haplotypes = np.zeros(N_ATTACK)
    for i in range(NUM_CHUNKS):
        scores_haplotypes_chunk, _ = score_haplotypes(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
        scores_haplotypes[chunks[i]] = scores_haplotypes_chunk

    score_haplotypes_against_self, percent_shared_haplotypes_against_self = score_haplotypes(interested_haplotypes, original_haplotypes, alt_af)
    highest_score_haplotypes = np.max(scores_haplotypes)
    haplotypes_g_to_gstar = score_haplotypes_against_self[0] - highest_score_haplotypes

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_scores.npy", scores_haplotypes)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_scores_self.npy", score_haplotypes_against_self)

    ########################
    # haplotype both equal
    ########################

    scores_haplotypes_both = np.zeros(N_ATTACK)
    for i in range(NUM_CHUNKS):
        scores_haplotypes_both_chunk, _ = score_both_haplotypes_equal(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
        scores_haplotypes_both[chunks[i]] = scores_haplotypes_both_chunk

    scores_haplotypes_both_against_self, percent_shared_haplotypes_both_against_self = score_both_haplotypes_equal(interested_haplotypes, original_haplotypes, alt_af)

    highest_score_haplotypes_both = np.max(scores_haplotypes_both)
    haplotypes_both_g_to_gstar = scores_haplotypes_both_against_self[0] - highest_score_haplotypes_both

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_both_scores.npy", scores_haplotypes_both)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_both_scores_self.npy", scores_haplotypes_both_against_self)

    #######################
    # No weight comparison
    #######################

    score_no_weight = np.zeros(N_ATTACK)
    for i in range(NUM_CHUNKS):
        score_no_weight_chunk, _ = score_no_weight_haplotype(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
        score_no_weight[chunks[i]] = score_no_weight_chunk

    scores_no_weight_against_self, _ = score_no_weight_haplotype(interested_haplotypes, original_haplotypes, alt_af)
    scores_no_weight_g_to_gstar = scores_no_weight_against_self[0] - np.max(score_no_weight)

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/no_weight_scores.npy", score_no_weight)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/no_weight_scores_self.npy", scores_no_weight_against_self)

    
    results = {
        "in_data": True,
        "genotype_g_to_gstar": float(genotypes_g_to_gstar),
        "genotype_highest_score": float(highest_score_genotypes),
        "genotypes_percent_shared": float(percent_shared_genotypes_against_self[0]),

        "haplotype_g_to_gstar": float(haplotypes_g_to_gstar),
        "haplotype_highest_score": float(highest_score_haplotypes),
        "haplotypes_percent_shared": float(percent_shared_haplotypes_against_self[0]),

        "haplotype_both_g_to_gstar": float(haplotypes_both_g_to_gstar),
        "haplotype_both__highest_score": float(highest_score_haplotypes_both),
        "haplotypes_both_percent_shared": float(percent_shared_haplotypes_both_against_self[0]),

        "no_weight_haplotype_g_to_gstar": float(scores_no_weight_g_to_gstar),
    }

    print(results)

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json", "w") as f:
        json.dump(results, f, indent=4)
    
if __name__ == '__main__':
    main()