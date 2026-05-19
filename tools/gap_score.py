import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data, add_command, get_chromosome_path, load_data_multichromosome, profile
import pandas
import json
import os
import matplotlib.pyplot as plt

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
    launch_job_multichromosome("gap_score", args, memory="120g", cpus="1", num_tasks=str(num_tasks))

def score_genotypes_remove_missing(individual_haplotypes, other_haplotypes, alt_af):
    # Per-site validity: all 4 alleles (2 individual + 2 other) must be non-missing.
    indiv_valid = (individual_haplotypes != -1).all(axis=-1)   # (n_sites,)
    other_valid = (other_haplotypes != -1).all(axis=-1)        # (n_others, n_sites)

    individual_genotypes = np.sum(individual_haplotypes, axis = 1)
    other_genotypes = np.sum(other_haplotypes, axis = 2)

    p = alt_af           # (n_loci,)
    q = 1.0 - p

    # Remove singletons: positions where p=0 or p=1 (all homozygous)
    non_singleton_mask = (p > 0) & (p < 1)
    print("number of singletons removed", np.sum(~non_singleton_mask))

    # Apply mask to filter out singletons
    p = p[non_singleton_mask]
    q = q[non_singleton_mask]
    individual_genotypes = individual_genotypes[non_singleton_mask]
    other_genotypes = other_genotypes[:, non_singleton_mask]
    indiv_valid = indiv_valid[non_singleton_mask]
    other_valid = other_valid[:, non_singleton_mask]

    p00 = q*q            # P(geno=0)
    p01 = 2*p*q          # P(geno=1)
    p11 = p*p

    p00, p01, p11 = p00[None, :], p01[None, :], p11[None, :]

    other_probs = np.where(other_genotypes==0, p00, np.where(other_genotypes==1, p01, p11))

    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)

    both_valid = indiv_valid[None, :] & other_valid
    equal_indices = (individual_genotypes == other_genotypes) & both_valid
    print("number of equal genotypes", np.sum(equal_indices))

    score = log_probs * equal_indices

    # make a histogram of the log_probs with different colors for equal and opposite
    plt.hist(log_probs[equal_indices], bins=100, alpha=0.5, label="equal")
    plt.hist(log_probs[~equal_indices], bins=100, alpha=0.5, label="opposite")
    plt.legend()
    plt.savefig("test_hist.png")
    plt.close()

    #print("mean opposite score: ", (log_probs * (1-equal_indices)).mean())
    #print("mean score", score.mean())

    scores = np.sum(score, axis = 1)
    snps_equal = np.mean(equal_indices, axis = 1)

    return scores, snps_equal


def build_cache(chromosome):
    AF_file = STARTING_DATA_PATH + f"/chr{chromosome}/{AFS_FREQ_NPY}"

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"
    
    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    haplotypes = np.load(chromosome_path + "pangenome.npy")

    attack_db_masked = np.load(chromosome_path + ONEK_MASKED_NPY)

    frequencies = np.load(AF_file)

    pangenome_mask = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy")

    #are these only snps?
    zero_one_frequencies = np.sum(frequencies[pangenome_mask][:, :2], axis=1)
    snp_indices = zero_one_frequencies == 1

    attack_db_haplotypes = attack_db_masked[:, snp_indices]

    assert np.all(attack_db_haplotypes < 3)

    alt_af = frequencies[pangenome_mask][snp_indices, 1]

    attack_db_haplotypes = attack_db_haplotypes.astype(np.float32)
    other_genotypes = attack_db_haplotypes.sum(axis=2)

    p = alt_af
    q = 1.0 - p
    
    p00 = q*q            # P(geno=0)
    p01 = 2*p*q          # P(geno=1)
    p11 = p*p 

    p00, p01, p11 = p00[None, :], p01[None, :], p11[None, :]

    other_probs = np.where(other_genotypes==0, p00, np.where(other_genotypes==1, p01, p11))
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)
    
    np.save(f"{BASE_PATH}/starting_data/chr{chromosome}/genotype_scores_log_probs.npy", log_probs)
    return

def score_genotypes_remove_missing_cached(individual_haplotypes, log_probs, other_genotypes):
    individual_genotypes = np.sum(individual_haplotypes, axis = 1)
    other_genotypes = np.sum(other_genotypes, axis = 2)
    
    equal_indices = (individual_genotypes == other_genotypes) & (individual_genotypes != -1)
    score = log_probs * equal_indices

    scores = np.sum(score, axis = 1)
    snps_equal = np.mean(equal_indices, axis = 1)
    
    return scores, snps_equal

def main():
    experiment_number_or_build_cache = sys.argv[1]
    if experiment_number_or_build_cache == "build_cache":
        chromosome = int(sys.argv[2])
        build_cache(chromosome)
        return
    
    experiment_number = int(experiment_number_or_build_cache)
    
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]

    does_cache_exist = os.path.exists(f"{BASE_PATH}/starting_data/chr{chromosome}/genotype_scores_log_probs.npy")

    #if os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json"):
    #    return

    AF_file = STARTING_DATA_PATH + f"/chr{chromosome}/{AFS_FREQ_NPY}"

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"
    
    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    haplotypes = np.load(chromosome_path + "pangenome.npy")

    attack_db_masked = np.load(chromosome_path + ONEK_MASKED_NPY)

    frequencies = np.load(AF_file)

    subject_index_pangenome = np.where(subjects == subject_name)[0][0]

    interested_haplotypes = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")

    pangenome_mask = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy")

    #are these only snps?
    zero_one_frequencies = np.sum(frequencies[pangenome_mask][:, :2], axis=1)
    snp_indices = zero_one_frequencies == 1

    # Restrict to bi-allelic SNP positions (REF=1bp, single ALT=1bp).
    biallelic_snp_mask_path = f"{BASE_PATH}/starting_data/chr{chromosome}/biallelic_snp_mask.npy"
    if not os.path.exists(biallelic_snp_mask_path):
        raise FileNotFoundError(
            f"{biallelic_snp_mask_path} missing — run starting_data/scripts/src/build_biallelic_snp_mask.py first"
        )
    biallelic_snp_full = np.load(biallelic_snp_mask_path)
    biallelic_aligned = biallelic_snp_full[pangenome_mask]
    snp_indices = snp_indices & biallelic_aligned
    print(f"Using {int(snp_indices.sum())} bi-allelic SNP sites after filtering", flush=True)

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

    if not does_cache_exist:
        for i in range(NUM_CHUNKS):
            scores_genotypes_chunk, _ = score_genotypes_remove_missing(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
            scores_genotypes[chunks[i]] = scores_genotypes_chunk
    else:
        log_probs = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/genotype_scores_log_probs.npy")
        scores_genotypes, _ = score_genotypes_remove_missing_cached(interested_haplotypes, log_probs, attack_db_haplotypes)
    
    score_genotypes_against_self, percent_shared_genotypes_against_self = score_genotypes_remove_missing(interested_haplotypes, original_haplotypes, alt_af)
    highest_score_genotypes = np.max(scores_genotypes)
    genotypes_g_to_gstar = score_genotypes_against_self[0] - highest_score_genotypes

    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores.npy", scores_genotypes)
    np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores_self.npy", score_genotypes_against_self)

    results = {
        "in_data": True,
        "genotype_g_to_gstar": float(genotypes_g_to_gstar),
        "genotype_highest_score": float(highest_score_genotypes),
        "genotypes_percent_shared": float(percent_shared_genotypes_against_self[0]),

        "haplotype_g_to_gstar": None,
        "haplotype_highest_score": None,
        "haplotypes_percent_shared": None,

        "haplotype_both_g_to_gstar": None,
        "haplotype_both__highest_score": None,
        "haplotypes_both_percent_shared": None,

        "no_weight_haplotype_g_to_gstar": None,
    }

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json", "w") as f:
        json.dump(results, f, indent=4)
    
if __name__ == '__main__':
    main()