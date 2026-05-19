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

DO_GENOTYPE_SCORING = True
DO_HAPLOTYPE_SCORING = True
DO_BOTH_HAPLOTYPE_SCORING = True
DO_NO_WEIGHT_SCORING = True

def diploid_gap_score_computer(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [experiment_number]
    return launch_job_multichromosome("diploid_gap_score", args, memory="150g", cpus="1", num_tasks=str(num_tasks), max_concurrent=16)

@profile
def score_haplotypes(individual_haplotypes, other_haplotypes, alt_af):
    ref_af = alt_af

    alt_af_broadcasted = alt_af[np.newaxis, :, np.newaxis]
    ref_af_broadcasted = ref_af[np.newaxis, :, np.newaxis]

    other_probs = other_haplotypes * alt_af_broadcasted + (1 - other_haplotypes) * ref_af_broadcasted
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)

    equal_indices = (individual_haplotypes == other_haplotypes) & (individual_haplotypes != -1) 

    score = log_probs * equal_indices

    scores = np.sum(score, axis=(1,2))
    haplotypes_same = np.mean(equal_indices, axis=(1,2))

    return scores, haplotypes_same


def score_genotypes_remove_missing(individual_haplotypes, other_haplotypes, alt_af):
    # Per-site validity: all 4 alleles (2 from individual + 2 from each other) must be non-missing.
    # individual_haplotypes: (n_sites, 2); other_haplotypes: (n_others, n_sites, 2)
    indiv_valid = (individual_haplotypes != -1).all(axis=-1)   # (n_sites,)
    other_valid = (other_haplotypes != -1).all(axis=-1)        # (n_others, n_sites)

    individual_genotypes = np.sum(individual_haplotypes, axis = 1)
    other_genotypes = np.sum(other_haplotypes, axis = 2)

    p = alt_af           # (n_loci,)
    q = 1.0 - p

    # Remove singletons: positions where p=0 or p=1 (all homozygous)
    non_singleton_mask = (p > 0) & (p < 1)

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

@profile
def score_no_weight_haplotype(individual_haplotypes, other_haplotypes, _):
    equal_indices = (individual_haplotypes == other_haplotypes) & (individual_haplotypes != -1)
    return np.sum(equal_indices, axis = (1,2)), None

@profile
def score_both_haplotypes_equal(individual_haplotypes, other_haplotypes, alt_af):
    ref_af = alt_af

    alt_af_broadcasted = alt_af[np.newaxis, :, np.newaxis]
    ref_af_broadcasted = ref_af[np.newaxis, :, np.newaxis]

    other_probs = other_haplotypes * alt_af_broadcasted + (1 - other_haplotypes) * ref_af_broadcasted
    
    log_probs = -1 * np.log(other_probs + 1e-10)
    log_probs = np.nan_to_num(log_probs, nan=0.0)

    equal_indices = (individual_haplotypes == other_haplotypes) & (individual_haplotypes != -1)
    equal_indices_both = np.logical_and(equal_indices[:, :, 0], equal_indices[:, :, 1])

    log_probs_both = np.sum(log_probs, axis = 2)

    score = log_probs_both * equal_indices_both

    scores = np.sum(score, axis=1)
    haplotypes_both = np.mean(equal_indices_both, axis=1)

    return scores, haplotypes_both

def main():
    from datetime import datetime
    def hb(msg):
        print(f"[{datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)

    experiment_number = sys.argv[1]
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    hb(f"task started: exp={experiment_number} chr={chromosome} row_id={row_id}")

    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]
    hb(f"subject={subject_name}")

    #if os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json"):
    #    return

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"

    hb("loading pangenome_subjects.npy")
    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    hb(f"loading pangenome.npy (mmap)")
    haplotypes = np.load(chromosome_path + "pangenome.npy", mmap_mode='r')
    hb(f"loaded pangenome haplotypes shape={haplotypes.shape}")

    # Use 30x (3202 samples) attack database with relaxed (pos, ref) mask
    posref_mask_30x_path = f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask_posref_30x.npy"
    posref_onek_30x_path = chromosome_path + "1000g_30x_phased_masked_posref.npy"

    # Legacy paths (Phase 3, 2504 samples)
    posref_mask_path = f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask_posref.npy"
    posref_onek_path = chromosome_path + "1000g_phased_masked_posref.npy"

    if os.path.exists(posref_mask_30x_path) and os.path.exists(posref_onek_30x_path):
        hb("loading posref_mask_30x and 1000g_30x_phased_masked_posref (mmap)")
        pangenome_mask = np.load(posref_mask_30x_path)
        attack_db_masked = np.load(posref_onek_30x_path, mmap_mode='r')
        print(f"Using 30x attack DB (3202 samples), relaxed (pos,ref) mask: {len(pangenome_mask)} sites", flush=True)
    elif os.path.exists(posref_mask_path) and os.path.exists(posref_onek_path):
        hb("loading posref_mask and 1000g_phased_masked_posref (mmap)")
        pangenome_mask = np.load(posref_mask_path)
        attack_db_masked = np.load(posref_onek_path, mmap_mode='r')
        print(f"Falling back to Phase 3 attack DB (2504 samples), relaxed mask: {len(pangenome_mask)} sites", flush=True)
    else:
        hb("loading pangenome_mask and ONEK_MASKED_NPY (mmap)")
        pangenome_mask = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy")
        attack_db_masked = np.load(chromosome_path + ONEK_MASKED_NPY, mmap_mode='r')
        print(f"Falling back to strict mask: {len(pangenome_mask)} sites", flush=True)
    hb(f"attack_db loaded shape={attack_db_masked.shape}")

    subject_index_pangenome = np.where(subjects == subject_name)[0][0]

    interested_haplotypes = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")

    # Recompute alt AF directly from 1000G attack database (binarized, so values are 0/1)
    # alt_af = fraction of 1000G haplotypes carrying the alt allele
    n_haplotypes_1000g = attack_db_masked.shape[0] * 2
    alt_counts = np.sum(attack_db_masked, axis=(0, 2))
    alt_af = alt_counts / n_haplotypes_1000g

    # Remove monomorphic sites in 1000G
    polymorphic = (alt_af > 0) & (alt_af < 1)

    # Restrict to bi-allelic SNP positions only (mask is over full pangenome variant array;
    # project through pangenome_mask to align with attack_db_masked / alt_af).
    biallelic_snp_mask_path = chromosome_path + "biallelic_snp_mask.npy"
    if not os.path.exists(biallelic_snp_mask_path):
        raise FileNotFoundError(
            f"{biallelic_snp_mask_path} missing — run starting_data/scripts/src/build_biallelic_snp_mask.py first"
        )
    biallelic_snp_full = np.load(biallelic_snp_mask_path)
    biallelic_snp_aligned = biallelic_snp_full[pangenome_mask]
    site_mask = polymorphic & biallelic_snp_aligned
    print(
        f"Using {int(site_mask.sum())} bi-allelic SNP sites "
        f"({int(polymorphic.sum())} polymorphic, {int(biallelic_snp_aligned.sum())} bi-allelic SNPs after pangenome_mask, of {len(alt_af)} aligned)",
        flush=True,
    )

    attack_db_haplotypes = attack_db_masked[:, site_mask]
    alt_af = alt_af[site_mask]

    # Binarize pangenome haplotypes: 0=ref, >0→1, keep -1 as missing
    interested_haplotypes = interested_haplotypes[pangenome_mask][site_mask]
    interested_haplotypes = np.where(interested_haplotypes > 0, 1, interested_haplotypes)

    print(f"Attack database: {attack_db_haplotypes.shape[0]} subjects")
    print(f"Using {np.sum(polymorphic)} polymorphic sites for scoring (AF from 1000G)")
    assert interested_haplotypes.shape[0] == attack_db_haplotypes.shape[1]
    assert interested_haplotypes.shape[0] == alt_af.shape[0]

    # Use the subject's haplotypes from the 1000G 30x dataset if available,
    # rather than from the pangenome, for a fair gap score comparison
    onek_30x_subjects_path = chromosome_path + "1000g_30x_phased_subjects.npy"
    original_from_attack_db = False

    if os.path.exists(onek_30x_subjects_path):
        onek_30x_subjects = np.load(onek_30x_subjects_path)
        if subject_name in onek_30x_subjects:
            subject_idx_30x = np.where(onek_30x_subjects == subject_name)[0][0]
            onek_30x_masked_path = chromosome_path + "1000g_30x_phased_masked_posref.npy"
            if os.path.exists(onek_30x_masked_path):
                onek_30x_masked = np.load(onek_30x_masked_path, mmap_mode='r')
                original_haplotypes = np.array(onek_30x_masked[subject_idx_30x])[site_mask]
                original_haplotypes = np.expand_dims(original_haplotypes, axis=0)
                original_from_attack_db = True
                print(f"Using original haplotypes from 1000G 30x for {subject_name} (idx {subject_idx_30x})")

    if not original_from_attack_db:
        # Fallback to pangenome
        original_haplotypes = haplotypes[subject_index_pangenome]
        original_haplotypes = np.expand_dims(original_haplotypes, axis=0)
        original_haplotypes = original_haplotypes[:, pangenome_mask][:, site_mask]
        original_haplotypes = np.where(original_haplotypes > 0, 1, original_haplotypes)
        print(f"Using original haplotypes from pangenome for {subject_name} (not found in 1000G 30x)")

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
    if DO_GENOTYPE_SCORING:
        scores_genotypes = np.zeros(N_ATTACK)

        for i in range(NUM_CHUNKS):
            scores_genotypes_chunk, _ = score_genotypes_remove_missing(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
            scores_genotypes[chunks[i]] = scores_genotypes_chunk
        
        score_genotypes_against_self, percent_shared_genotypes_against_self = score_genotypes_remove_missing(interested_haplotypes, original_haplotypes, alt_af)
        highest_score_genotypes = np.max(scores_genotypes)
        genotypes_g_to_gstar = score_genotypes_against_self[0] - highest_score_genotypes

        print(f"Self score: {score_genotypes_against_self[0]:.2f}")
        print(f"Highest 1000G score: {highest_score_genotypes:.2f}")
        print(f"Gap (self - best): {genotypes_g_to_gstar:.2f}")

        # Debug: check target subject position in attack DB
        # Try 30x first, fall back to phase 3
        attack_subjects_path = chromosome_path + "1000g_30x_phased_subjects.npy"
        if not os.path.exists(attack_subjects_path):
            attack_subjects_path = chromosome_path + "1000g_phased_subjects.npy"
        if os.path.exists(attack_subjects_path):
            attack_subjects = np.load(attack_subjects_path)
            best_idx = np.argmax(scores_genotypes)
            if subject_name in attack_subjects:
                subject_idx_attack = np.where(attack_subjects == subject_name)[0][0]
                subject_score = scores_genotypes[subject_idx_attack]
                print(f"DEBUG: {subject_name} found in attack DB at index {subject_idx_attack}, score={subject_score:.2f}")
                print(f"DEBUG: Gap (self - subject in attack DB): {score_genotypes_against_self[0] - subject_score:.2f}")
                print(f"DEBUG: Rank of {subject_name}: {int(np.sum(scores_genotypes >= subject_score))}/{len(scores_genotypes)}")
            else:
                print(f"DEBUG: {subject_name} NOT in attack DB.")
            print(f"DEBUG: Best match: {attack_subjects[best_idx]} (score={highest_score_genotypes:.2f})")

        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores.npy", scores_genotypes)
        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/genotypes_scores_self.npy", score_genotypes_against_self)
    else:
        print("skipping genotype scoring")
    
    ########################
    # haplotype comparison
    ########################
    if DO_HAPLOTYPE_SCORING:
        scores_haplotypes = np.zeros(N_ATTACK)
        for i in range(NUM_CHUNKS):
            scores_haplotypes_chunk, _ = score_haplotypes(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
            scores_haplotypes[chunks[i]] = scores_haplotypes_chunk

        score_haplotypes_against_self, percent_shared_haplotypes_against_self = score_haplotypes(interested_haplotypes, original_haplotypes, alt_af)
        highest_score_haplotypes = np.max(scores_haplotypes)
        haplotypes_g_to_gstar = score_haplotypes_against_self[0] - highest_score_haplotypes

        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_scores.npy", scores_haplotypes)
        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_scores_self.npy", score_haplotypes_against_self)
    else:
        print("skipping haplotype scoring")

    ########################
    # haplotype both equal
    ########################

    if DO_BOTH_HAPLOTYPE_SCORING:
        scores_haplotypes_both = np.zeros(N_ATTACK)
        for i in range(NUM_CHUNKS):
            scores_haplotypes_both_chunk, _ = score_both_haplotypes_equal(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
            scores_haplotypes_both[chunks[i]] = scores_haplotypes_both_chunk

        scores_haplotypes_both_against_self, percent_shared_haplotypes_both_against_self = score_both_haplotypes_equal(interested_haplotypes, original_haplotypes, alt_af)

        highest_score_haplotypes_both = np.max(scores_haplotypes_both)
        haplotypes_both_g_to_gstar = scores_haplotypes_both_against_self[0] - highest_score_haplotypes_both

        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_both_scores.npy", scores_haplotypes_both)
        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/haplotypes_both_scores_self.npy", scores_haplotypes_both_against_self)
    else:
        print("skipping both haplotype scoring")

    #######################
    # No weight comparison
    #######################

    if DO_NO_WEIGHT_SCORING:
        score_no_weight = np.zeros(N_ATTACK)
        for i in range(NUM_CHUNKS):
            score_no_weight_chunk, _ = score_no_weight_haplotype(interested_haplotypes, attack_db_haplotypes[chunks[i]], alt_af)
            score_no_weight[chunks[i]] = score_no_weight_chunk

        scores_no_weight_against_self, _ = score_no_weight_haplotype(interested_haplotypes, original_haplotypes, alt_af)
        scores_no_weight_g_to_gstar = scores_no_weight_against_self[0] - np.max(score_no_weight)

        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/no_weight_scores.npy", score_no_weight)
        np.save(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/no_weight_scores_self.npy", scores_no_weight_against_self)
    else:
        print("skipping no weight scoring")

    results = {
        "in_data": True,
        "genotype_g_to_gstar": float(genotypes_g_to_gstar) if DO_GENOTYPE_SCORING else None,
        "genotype_highest_score": float(highest_score_genotypes) if DO_GENOTYPE_SCORING else None,
        "genotypes_percent_shared": float(percent_shared_genotypes_against_self[0]) if DO_GENOTYPE_SCORING else None,

        "haplotype_g_to_gstar": float(haplotypes_g_to_gstar) if DO_HAPLOTYPE_SCORING else None,
        "haplotype_highest_score": float(highest_score_haplotypes) if DO_HAPLOTYPE_SCORING else None,
        "haplotypes_percent_shared": float(percent_shared_haplotypes_against_self[0]) if DO_HAPLOTYPE_SCORING else None,

        "haplotype_both_g_to_gstar": float(haplotypes_both_g_to_gstar) if DO_BOTH_HAPLOTYPE_SCORING else None,
        "haplotype_both__highest_score": float(highest_score_haplotypes_both) if DO_BOTH_HAPLOTYPE_SCORING else None,
        "haplotypes_both_percent_shared": float(percent_shared_haplotypes_both_against_self[0]) if DO_BOTH_HAPLOTYPE_SCORING else None,

        "no_weight_haplotype_g_to_gstar": float(scores_no_weight_g_to_gstar) if DO_NO_WEIGHT_SCORING else None,
    }

    print(results)

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/gap_score_all.json", "w") as f:
        json.dump(results, f, indent=4)

    import resource
    print(f"Peak RSS: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB")

if __name__ == '__main__':
    main()