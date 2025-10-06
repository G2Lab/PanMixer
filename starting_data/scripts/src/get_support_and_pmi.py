import sys
import json
import numpy as np
import pickle
import os

from hmm_1000g import HaplotypeHMM

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

base_path = f"{BASE_PATH}/starting_data/"
INFINITY = 1e200

def get_pmi(chromosome, subject_id):
    blocks = json.load(open(base_path + f"chr{chromosome}/blocks_dict.json"))
    pangenome = np.load(base_path + f"chr{chromosome}/pangenome.npy")
    population_af = np.load(base_path + f"chr{chromosome}/allele_frequencies.npy")
    pangenome_positions = np.load(base_path + f"chr{chromosome}/pangenome_positions.npy")
        
    hmm = HaplotypeHMM(chromosome, subject_id)

    haplotypes = pangenome[subject_id]

    # one hot encode the haplotypes
    haplotype_one_hot_1 = np.zeros((haplotypes.shape[0], population_af.shape[1]), dtype=int)

    for j,haplotype in enumerate(haplotypes):
        if haplotype[0] != -1:
            haplotype_one_hot_1[j, haplotype[0]] = 1

    
    haplotype_one_hot_2 = np.zeros((haplotypes.shape[0], population_af.shape[1]), dtype=int)

    for j,haplotype in enumerate(haplotypes):
        if haplotype[1] != -1:
            haplotype_one_hot_2[j, haplotype[1]] = 1

    p = population_af

    #get locations where allele is not seen in the population

    predicate = np.sum(haplotype_one_hot_1 * p, axis=1) == 0

    pmi = -np.log(p)
    pmi = np.nan_to_num(pmi)

    pmi_1 = np.sum(haplotype_one_hot_1 * pmi, axis=1)
    pmi_1[haplotypes[:, 0] == -1] = 0

    pmi_2 = np.sum(haplotype_one_hot_2 * pmi, axis=1)
    pmi_2[haplotypes[:, 1] == -1] = 0

    pmi_subject = np.array([pmi_1, pmi_2]).T

    assert np.all(pmi_subject >= 0)
    assert np.all(np.isnan(pmi_subject) == False)

    way_too_high_idx = np.where(pmi_subject > INFINITY)[0]

    if len(way_too_high_idx) > 0:
        mappings_pangenome_to_thousand_g_phased = pickle.load(open(base_path + f"chr{chromosome}/pangenome_to_thousand_g_phased.pickle", "rb"))
        mappings_pangenome_to_thousand_g_alignments = pickle.load(open(base_path + f"chr{chromosome}/pangenome_to_thousand_g_alignments.pickle", "rb"))

        for i in way_too_high_idx:
            if i in mappings_pangenome_to_thousand_g_alignments:
                print(f"Position {i} is way too high in pmi_subject, but is in thousand g alignments")
            elif i in mappings_pangenome_to_thousand_g_phased:
                print(f"Position {i} is way too high in pmi_subject, but is in thousand g phased")
            else:
                print(f"Way too high pmi_subject: {len(way_too_high_idx)}")
    
    #FIX ME
    pmi_subject[way_too_high_idx] = 0

    keys = blocks.keys()

    pmi_blocks = np.zeros((len(keys), 2), dtype=float)

    for j, key in enumerate(keys):
        #print(f"Completed {j} out of {len(keys)} blocks", end="\r")
        if len(blocks[key][1]) == 0:
            continue
        if len(blocks[key][1]) == 1:
            for pos in blocks[key][1]:
                pmi_blocks[j] += pmi_subject[pos]
        elif len(hmm.get_anchor_snps(key)[0]) == 0:
            for pos in blocks[key][1]:
                pmi_blocks[j] += pmi_subject[pos]
        else:
            log_prob_hap_1 = hmm.forward_algorithm(haplotypes[:, 0], key)
            log_prob_hap_2 = hmm.forward_algorithm(haplotypes[:, 1], key)

            pmi_blocks[j][0] = -1 * log_prob_hap_1
            pmi_blocks[j][1] = -1 * log_prob_hap_2
    return pmi_blocks

def get_support_and_pmi(chromosome, subject_id):
    #1 compute utility vector

    pangenome = np.load(base_path + f"chr{chromosome}/pangenome.npy")
    not_missing = np.sum(pangenome != -1, axis=(0,2))

    valid_positions = np.where(not_missing > 0)[0]

    utility_loss = np.zeros(len(not_missing))
    utility_loss[valid_positions] = 1 / not_missing[valid_positions]

    maximum_utility_loss = np.sum(utility_loss)
    print(f"Maximum utility loss: {maximum_utility_loss}")

    #2 compute pmi for each block
    pmi_blocks = get_pmi(chromosome, subject_id)

    assert np.sum(pmi_blocks) < INFINITY, "PMI blocks are way too high likely due to missing data"

    # Save PMI blocks
    os.makedirs(base_path + f"chr{chromosome}/subjects/", exist_ok=True)
    os.makedirs(base_path + f"chr{chromosome}/subjects/{subject_id}/", exist_ok=True)

    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/pmi_blocks.npy", pmi_blocks)
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/utility_loss.npy", utility_loss)

    print(f"PMI blocks saved for chromosome {chromosome} and subject {subject_id}")

    #3 find a resampling for each block
    original_haplotypes = pangenome[subject_id]
    new_haplotypes = np.zeros_like(pangenome[subject_id])
    new_haplotypes[:] = -1

    hmm = HaplotypeHMM(chromosome, subject_id)
    blocks = json.load(open(base_path + f"chr{chromosome}/blocks_dict.json"))
    population_af = np.load(base_path + f"chr{chromosome}/allele_frequencies.npy")

    for j in range(2):
        for i in blocks.keys():
            block_indices = blocks[i][1]
            pos_indicies = blocks[i][0]
            anchor_snps_idx = hmm.get_anchor_snps(i)[0]

            if len(block_indices) == 1 or len(anchor_snps_idx) <= 1:
                #print(f"Process block {i} with {len(block_indices)} indices and {len(anchor_snps_idx)} anchor snps")
                for block_index, pos_index in zip(block_indices, pos_indicies):
                    original_variant = original_haplotypes[block_index, j]
                    allele_frequencies = population_af[block_index]
                    allele_frequencies = np.nan_to_num(allele_frequencies)
                    new_frequencies = allele_frequencies.copy()
                    new_frequencies_without_allele = [new_frequencies[allele] for allele in range(len(allele_frequencies)) if allele != original_variant]
                    
                    if len(new_frequencies_without_allele) == 0 or np.sum(new_frequencies_without_allele) == 0 or np.any(np.isnan(new_frequencies_without_allele)):
                        new_haplotypes[block_index, j] = original_variant
                        continue
                    new_frequencies_normalized = [f / sum(new_frequencies_without_allele) for f in new_frequencies_without_allele]

                    # pick random new allele based on normalized frequencies
                    new_allele = np.random.choice([allele for allele in range(len(allele_frequencies)) if allele != original_variant], 
                                                p=new_frequencies_normalized)
                    # Update the variant at the position
                    new_haplotypes[block_index, j] = new_allele
            else:
                sample_new_block = hmm.sample_block_prior(i)

                print(f"Diff between original and new haplotypes: {np.mean(sample_new_block != new_haplotypes[block_indices, j])}")
                
                hmm.fill_block(new_haplotypes[:, j], i, sample_new_block)

    #save new haplotypes
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy", new_haplotypes)

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
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/support_blocks.npy", support_blocks)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_mappings.py <chromosome> <subject_id>")
        sys.exit(1)
    
    chromosome = sys.argv[1]
    subject_id = int(sys.argv[2])
    get_support_and_pmi(chromosome, subject_id)