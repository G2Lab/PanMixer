import sys
import json
import numpy as np
import pickle
import os

from hmm import HaplotypeHMM

from constants import BASE_PATH
base_path = BASE_PATH + "/starting_data/"
INFINITY = 1e200

def get_best(chromosome, subject_id):

    pangenome = np.load(base_path + f"chr{chromosome}/pangenome.npy")

    possible_alleles = np.load(base_path + f"chr{chromosome}/num_alleles.npy")
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
                    assert new_haplotypes[block_index, j] < possible_alleles[block_index], f"New allele {new_haplotypes[block_index, j]} is not in the possible alleles {possible_alleles[block_index]} with frequencies {allele_frequencies}"
            else:
                sample_new_block = hmm.sample_block_prior(i)

                #print(f"Diff between original and new haplotypes: {np.mean(sample_new_block != new_haplotypes[block_indices, j])}")
                
                hmm.fill_block(new_haplotypes[:, j], i, sample_new_block)
                
                for block_index, pos_index in zip(block_indices, pos_indicies):
                    assert new_haplotypes[block_index, j] < possible_alleles[block_index], f"New allele {new_haplotypes[block_index, j]} is not in the possible alleles {possible_alleles[block_index]}"

    #save new haplotypes
    #np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy", new_haplotypes)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_mappings.py <chromosome> <subject_id>")
        sys.exit(1)
    
    chromosome = sys.argv[1]
    subject_id = int(sys.argv[2])
    get_best(chromosome, subject_id)