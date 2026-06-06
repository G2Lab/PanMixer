import sys
import json
import numpy as np
import pickle
import os
import time
import tracemalloc

from hmm import HaplotypeHMM

from paths import BASE_PATH

base_path = f"{BASE_PATH}/starting_data/"
INFINITY = 1e200

def get_mem_mb():
    current, peak = tracemalloc.get_traced_memory()
    return current / 1024 / 1024, peak / 1024 / 1024

def get_pmi(chromosome, subject_id):
    t_start = time.time()
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

        pass
    
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
    cur, peak = get_mem_mb()
    print(f"[TIMING] get_pmi completed in {time.time() - t_start:.2f}s | [MEMORY] current: {cur:.1f}MB, peak: {peak:.1f}MB")
    return pmi_blocks

def get_support_and_pmi(chromosome, subject_id):
    tracemalloc.start()
    t_total = time.time()

    # --- 1. Compute utility vector ---
    t_step = time.time()
    pangenome = np.load(base_path + f"chr{chromosome}/pangenome.npy")
    not_missing = np.sum(pangenome != -1, axis=(0,2))

    valid_positions = np.where(not_missing > 0)[0]

    utility_loss = np.zeros(len(not_missing))
    utility_loss[valid_positions] = 1 / not_missing[valid_positions]

    maximum_utility_loss = np.sum(utility_loss)

    cur, peak = get_mem_mb()
    print(f"[TIMING] Utility vector computed in {time.time() - t_step:.2f}s | [MEMORY] current: {cur:.1f}MB, peak: {peak:.1f}MB")

    # --- 2. Compute PMI for each block ---
    pmi_blocks = get_pmi(chromosome, subject_id)

    # INFINITY check logic assumed here
    #assert np.sum(pmi_blocks) < INFINITY...

    # Save PMI blocks
    os.makedirs(base_path + f"chr{chromosome}/subjects/{subject_id}/", exist_ok=True)
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/pmi_blocks.npy", pmi_blocks)
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/utility_loss.npy", utility_loss)

    # --- 3. Find a resampling for each block ---
    original_haplotypes = pangenome[subject_id]
    
    # Initialize new haplotypes to -1
    new_haplotypes = np.full_like(original_haplotypes, -1)

    hmm = HaplotypeHMM(chromosome, subject_id)
    blocks = json.load(open(base_path + f"chr{chromosome}/blocks_dict.json"))
    population_af = np.load(base_path + f"chr{chromosome}/allele_frequencies.npy")

    cur, peak = get_mem_mb()
    print(f"[TIMING] PMI + data loading completed in {time.time() - t_step:.2f}s | [MEMORY] current: {cur:.1f}MB, peak: {peak:.1f}MB")

    # Track changes for real-time logging
    t_step = time.time()
    total_sites_processed = 0
    total_sites_changed = 0

    same_blocks = 0
    different_blocks = 0
    total_blocks = 0

    for j in range(2): # For each haplotype strand
        for i in blocks.keys():
            block_indices = blocks[i][1]
            pos_indicies = blocks[i][0]
            anchor_snps_idx = hmm.get_anchor_snps(i)[0]

            # Logic for Small Blocks (Random Sampling)
            if len(block_indices) == 1 or len(anchor_snps_idx) <= 1:
                for block_index, pos_index in zip(block_indices, pos_indicies):
                    original_variant = original_haplotypes[block_index, j]
                    
                    # Fetch frequencies
                    allele_frequencies = np.nan_to_num(population_af[block_index])
                    new_frequencies = allele_frequencies.copy()
                    new_frequencies_normalized = new_frequencies / np.nansum(new_frequencies)

                    if np.isnan(new_frequencies_normalized).any():
                        new_haplotypes[block_index, j] = 0
                        continue

                    # Pick random new allele
                    new_allele = np.random.choice(
                        range(len(allele_frequencies)), 
                        p=new_frequencies_normalized
                    )
                    
                    # Update variant
                    new_haplotypes[block_index, j] = new_allele 

            # Logic for Large Blocks (HMM Sampling)
            else:
                #print(f"Block {i} (Hap {j}) is a large block")
                sample_new_block = hmm.sample_block_prior(i)
                
                # --- VERIFICATION 1: Compare against ORIGINAL, not NEW ---
                original_segment = original_haplotypes[block_indices, j]
                
                # Mask out missing data (-1) from the comparison to ensure fairness
                valid_mask = (original_segment != -1)

                
                if np.sum(valid_mask) > 0:
                    diff = np.mean(sample_new_block[valid_mask] != original_segment[valid_mask])
                    
                    # Only print if the difference is strangely low (e.g., 0%)
                    if diff == 0.0:
                        #print(f"WARNING: Block {i} (Hap {j}) resulted in identical resampling.")

                        anchor_snp_idx, anchor_snps_within = hmm.get_anchor_snps(i)
                        anchor_snp_idx = anchor_snp_idx[:, 0]
                        block_across_pop = pangenome[:, anchor_snp_idx]

                        # flatten the haplotypes across the two each strand
                        #block_across_pop.shape = (subjects, anchors, 2)
                        # we want to flatten the last axis
                        #block_across_pop.shape = (subjects * 2, anchors)
                        # then we can get unique haplotypes
                        flat_haplotypes = block_across_pop.transpose(0, 2, 1).reshape(-1, len(anchor_snp_idx))
                        # 4. Get Uniques and Counts
                        uniques, counts = np.unique(flat_haplotypes, axis=0, return_counts=True)

                        anchor_snps_selected = sample_new_block[anchor_snps_within]
                        
                        # 5. Sort by frequency (descending)
                        sorted_idx = np.argsort(-counts)
                        uniques = uniques[sorted_idx]
                        counts = counts[sorted_idx]
                        total_haplotypes = np.sum(counts)

                        assert np.array_equal(original_segment, pangenome[subject_id, block_indices, j]), "Original segment not found in pangenome"

                        original_segment_within = original_segment[anchor_snps_within]


                        index_selected = np.where((uniques == anchor_snps_selected).all(axis=1))[0]
                        if len(index_selected) == 0:
                            index_selected = -1
                        else:
                            index_selected = index_selected[0]   
                        #print(f"Total Haplotypes in Pool: {total_haplotypes}")
                        #print(f"Number of Unique Patterns: {len(uniques)}")
                        #print(f"Top 10 Haplotypes in this block:")
                        #print(f"{'Count':<10} {'Freq (%)':<10} {'Pattern (First 10 SNPs)'}")
                        #print("-" * 50)

                        #for l, (u, c) in enumerate(zip(uniques[:90], counts[:90])):
                        #    freq = (c / total_haplotypes) * 100
                        #    # Convert array to string for cleaner printing (truncate if long)
                        #    pattern_str = "".join(str(x) for x in u[:10]) 
                        #    if len(u) > 10: pattern_str += "..."
                        #    
                        #    # Highlight if this matches our "stuck" segment (need to match shapes carefully)
                        #    # (Optional: simply printing the table is usually enough to see the problem)
                        #    if index_selected == l:
                        #        #make this green
                        #        print(f"\033[92m{c:<10} {freq:<10.2f} {pattern_str}\033[0m")
                        #    else:
                        #        print(f"{c:<10} {freq:<10.2f} {pattern_str}")
                        #
                        #print("---------------------------------------\n")
                        
                    
                    total_blocks += 1
                    if diff == 0.0:
                        same_blocks += 1
                    else:
                        different_blocks += 1
                    

                    #print(f"Running stats: {same_blocks}/{total_blocks} blocks are the same.")
                        
                
                hmm.fill_block(new_haplotypes[:, j], i, sample_new_block)

    cur, peak = get_mem_mb()
    print(f"[TIMING] Resampling completed in {time.time() - t_step:.2f}s | [MEMORY] current: {cur:.1f}MB, peak: {peak:.1f}MB")

    # --- 4. GLOBAL VERIFICATION CHECKS ---
    valid_mask = (original_haplotypes != -1) & (new_haplotypes != -1)
    n_valid = np.sum(valid_mask)

    if n_valid == 0:
        raise ValueError("CRITICAL: No valid positions found to compare after processing.")

    changed_mask = (original_haplotypes != new_haplotypes) & valid_mask
    n_changed = np.sum(changed_mask)
    percent_changed = (n_changed / n_valid) * 100

    if percent_changed == 0:
        raise RuntimeError("OBSFUCATION FAILED: New haplotypes are identical to the original.")

    #save new haplotypes
    np.save(base_path + f"chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy", new_haplotypes)

    #4 compute support
    t_step = time.time()
    support_blocks = np.zeros((len(blocks), 2), dtype=float)
    for j, key in enumerate(blocks.keys()):
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

    cur, peak = get_mem_mb()
    print(f"[TIMING] Support blocks computed in {time.time() - t_step:.2f}s | [MEMORY] current: {cur:.1f}MB, peak: {peak:.1f}MB")
    print(f"[TIMING] Total get_support_and_pmi completed in {time.time() - t_total:.2f}s | [MEMORY] peak: {peak:.1f}MB")
    tracemalloc.stop()


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4, 5):
        print("Usage: python get_support_and_pmi.py <chromosome> <subject_id> [--seed SEED]")
        sys.exit(1)
    
    chromosome = sys.argv[1]
    subject_id = int(sys.argv[2])
    if len(sys.argv) == 5 and sys.argv[3] == "--seed":
        seed = int(sys.argv[4])
    elif len(sys.argv) == 4:
        seed = None if sys.argv[3] == "" else int(sys.argv[3])
    elif len(sys.argv) == 3:
        seed = None
    else:
        print("Usage: python get_support_and_pmi.py <chromosome> <subject_id> [--seed SEED]")
        sys.exit(1)
    if seed is not None:
        np.random.seed(seed)
        print(f"Using random seed {seed}")
    get_support_and_pmi(chromosome, subject_id)
