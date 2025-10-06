import numpy as np
import json
import pickle

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

#make profile timing decorator
def profile(func):
    def wrapper(*args, **kwargs):
        import time
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Function {func.__name__} took {end_time - start_time:.4f} seconds to run")
        return result
    return wrapper

class HaplotypeHMM:
    def __init__(self, chromosome, subject_id):
        self.chromosome = chromosome
        self.subject_id = subject_id

        self.thousand_g = np.load(f"{base_path}/chr{chromosome}/1000g.npy")
        self.thousand_g_positions = np.load(f"{base_path}/chr{chromosome}/1000g_positions.npy")
        self.thousand_g_subjects = np.load(f"{base_path}/chr{chromosome}/1000g_subjects.npy")

        self.thousand_g_haplotypes = np.concatenate(
            [self.thousand_g[:, :, 0], self.thousand_g[:, :, 1]], axis=0
        )

        self.snp_positions = pickle.load(open(f"{base_path}/chr{chromosome}/pangenome_to_thousand_g_alignments.pickle", "rb")) #all alignments are snps
        self.snp_positions_set = set(self.snp_positions.keys())

        self.blocks = json.load(open(f"{base_path}/chr{chromosome}/blocks_dict.json"))

        self.ne = self.thousand_g.shape[0]
        self.r = 1.26

        self.num_states = self.thousand_g_haplotypes.shape[0]
        print("Loaded HMM data for chromosome", chromosome)

    def get_anchor_snps(self, block_idx):
        block = self.blocks[str(block_idx)]
        block_anchor_snps_positions = [(snp, block[1][j]) for j,snp in enumerate(block[1]) if snp in self.snp_positions_set]
        anchor_snp_idx_within_block = [j for j,snp in enumerate(block[1]) if snp in self.snp_positions_set]
        return np.array(block_anchor_snps_positions), np.array(anchor_snp_idx_within_block)

    def get_transition_matrix(self, pos_1, pos_2):
        d = (1 / 1_000_000) * (pos_2 - pos_1) * self.ne * self.r
        exp_term = np.exp(-d / self.num_states)
        p = exp_term + (1 - exp_term) / self.num_states
        q = (1 - exp_term) / self.num_states

        log_p = np.log(p)
        log_q = np.log(q)

        # Create matrix filled with q, then fill diagonal with p
        trans = np.full((self.num_states, self.num_states), log_q)
        np.fill_diagonal(trans, log_p)
        
        #assert np.allclose(np.sum(trans, axis=1), 1)
        #assert np.allclose(np.sum(trans, axis = 0), 1)
        return trans

    @profile
    def get_transition_probabilities(self, anchor_snps):
        """
        Compute the forward algorithm for a given haplotype and block index.
        """
        # Initialize forward probabilities
        transitions_probs = []
        anchor_snps_positions = anchor_snps[:, 0]

        for i in range(len(anchor_snps_positions) - 1):
            snp1_pos = anchor_snps_positions[i]
            snp2_pos = anchor_snps_positions[i + 1]
            transition_prob = self.get_transition_matrix(snp1_pos, snp2_pos)
            transitions_probs.append(transition_prob)
        return transitions_probs
    
    @profile
    def forward_algorithm(self, haplotype, block_idx):
        """
        Perform the forward algorithm to compute the log-probability of observing a haplotype
        through the HMM for a given block.
        """
        anchor_snps,_ = self.get_anchor_snps(block_idx)

        # Extract anchor SNP indices and positions
        anchor_indices = anchor_snps[:, 1].astype(int)

        # Observed haplotype SNPs at anchor positions
        observed = haplotype[anchor_indices]  # shape: (num_snps,)

        # thousand_g haplotypes at anchor positions
        all_haplotypes = self.thousand_g_haplotypes[:, anchor_indices]  # shape: (num_states, num_snps)

        # Number of states and observations
        num_states = self.num_states
        num_snps = len(anchor_indices)

        # Compute emission log-probabilities: log(1-eps) if match, log(eps) if mismatch
        eps = 1e-4  # small error probability for sequencing/mutation
        emissions = np.where(all_haplotypes == observed, np.log(1 - eps), np.log(eps))  # shape: (num_states, num_snps)

        # Initialize alpha (log-space)
        alpha = np.zeros((num_states, num_snps))
        alpha[:, 0] = emissions[:, 0] - np.log(num_states)  # uniform initial probability

        # Compute transition matrices between anchor SNPs
        transition_matrices = self.get_transition_probabilities(anchor_snps)

        for t in range(1, num_snps):
            trans = transition_matrices[t - 1]  # log-transitions between t-1 and t

            # log-sum-exp trick for stability
            for j in range(num_states):
                log_probs = alpha[:, t - 1] + trans[:, j]
                max_log_prob = np.max(log_probs)
                alpha[j, t] = emissions[j, t] + max_log_prob + np.log(np.sum(np.exp(log_probs - max_log_prob)))

        # Final log-probability (sum over last column of alpha)
        final_log_probs = alpha[:, -1]
        max_final = np.max(final_log_probs)
        total_log_prob = max_final + np.log(np.sum(np.exp(final_log_probs - max_final)))

        return total_log_prob

    def sample_block_prior(self, block_idx):
        """
        Sample a haplotype block from the prior (unconditioned HMM path).
        Returns a haplotype array of just the SNPs in the block (anchor + non-anchor), in block order.
        """
        block = self.blocks[str(block_idx)]
        anchor_snps, anchor_snp_idx_with_block = self.get_anchor_snps(block_idx)
        anchor_positions = anchor_snps[:, 0].astype(int)
        anchor_indices = anchor_snps[:, 1].astype(int)

        num_snps = len(anchor_indices)
        num_states = self.num_states_without_subject

        # Flatten the thousand_g: 2N haplotypes x all SNPs
        all_haplotypes_full = self.thousand_g_haplotypes[:, anchor_indices]  # shape: (num_states, num_snps)
        all_haplotypes_anchors = all_haplotypes_full[:, anchor_indices]  # shape: (num_states, num_anchor_snps)

        # Sample state path through anchor SNPs
        sampled_states = np.zeros(num_snps, dtype=int)
        sampled_states[0] = np.random.choice(num_states)

        transition_matrices = self.get_transition_matrix(anchor_snps)
        for t in range(1, num_snps):
            trans = transition_matrices[t - 1]  # log-space
            probs = np.exp(trans[sampled_states[t - 1]])
            probs /= probs.sum()
            sampled_states[t] = np.random.choice(num_states, p=probs)

        # Create block-local haplotype
        block_positions = block[0]
        block_indices = block[1]
        haplotype_block = np.zeros(len(block_positions), dtype=int)

        # Fill anchor SNPs
        for i,idx in enumerate(anchor_snp_idx_with_block):
            haplotype_block[idx] = all_haplotypes_anchors[sampled_states[i], i]

        # Fill non-anchor SNPs
            # Fill non-anchor SNPs
        anchor_pos_to_state = dict(zip(anchor_positions, sampled_states))
        anchor_pos_array = np.array(anchor_positions)

        for i, (pos, global_idx) in enumerate(zip(block_positions, block_indices)):
            if pos in anchor_pos_to_state:
                continue  # already filled
            closest_anchor_idx = np.argmin(np.abs(anchor_pos_array - pos))
            closest_anchor_pos = anchor_pos_array[closest_anchor_idx]
            hap_idx = anchor_pos_to_state[closest_anchor_pos]
            haplotype_block[i] = all_haplotypes_full[hap_idx, global_idx]

        return haplotype_block

    def fill_block(self, haplotype, block_idx, block_haplotype):
        """
        Fill a block with the given haplotype.
        """
        block = self.blocks[str(block_idx)]
        block_idx = block[1]

        haplotype[block_idx] = block_haplotype