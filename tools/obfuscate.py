"""
Standalone obfuscation tool for a single subject on a single chromosome.
Runs the full pipeline: PMI/utility computation -> optimize -> stack -> VCF output.

Usage:
    python tools/obfuscate.py <data_dir> <subject_name> <capacity> <chromosome> <output_dir> [--seed SEED]

    data_dir:      path to starting_data directory (must contain chr{N}/ subdirs)
    subject_name:  sample name as it appears in the VCF (e.g. HG00438)
    capacity:      utility loss budget as fraction 0-1 (e.g. 0.1 = allow 10% utility loss)
    chromosome:    chromosome number (1-22)
    output_dir:    output directory (will create chr{N}/ subdir)
    seed:          optional random seed for HMM/allele-frequency sampling

Each chromosome is fully independent so this can be run as a SLURM array job.

Required files per chromosome in data_dir/chr{N}/:
    pangenome.npy, pangenome_positions.npy, pangenome_subjects.npy
    blocks_dict.json
    allele_frequencies.npy
    pangenome_to_thousand_g_alignments.pickle  (anchor SNP mappings for HMM)
    pangenome.vcf.gz  (template VCF for output)
"""

import numpy as np
import json
import pickle
import sys
import os
import time
import gzip
import resource
from scipy.special import logsumexp
from ortools.linear_solver import pywraplp


# ---------------------------------------------------------------------------
# HMM (from starting_data/scripts/src/hmm.py)
# ---------------------------------------------------------------------------

EFFECTIVE_N = 1.0 / 10_000.0
INFINITY = 1e200


def logsubexp(a, b):
    if np.any(b > a):
        raise ValueError("logsubexp requires a > b elementwise")
    return a + np.log1p(-np.exp(b - a))


class HaplotypeHMM:
    def __init__(self, chr_dir, pangenome, pangenome_positions, subject_id):
        self.pangenome = pangenome
        self.pangenome_positions = pangenome_positions
        self.subject_id = subject_id

        self.pangenome_haplotypes = np.concatenate(
            [pangenome[:, :, 0], pangenome[:, :, 1]], axis=0
        )

        self.snp_positions = pickle.load(
            open(f"{chr_dir}/pangenome_to_thousand_g_alignments.pickle", "rb")
        )
        self.snp_positions_set = set(self.snp_positions.keys())

        self.blocks = json.load(open(f"{chr_dir}/blocks_dict.json"))

        self.ne = pangenome.shape[0]
        self.r = 1.26
        self.num_states = self.pangenome_haplotypes.shape[0]

        self.pangenome_without_subject = np.delete(pangenome, subject_id, axis=0)
        self.pangenome_haplotypes_without_subject = (
            self.pangenome_without_subject
            .transpose(0, 2, 1)
            .reshape(-1, self.pangenome_without_subject.shape[1])
        )
        self.r_without_subject = 1.26
        self.num_states_without_subject = self.pangenome_haplotypes_without_subject.shape[0]

    def get_anchor_snps(self, block_idx):
        block = self.blocks[str(block_idx)]
        block_anchor_snps_positions = [
            (snp, block[1][j])
            for j, snp in enumerate(block[1])
            if snp in self.snp_positions_set
        ]
        anchor_snp_idx_within_block = [
            j for j, snp in enumerate(block[1]) if snp in self.snp_positions_set
        ]
        return np.array(block_anchor_snps_positions), np.array(anchor_snp_idx_within_block)

    def get_transitions(self, pos_1, pos_2):
        d = (1.0 / 1_000_000.0) * (pos_2 - pos_1) * EFFECTIVE_N * self.r_without_subject
        exp_term = np.exp(-d / self.num_states_without_subject)
        p = exp_term + (1 - exp_term) / self.num_states_without_subject
        q = (1 - exp_term) / self.num_states_without_subject
        return (np.log(p), np.log(q))

    def get_transition_matrix_without_subject(self, pos_1, pos_2):
        d = (1 / 1_000_000) * (pos_2 - pos_1) * 4 * EFFECTIVE_N * self.r
        exp_term = np.exp(-d / self.num_states_without_subject)
        p = exp_term + (1 - exp_term) / self.num_states_without_subject
        q = (1 - exp_term) / self.num_states_without_subject
        trans = np.full(
            (self.num_states_without_subject, self.num_states_without_subject), q
        )
        np.fill_diagonal(trans, p)
        return trans

    def get_transition_probabilities(self, anchor_snps):
        transitions = []
        positions = anchor_snps[:, 0]
        for i in range(len(positions) - 1):
            transitions.append(self.get_transitions(positions[i], positions[i + 1]))
        return transitions

    def get_transition_probabilities_without_subject(self, anchor_snps):
        transitions = []
        positions = anchor_snps[:, 0]
        for i in range(len(positions) - 1):
            transitions.append(
                self.get_transition_matrix_without_subject(positions[i], positions[i + 1])
            )
        return transitions

    def forward_algorithm(self, haplotype, block_idx):
        anchor_snps, _ = self.get_anchor_snps(block_idx)
        anchor_indices = anchor_snps[:, 1].astype(int)
        observed = haplotype[anchor_indices]
        all_haplotypes = self.pangenome_haplotypes[:, anchor_indices]

        num_states = self.num_states
        num_snps = len(anchor_indices)

        eps = 1e-4
        emissions = np.where(
            all_haplotypes == observed, np.log(1 - eps), np.log(eps)
        )

        log_alpha = np.log(np.ones(num_states) / num_states) + emissions[:, 0]

        transition_probs = self.get_transition_probabilities(anchor_snps)
        for t in range(1, num_snps):
            p, q = transition_probs[t - 1]
            p_alpha = log_alpha + p
            q_alpha = log_alpha + q
            log_alpha_sum_q = logsumexp(log_alpha) + q
            log_alpha = np.logaddexp(p_alpha, log_alpha_sum_q)
            log_alpha = logsubexp(log_alpha, q_alpha)
            log_alpha += emissions[:, t]

        return logsumexp(log_alpha)

    def sample_block_prior(self, block_idx):
        block = self.blocks[str(block_idx)]
        anchor_snps, anchor_snp_idx_within = self.get_anchor_snps(block_idx)
        anchor_positions = anchor_snps[:, 0].astype(int)
        anchor_indices = anchor_snps[:, 1].astype(int)

        num_snps = len(anchor_indices)
        num_states = self.num_states_without_subject

        all_hap_full = self.pangenome_haplotypes_without_subject
        all_hap_anchors = all_hap_full[:, anchor_indices]

        sampled_states = np.zeros(num_snps, dtype=int)
        sampled_states[0] = np.random.choice(num_states)

        transition_matrices = self.get_transition_probabilities_without_subject(anchor_snps)
        for t in range(1, num_snps):
            probs = transition_matrices[t - 1][sampled_states[t - 1]]
            probs = probs / probs.sum()
            sampled_states[t] = np.random.choice(num_states, p=probs)

        block_positions = block[0]
        block_indices = block[1]
        haplotype_block = np.zeros(len(block_positions), dtype=int)

        for i, idx in enumerate(anchor_snp_idx_within):
            haplotype_block[idx] = all_hap_anchors[sampled_states[i], i]

        anchor_pos_to_state = dict(zip(anchor_positions, sampled_states))
        anchor_pos_array = np.array(anchor_positions)

        for i, (pos, global_idx) in enumerate(zip(block_positions, block_indices)):
            if pos in anchor_pos_to_state:
                continue
            closest = np.argmin(np.abs(anchor_pos_array - pos))
            hap_idx = anchor_pos_to_state[anchor_pos_array[closest]]
            haplotype_block[i] = all_hap_full[hap_idx, global_idx]

        return haplotype_block

    def fill_block(self, haplotype, block_idx, block_haplotype):
        block = self.blocks[str(block_idx)]
        haplotype[block[1]] = block_haplotype


# ---------------------------------------------------------------------------
# Step 1: Compute PMI and utility (from get_support_and_pmi.py)
# ---------------------------------------------------------------------------

def compute_utility_and_pmi(chr_dir, pangenome, population_af, subject_id):
    """Compute utility_loss, pmi_blocks, support_blocks, and resampled haplotypes."""
    t0 = time.time()

    # Utility vector: 1 / (number of non-missing alleles at each position)
    not_missing = np.sum(pangenome != -1, axis=(0, 2))
    valid = np.where(not_missing > 0)[0]
    utility_loss = np.zeros(len(not_missing))
    utility_loss[valid] = 1.0 / not_missing[valid]
    print(f"  [1/4] Utility vector: {time.time()-t0:.1f}s", flush=True)

    # PMI per block
    t1 = time.time()
    pangenome_positions = np.load(f"{chr_dir}/pangenome_positions.npy")
    hmm = HaplotypeHMM(chr_dir, pangenome, pangenome_positions, subject_id)
    blocks = hmm.blocks
    keys = list(blocks.keys())
    haplotypes = pangenome[subject_id]

    # One-hot encode haplotypes for PMI
    num_alleles = population_af.shape[1]
    hap_oh_1 = np.zeros((haplotypes.shape[0], num_alleles), dtype=int)
    hap_oh_2 = np.zeros((haplotypes.shape[0], num_alleles), dtype=int)
    for j, hap in enumerate(haplotypes):
        if hap[0] != -1:
            hap_oh_1[j, hap[0]] = 1
        if hap[1] != -1:
            hap_oh_2[j, hap[1]] = 1

    pmi = -np.log(population_af)
    pmi = np.nan_to_num(pmi)

    pmi_1 = np.sum(hap_oh_1 * pmi, axis=1)
    pmi_1[haplotypes[:, 0] == -1] = 0
    pmi_2 = np.sum(hap_oh_2 * pmi, axis=1)
    pmi_2[haplotypes[:, 1] == -1] = 0
    pmi_subject = np.array([pmi_1, pmi_2]).T
    pmi_subject[pmi_subject > INFINITY] = 0

    pmi_blocks = np.zeros((len(keys), 2), dtype=float)
    for j, key in enumerate(keys):
        if len(blocks[key][1]) == 0:
            continue
        if len(blocks[key][1]) == 1 or len(hmm.get_anchor_snps(key)[0]) == 0:
            for pos in blocks[key][1]:
                pmi_blocks[j] += pmi_subject[pos]
        else:
            pmi_blocks[j][0] = -1 * hmm.forward_algorithm(haplotypes[:, 0], key)
            pmi_blocks[j][1] = -1 * hmm.forward_algorithm(haplotypes[:, 1], key)
    print(f"  [2/4] PMI blocks: {time.time()-t1:.1f}s", flush=True)

    # Resample haplotypes
    t2 = time.time()
    original_haplotypes = pangenome[subject_id]
    new_haplotypes = np.full_like(original_haplotypes, -1)

    for j in range(2):
        for i in keys:
            block_indices = blocks[i][1]
            anchor_snps_idx = hmm.get_anchor_snps(i)[0]

            if len(block_indices) == 1 or len(anchor_snps_idx) <= 1:
                for block_index in block_indices:
                    af = np.nan_to_num(population_af[block_index])
                    af_norm = af / np.nansum(af)
                    if np.isnan(af_norm).any():
                        new_haplotypes[block_index, j] = 0
                    else:
                        new_haplotypes[block_index, j] = np.random.choice(len(af), p=af_norm)
            else:
                sample = hmm.sample_block_prior(i)
                hmm.fill_block(new_haplotypes[:, j], i, sample)
    print(f"  [3/4] Resampling: {time.time()-t2:.1f}s", flush=True)

    # Compute support blocks
    t3 = time.time()
    support_blocks = np.zeros((len(keys), 2), dtype=float)
    for j, key in enumerate(keys):
        if len(blocks[key][1]) == 0:
            continue
        if len(blocks[key][1]) == 1:
            support_blocks[j] = utility_loss[blocks[key][1]]
        else:
            new_block = new_haplotypes[blocks[key][1]]
            orig_block = original_haplotypes[blocks[key][1]]
            ul_block = utility_loss[blocks[key][1]]
            ul_stacked = np.column_stack((ul_block, ul_block))
            support_blocks[j] = np.sum(ul_stacked * (new_block != orig_block), axis=0)
    print(f"  [4/4] Support blocks: {time.time()-t3:.1f}s", flush=True)

    return utility_loss, pmi_blocks, support_blocks, new_haplotypes


# ---------------------------------------------------------------------------
# Step 2: Optimizer (from tools/optimizer.py)
# ---------------------------------------------------------------------------

def optimize(weights, values, capacity):
    """LP solver: maximize sum(x * values) subject to sum(x * weights) <= capacity."""
    weights_flat = weights.flatten()
    values_flat = values.flatten()
    n = len(weights_flat)

    if capacity <= 0:
        return np.zeros(n, dtype=int).reshape(weights.shape)

    solver = pywraplp.Solver.CreateSolver("GLOP")
    x = [solver.NumVar(0, 1, f"x[{i}]") for i in range(n)]

    objective = solver.Objective()
    for i in range(n):
        objective.SetCoefficient(x[i], values_flat[i])
    objective.SetMaximization()

    constraint = solver.Constraint(-solver.infinity(), capacity)
    for i in range(n):
        constraint.SetCoefficient(x[i], weights_flat[i])

    solver.Solve()
    x_sol = np.array([x[i].solution_value() for i in range(n)])
    return x_sol.round().astype(int).reshape(weights.shape)


# ---------------------------------------------------------------------------
# Step 3: Stacker (from tools/stacker.py)
# ---------------------------------------------------------------------------

def stack(original_haplotypes, new_haplotypes, xsol, blocks):
    """Apply optimizer decisions: replace blocks where xsol==1 with resampled haplotypes."""
    result = original_haplotypes.copy()
    keys = list(blocks.keys())

    true_indices = np.argwhere(xsol == 1)
    for i, j in true_indices:
        for bi in blocks[keys[i]][1]:
            result[bi][j] = new_haplotypes[bi][j]

    return result


# ---------------------------------------------------------------------------
# Step 4: VCF output (from tools/convert_2_vcf.py)
# ---------------------------------------------------------------------------

BATCHED_WRITES = 1000


def write_vcf(genotypes, master_vcf, subject_name, output_path):
    """Replace subject's genotypes in the master VCF and write to output."""
    genotypes = np.round(genotypes).astype(int)

    open_func = gzip.open if master_vcf.endswith(".gz") else open
    subject_index = -1

    with open_func(master_vcf, "rt") as vcf_in, open(output_path, "w") as vcf_out:
        line_num = 0
        batch = []
        for line in vcf_in:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    subject_index = line.strip().split("\t").index(subject_name)
                vcf_out.write(line)
                continue

            if line_num % BATCHED_WRITES == 0:
                vcf_out.writelines(batch)
                batch = []

            fields = line.strip().split("\t")

            g0 = "." if genotypes[line_num, 0] == -1 else str(genotypes[line_num, 0])
            g1 = "." if genotypes[line_num, 1] == -1 else str(genotypes[line_num, 1])
            fields[subject_index] = f"{g0}|{g1}"

            batch.append("\t".join(fields) + "\n")
            line_num += 1

        if batch:
            vcf_out.writelines(batch)

    assert subject_index != -1, f"Subject {subject_name} not found in VCF header"
    assert line_num == genotypes.shape[0], (
        f"VCF lines ({line_num}) != genotypes rows ({genotypes.shape[0]})"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def get_rss_mb():
    """Current peak RSS in MB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def main():
    if len(sys.argv) not in (6, 7, 8):
        print(__doc__)
        sys.exit(1)

    data_dir = sys.argv[1]
    subject_name = sys.argv[2]
    capacity = float(sys.argv[3])
    chromosome = int(sys.argv[4])
    output_dir = sys.argv[5]
    if len(sys.argv) == 8 and sys.argv[6] == "--seed":
        seed = int(sys.argv[7])
    elif len(sys.argv) == 7:
        seed = None if sys.argv[6] == "" else int(sys.argv[6])
    elif len(sys.argv) == 6:
        seed = None
    else:
        print(__doc__)
        sys.exit(1)
    if seed is not None:
        np.random.seed(seed)

    chr_dir = f"{data_dir}/chr{chromosome}"
    out_chr_dir = f"{output_dir}/chr{chromosome}"
    os.makedirs(out_chr_dir, exist_ok=True)

    seed_msg = "unseeded" if seed is None else str(seed)
    print(f"=== Obfuscating {subject_name} | chr{chromosome} | capacity={capacity} | seed={seed_msg} ===", flush=True)
    t_total = time.time()
    timings = {}

    # --- Load shared data ---
    t_load = time.time()
    pangenome = np.load(f"{chr_dir}/pangenome.npy")
    subjects = np.load(f"{chr_dir}/pangenome_subjects.npy")
    population_af = np.load(f"{chr_dir}/allele_frequencies.npy")

    subject_id = int(np.where(subjects == subject_name)[0][0])
    timings["load_data_s"] = round(time.time() - t_load, 2)
    timings["load_data_rss_mb"] = round(get_rss_mb(), 1)
    print(f"Loaded data: {timings['load_data_s']}s | RSS={timings['load_data_rss_mb']:.0f}MB", flush=True)

    # --- Step 1: Compute PMI, utility, resample ---
    print("Step 1: Computing utility, PMI, and resampling...", flush=True)
    t_pmi = time.time()
    utility_loss, pmi_blocks, support_blocks, resampled = compute_utility_and_pmi(
        chr_dir, pangenome, population_af, subject_id
    )
    timings["graph_prep_s"] = round(time.time() - t_pmi, 2)
    timings["graph_prep_rss_mb"] = round(get_rss_mb(), 1)
    print(f"  Graph prep total: {timings['graph_prep_s']}s | RSS={timings['graph_prep_rss_mb']:.0f}MB", flush=True)

    # --- Step 2: Optimize ---
    print("Step 2: Optimizing...", flush=True)
    t_opt = time.time()
    max_utility_loss = 2 * np.sum(utility_loss)
    target = capacity * max_utility_loss
    xsol = optimize(support_blocks, pmi_blocks, target)

    pmi_gain = np.sum(xsol * pmi_blocks)
    util_loss = np.sum(xsol * support_blocks)
    timings["optimizer_s"] = round(time.time() - t_opt, 2)
    timings["optimizer_rss_mb"] = round(get_rss_mb(), 1)
    print(f"  Optimizer: {timings['optimizer_s']}s | RSS={timings['optimizer_rss_mb']:.0f}MB | "
          f"moves={np.sum(xsol)} | "
          f"utility_loss={util_loss:.2f}/{max_utility_loss:.2f} "
          f"({util_loss/max_utility_loss*100:.1f}%) | "
          f"pmi_gain={pmi_gain:.2f}", flush=True)

    # --- Step 3: Stack ---
    print("Step 3: Stacking...", flush=True)
    t_stack = time.time()
    blocks = json.load(open(f"{chr_dir}/blocks_dict.json"))
    original = pangenome[subject_id]
    final_haplotypes = stack(original, resampled, xsol, blocks)

    valid_mask = (original != -1) & (final_haplotypes != -1)
    changed = np.sum(original[valid_mask] != final_haplotypes[valid_mask])
    total_valid = np.sum(valid_mask)
    timings["stacker_s"] = round(time.time() - t_stack, 2)
    timings["stacker_rss_mb"] = round(get_rss_mb(), 1)
    print(f"  Stacker: {timings['stacker_s']}s | RSS={timings['stacker_rss_mb']:.0f}MB | "
          f"changed={changed}/{total_valid} alleles "
          f"({changed/total_valid*100:.2f}%)", flush=True)

    # --- Step 4: Write VCF ---
    print("Step 4: Writing VCF...", flush=True)
    t_vcf = time.time()
    vcf_path = f"{out_chr_dir}/obfuscated.vcf"
    master_vcf = f"{chr_dir}/pangenome.vcf.gz"
    write_vcf(final_haplotypes, master_vcf, subject_name, vcf_path)
    import shutil
    bgzip = shutil.which("bgzip") or os.path.join(os.path.dirname(sys.executable), "bgzip")
    bcftools = shutil.which("bcftools") or os.path.join(os.path.dirname(sys.executable), "bcftools")
    os.system(f"{bgzip} -f {vcf_path}")
    os.system(f"{bcftools} index {vcf_path}.gz")
    timings["vcf_write_s"] = round(time.time() - t_vcf, 2)
    timings["vcf_write_rss_mb"] = round(get_rss_mb(), 1)
    print(f"  VCF: {timings['vcf_write_s']}s | RSS={timings['vcf_write_rss_mb']:.0f}MB", flush=True)

    # --- Totals ---
    timings["total_s"] = round(time.time() - t_total, 2)
    timings["peak_rss_mb"] = round(get_rss_mb(), 1)

    # --- Save intermediates ---
    np.save(f"{out_chr_dir}/xsol.npy", xsol)
    np.save(f"{out_chr_dir}/new_haplotypes.npy", final_haplotypes)
    json.dump({
        "subject_name": subject_name,
        "subject_id": int(subject_id),
        "chromosome": chromosome,
        "capacity": capacity,
        "seed": seed,
        "max_utility_loss": float(max_utility_loss),
        "utility_loss": float(util_loss),
        "utility_loss_frac": float(util_loss / max_utility_loss),
        "pmi_gain": float(pmi_gain),
        "max_pmi_gain": float(np.sum(pmi_blocks)),
        "pmi_gain_frac": float(pmi_gain / np.sum(pmi_blocks)) if np.sum(pmi_blocks) > 0 else 0,
        "moves": int(np.sum(xsol)),
        "changed_alleles": int(changed),
        "total_alleles": int(total_valid),
        "changed_frac": float(changed / total_valid),
        "timings": timings,
    }, open(f"{out_chr_dir}/stats.json", "w"), indent=2)

    print(f"\n=== Done chr{chromosome} in {timings['total_s']}s | "
          f"peak RSS={timings['peak_rss_mb']:.0f}MB ===", flush=True)
    print(f"  load={timings['load_data_s']}s  "
          f"graph_prep={timings['graph_prep_s']}s  "
          f"optimizer={timings['optimizer_s']}s  "
          f"stacker={timings['stacker_s']}s  "
          f"vcf={timings['vcf_write_s']}s", flush=True)


if __name__ == "__main__":
    main()
