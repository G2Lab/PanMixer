import numpy as np
from tools.utils import load_data_multichromosome, load_data
import json
import sys
from tools.slurm_helper import launch_job_multichromosome

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
)

def get_all_af_subsets(chromosome):

    snp_indices = np.load(STARTING_DATA_PATH + f"/chr{chromosome}/pangenome_mask.npy")

    global_af = np.load(STARTING_DATA_PATH + f"/chr{chromosome}/allele_frequencies.npy")
    
    snp_mask = np.zeros(global_af.shape[0], dtype=bool)
    snp_mask[snp_indices] = True
    
    print(snp_mask.shape, global_af.shape)
    # Compute minor allele frequency
    global_maf = np.zeros(global_af.shape[0])
    for i in range(global_af.shape[0]):
        if np.any(np.isnan(global_af[i])):
            global_maf[i] = 1
            continue
        mininum_non_zero = np.min(global_af[i][global_af[i] > 0])
        maximum_non_zero = np.max(global_af[i][global_af[i] > 0])
        if mininum_non_zero < 0.5:
            global_maf[i] = mininum_non_zero
        else:
            global_maf[i] = maximum_non_zero

    # Mask non-SNP positions by setting their MAF to 1 (out of consideration for subsets)
    global_maf[~snp_mask] = 1

    # Define subsets based on MAF thresholds
    snps_maf_0_01 = global_maf < 0.01
    snps_maf_0_05 = (global_maf >= 0.01) & (global_maf < 0.05)
    snps_maf_0_1 = (global_maf >= 0.05) & (global_maf < 0.1)
    snps_maf_0_5 = (global_maf >= 0.1) & (global_maf <= 0.5)  # Ensure 1 is ignored

    # Sanity check: all SNPs must belong to one of the MAF subsets
    #assert np.equal(subset_indices, snps_maf_0_01 | snps_maf_0_05 | snps_maf_0_1 | snps_maf_0_5).all(), \
    #    "Not all SNPs are covered by the MAF subsets!"

    # Sanity check: subsets must be disjoint
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_01 & snps_maf_0_05).all(), \
        "MAF subsets are not disjoint (0.01 and 0.05 overlap)!"
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_01 & snps_maf_0_1).all(), \
        "MAF subsets are not disjoint (0.01 and 0.1 overlap)!"
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_01 & snps_maf_0_5).all(), \
        "MAF subsets are not disjoint (0.01 and 0.5 overlap)!"
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_05 & snps_maf_0_1).all(), \
        "MAF subsets are not disjoint (0.05 and 0.1 overlap)!"
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_05 & snps_maf_0_5).all(), \
        "MAF subsets are not disjoint (0.05 and 0.5 overlap)!"
    assert np.equal(np.zeros_like(snp_mask), snps_maf_0_1 & snps_maf_0_5).all(), \
        "MAF subsets are not disjoint (0.1 and 0.5 overlap)!"

    subsets = {
        "complex": ~snp_mask,
        "snp_only": snp_mask,
        "snp_only_maf_0.01": snps_maf_0_01,
        "snp_only_maf_0.05": snps_maf_0_05,
        "snp_only_maf_0.1": snps_maf_0_1,
        "snp_only_maf_0.5": snps_maf_0_5,
    }
    return subsets

def compute_l1_loss_af_complex(original_haplotypes, new_haplotypes):
    """
    Computes the L1 loss for allele frequencies (AF) between the original and new haplotypes,
    including support for multiallelic SNPs.

    Parameters:
        original_haplotypes (np.ndarray): Original haplotypes (n_variants x n_haplotypes).
        new_haplotypes (np.ndarray): New haplotypes (n_variants x n_haplotypes).

    Returns:
        float: The L1 loss for allele frequencies.
    """
    # Ensure the same number of variants
    assert original_haplotypes.shape == new_haplotypes.shape, "Shapes must match"
    
    # Initialize total L1 loss
    l1_loss = np.zeros(original_haplotypes.shape[1])

    count = np.zeros(original_haplotypes.shape[1])
    
    # Loop over each SNP (row in haplotypes matrix)
    for variant_idx in range(original_haplotypes.shape[1]):
        # Extract alleles for the current SNP
        original_alleles = original_haplotypes[:, variant_idx]
        new_alleles = new_haplotypes[:, variant_idx]
        
        # Get unique alleles and counts for original and new haplotypes
        unique_original, counts_original = np.unique(original_alleles, return_counts=True)
        unique_new, counts_new = np.unique(new_alleles, return_counts=True)
        
        # Normalize counts to allele frequencies
        total_original = counts_original.sum()
        total_new = counts_new.sum()
        
        af_original = dict(zip(unique_original, counts_original / total_original))
        af_new = dict(zip(unique_new, counts_new / total_new))
        
        # Compute L1 loss for this SNP
        all_alleles = set(af_original.keys()).union(set(af_new.keys()))
        for allele in all_alleles:
            af_orig = af_original.get(allele, 0.0)
            af_new_val = af_new.get(allele, 0.0)
            l1_loss[variant_idx] += abs(af_orig - af_new_val)
            count[variant_idx] += 1     
    
    return l1_loss, count


def af_loss_computer(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])
    
    args = [experiment_number]
    launch_job_multichromosome("af_loss", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

def main():
    experiment_number = sys.argv[1]
    chromosome = sys.argv[2]
    row_id = int(sys.argv[3])

    get_af_loss(experiment_number, chromosome, row_id)

def get_af_loss(experiment_number, chromosome, row_id):
    data_df = load_data(experiment_number, chromosome)

    original_subjects = np.load(STARTING_DATA_PATH + f"/chr{chromosome}/pangenome_subjects.npy")
    original_haplotypes = np.load(STARTING_DATA_PATH + f"/chr{chromosome}/pangenome.npy")

    num_tasks = len(data_df)

    subsets_af = get_all_af_subsets(chromosome)

    print(f"Processing row {row_id} of {num_tasks}")
    row = data_df.iloc[row_id]
    subject = row["subject"]
    capacity = row["capacity"]

    interested_variants = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")

    subject_index = np.where(original_subjects == subject)[0][0]
    subject_haplotype = original_haplotypes[subject_index]

    new_pangenome = original_haplotypes.copy()
    new_pangenome[subject_index] = interested_variants


    all_loss,all_count = compute_l1_loss_af_complex(original_haplotypes, new_pangenome)
    results = {}

    results["all"] = float(np.sum(all_loss) / np.sum(all_count))

    for key, value in subsets_af.items():
        subset_all_loss = all_loss[value]
        subset_all_count = all_count[value]
        results[key] = float(np.sum(subset_all_loss) / np.sum(subset_all_count))

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/af_loss.json", "w") as f:
        json.dump(results, f, indent=4)

if __name__ == "__main__":
    main()