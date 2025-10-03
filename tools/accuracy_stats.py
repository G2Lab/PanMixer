import numpy as np
from tools.utils import load_data_multichromosome, load_data
from tools.VCFtoNP import get_numpy_matrices
import json
import sys
from tools.slurm_helper import launch_job

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
)

def align_beagle_to_positions(beagle_positions, beagle_haplotypes, target_positions, fill_value=-1):
    # Create mapping from Beagle positions to index
    beagle_pos_to_idx = {pos: i for i, pos in enumerate(beagle_positions)}

    matched_indices = []
    missing_count = 0
    for pos in target_positions:
        if pos in beagle_pos_to_idx:
            matched_indices.append(beagle_pos_to_idx[pos])
        else:
            matched_indices.append(None)
            missing_count += 1

    if missing_count > 0:
        print(f"⚠️ Warning: {missing_count} out of {len(target_positions)} target positions were not found in Beagle output.")

    # Build aligned array with placeholder for missing
    aligned = np.full((len(target_positions), 2), fill_value, dtype=beagle_haplotypes.dtype)

    for i, idx in enumerate(matched_indices):
        if idx is not None:
            aligned[i] = beagle_haplotypes[idx]

    return aligned

def get_all_variant_subsets(chromosome):

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
        "all": np.ones_like(snp_mask, dtype=bool),
        "complex": ~snp_mask,
        "snp_only": snp_mask,
        "snp_only_maf_0.01": snps_maf_0_01,
        "snp_only_maf_0.05": snps_maf_0_05,
        "snp_only_maf_0.1": snps_maf_0_1,
        "snp_only_maf_0.5": snps_maf_0_5,
    }
    return subsets

def compute_accuracy(original_haplotypes, obfuscated_haplotypes, beagle_haplotypes, subset, allele_frequency):
    assert original_haplotypes.shape == obfuscated_haplotypes.shape == beagle_haplotypes.shape, \
        f"Haplotype shapes do not match! Original: {original_haplotypes.shape}, Obfuscated: {obfuscated_haplotypes.shape}, Beagle: {beagle_haplotypes.shape}"
    
    # Compute the number of SNPs in the subset
    print(original_haplotypes[subset].shape)
    print(obfuscated_haplotypes[subset].shape)
    print(beagle_haplotypes[subset].shape)

    print(original_haplotypes[subset])
    print(obfuscated_haplotypes[subset])
    print(beagle_haplotypes[subset])
    accu_beagle = np.mean(original_haplotypes[subset] == beagle_haplotypes[subset])
    accu_obfuscated = np.mean(original_haplotypes[subset] == obfuscated_haplotypes[subset])

    non_missing_original = original_haplotypes[subset] != -1
    non_missing_obfuscated = obfuscated_haplotypes[subset] != -1
    non_missing_beagle = beagle_haplotypes[subset] != -1

    print("non_missing_original", non_missing_original.shape)
    print("non_missing_obfuscated", non_missing_obfuscated.shape)
    print("non_missing_beagle", non_missing_beagle.shape)

    non_missing_all = np.logical_and(np.logical_and(non_missing_original, non_missing_beagle), non_missing_obfuscated)

    obfuscated_different = obfuscated_haplotypes[subset][non_missing_all] != original_haplotypes[subset][non_missing_all]

    reconstruction_accuracy = np.mean(original_haplotypes[subset][non_missing_all][obfuscated_different] == beagle_haplotypes[subset][non_missing_all][obfuscated_different])
    print("reconstruction accuracy", reconstruction_accuracy)

    print("allele_freq_shape",allele_frequency.shape)
    major_allele = np.argmax(allele_frequency, axis=1)
    major_allele_stacked = np.column_stack((major_allele, major_allele))
    print(major_allele_stacked.shape)
    major_allele_accuracy = np.mean(original_haplotypes[subset] == major_allele_stacked[subset])
    # Compute the accuracy of the major allele

    return {
        "obfuscated_accuracy": accu_obfuscated,
        "beagle_accuracy": accu_beagle,
        "major_allele_accuracy": major_allele_accuracy,
        "reconstruction_accuracy": reconstruction_accuracy,
    }

def accuracy_computer(experiment_number, chromosome, row_id):
    data_df = load_data(experiment_number, chromosome)

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"

    original_subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    original_positions = np.load(chromosome_path + "pangenome_positions.npy")
    original_haplotypes = np.load(chromosome_path + "pangenome.npy")

    subsets = get_all_variant_subsets(chromosome)

    row = data_df.iloc[row_id]
    subject = row["subject"]

    blocks_for_pop = json.load(open(chromosome_path + "blocks_dict.json"))
    keys = list(blocks_for_pop.keys())

    obfuscated_haplotypes = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")

    xsol = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/xsol.npy")

    true_indices = np.argwhere(xsol == 1)

    xsol_expanded = np.zeros_like(obfuscated_haplotypes, dtype=bool)

    allele_frequency = np.load(chromosome_path + "allele_frequencies.npy")

    for i,j in true_indices:
        block_indices = blocks_for_pop[keys[i]][1]
        for block_index in block_indices:
            xsol_expanded[block_index, j] = True

    subsets_for_just_obfuscated = {}
    for subset_name, subset in subsets.items():
        subset_stacked = np.column_stack((subset, subset))
        subsets_for_just_obfuscated[subset_name] = np.logical_and(subset_stacked, xsol_expanded)
        print(np.sum(subsets_for_just_obfuscated[subset_name]), subset_name)

    subject_index = np.where(original_subjects == subject)[0][0]
    original_haplotypes = original_haplotypes[subject_index]

    beagle_genotypes, beagle_subjects, beagle_positions = get_numpy_matrices(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/beagle_refined.vcf.gz.vcf.gz", phased=True)

    beagle_subject_index = list(beagle_subjects).index(subject)
    beagle_haplotypes = beagle_genotypes[beagle_subject_index]    

    # Align Beagle haplotypes to original haplotypes
    beagle_haplotypes = align_beagle_to_positions(beagle_positions, beagle_haplotypes, original_positions)

    results = {}
    for subset_name, subset in subsets_for_just_obfuscated.items():
        results[subset_name] = compute_accuracy(original_haplotypes, obfuscated_haplotypes, beagle_haplotypes, subset, allele_frequency)

    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/accuracy_results.json", "w") as file:
        json.dump(results, file)
    
if __name__ == "__main__":
    experiment_number = sys.argv[1]
    chromosome = sys.argv[2]
    row_id = int(sys.argv[3])

    accuracy_computer(experiment_number, chromosome, row_id)

def accuracy_stats(experiment_number):
    data_df = load_data(experiment_number, 21)
    num_tasks = len(data_df)
    launch_job("accuracy_stats", [experiment_number, 21], memory="32g", cpus="1", num_tasks=str(num_tasks))