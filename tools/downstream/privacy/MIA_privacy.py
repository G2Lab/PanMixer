import numpy as np
import sys
import json
import os
from tools.common.slurm_helper import launch_job_multichromosome
from tools.common.utils import load_data, load_data_multichromosome, profile

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
    BASE_PATH,
)


def MIA_privacy_computer(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [experiment_number]
    return launch_job_multichromosome("MIA_privacy", args, memory="64g", cpus="1", num_tasks=str(num_tasks))


def binarize(x):
    """0=ref, >0 becomes 1, -1 stays -1."""
    return np.where(x > 0, 1, x)


def haplotype_match_scores(individual_haps, other_haps_3d):
    """Count matching haplotype alleles (exact match, ignoring missing).
    individual_haps: (n_sites, 2)
    other_haps_3d: (n_subjects, n_sites, 2)
    Returns: (n_subjects,) match counts, (n_subjects,) valid counts
    """
    ind = individual_haps[None, :, :]
    valid = (ind != -1) & (other_haps_3d != -1)
    matching = (ind == other_haps_3d) & valid
    return np.sum(matching, axis=(1, 2)), np.sum(valid, axis=(1, 2))


def genotype_match_scores(individual_haps, other_haps_3d):
    """Count matching genotypes (sum of two haplotypes per site).
    Returns: (n_subjects,) match counts, (n_subjects,) valid counts
    """
    ind_geno = np.sum(individual_haps, axis=1)
    other_geno = np.sum(other_haps_3d, axis=2)
    valid = (ind_geno[None, :] >= 0) & (other_geno >= 0)
    matching = (ind_geno[None, :] == other_geno) & valid
    return np.sum(matching, axis=1), np.sum(valid, axis=1)


@profile
def main():
    experiment_number = sys.argv[1]
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])

    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"

    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    pangenome = np.load(chromosome_path + "pangenome.npy")

    subject_index = np.where(subjects == subject_name)[0][0]

    # Load obfuscated haplotypes
    obfuscated = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")
    original = pangenome[subject_index]
    pg_others = np.delete(pangenome, subject_index, axis=0)
    pg_other_names = np.delete(subjects, subject_index)

    n_positions = pangenome.shape[1]

    # --- 1. Full allele comparison (all positions, multiallelic) ---
    hap_match_orig, hap_total_orig = haplotype_match_scores(obfuscated, original[None, :, :])
    hap_match_pg, hap_total_pg = haplotype_match_scores(obfuscated, pg_others)
    geno_match_orig, geno_total_orig = genotype_match_scores(obfuscated, original[None, :, :])
    geno_match_pg, geno_total_pg = genotype_match_scores(obfuscated, pg_others)

    # --- 2. Binarized comparison at 1000G-overlapping sites ---
    posref_mask_path = chromosome_path + "pangenome_mask_posref.npy"
    onek_masked_path = chromosome_path + "1000g_phased_masked_posref.npy"

    has_1000g = os.path.exists(posref_mask_path) and os.path.exists(onek_masked_path)

    if has_1000g:
        posref_mask = np.load(posref_mask_path)
        attack_db = np.load(onek_masked_path).astype(np.float32)

        obfuscated_bin = binarize(obfuscated[posref_mask]).astype(np.float32)
        original_bin = binarize(original[posref_mask]).astype(np.float32)
        pg_others_bin = binarize(pg_others[:, posref_mask]).astype(np.float32)

        hap_match_orig_bin, hap_total_orig_bin = haplotype_match_scores(obfuscated_bin, original_bin[None, :, :])
        hap_match_pg_bin, hap_total_pg_bin = haplotype_match_scores(obfuscated_bin, pg_others_bin)
        hap_match_1000g, hap_total_1000g = haplotype_match_scores(obfuscated_bin, attack_db)

        geno_match_orig_bin, geno_total_orig_bin = genotype_match_scores(obfuscated_bin, original_bin[None, :, :])
        geno_match_pg_bin, geno_total_pg_bin = genotype_match_scores(obfuscated_bin, pg_others_bin)
        geno_match_1000g, geno_total_1000g = genotype_match_scores(obfuscated_bin, attack_db)

    # --- Save per-chromosome results ---
    output_dir = EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}"
    os.makedirs(output_dir, exist_ok=True)

    # Save arrays for cross-chromosome aggregation
    np.savez(
        f"{output_dir}/MIA_privacy.npz",
        # Full allele comparison
        hap_match_orig=hap_match_orig,
        hap_total_orig=hap_total_orig,
        hap_match_pg=hap_match_pg,
        hap_total_pg=hap_total_pg,
        geno_match_orig=geno_match_orig,
        geno_total_orig=geno_total_orig,
        geno_match_pg=geno_match_pg,
        geno_total_pg=geno_total_pg,
        pg_other_names=pg_other_names,
        # Binarized at overlapping sites
        hap_match_orig_bin=hap_match_orig_bin if has_1000g else np.array([]),
        hap_total_orig_bin=hap_total_orig_bin if has_1000g else np.array([]),
        hap_match_pg_bin=hap_match_pg_bin if has_1000g else np.array([]),
        hap_total_pg_bin=hap_total_pg_bin if has_1000g else np.array([]),
        hap_match_1000g=hap_match_1000g if has_1000g else np.array([]),
        hap_total_1000g=hap_total_1000g if has_1000g else np.array([]),
        geno_match_orig_bin=geno_match_orig_bin if has_1000g else np.array([]),
        geno_total_orig_bin=geno_total_orig_bin if has_1000g else np.array([]),
        geno_match_pg_bin=geno_match_pg_bin if has_1000g else np.array([]),
        geno_total_pg_bin=geno_total_pg_bin if has_1000g else np.array([]),
        geno_match_1000g=geno_match_1000g if has_1000g else np.array([]),
        geno_total_1000g=geno_total_1000g if has_1000g else np.array([]),
    )

    # Human-readable summary for this chromosome
    hap_frac_orig = hap_match_orig[0] / hap_total_orig[0]
    best_pg = np.argmax(hap_match_pg / hap_total_pg)
    hap_frac_best_pg = hap_match_pg[best_pg] / hap_total_pg[best_pg]
    hap_rank = int(np.sum((hap_match_pg / hap_total_pg) >= hap_frac_orig)) + 1

    summary = {
        "subject": subject_name,
        "chromosome": chromosome,
        "n_positions": int(n_positions),
        "hap_match_frac_vs_original": float(hap_frac_orig),
        "hap_match_frac_vs_best_pangenome": float(hap_frac_best_pg),
        "hap_best_pangenome_subject": str(pg_other_names[best_pg]),
        "hap_original_rank": hap_rank,
        "hap_original_rank_out_of": len(pg_other_names) + 1,
    }

    if has_1000g:
        n_sites = len(posref_mask)
        best_1000g = np.argmax(hap_match_1000g / hap_total_1000g)
        summary["n_overlapping_sites"] = int(n_sites)
        summary["hap_bin_match_frac_vs_original"] = float(hap_match_orig_bin[0] / hap_total_orig_bin[0])
        summary["hap_bin_match_frac_vs_best_1000g"] = float(hap_match_1000g[best_1000g] / hap_total_1000g[best_1000g])
        summary["hap_bin_match_frac_vs_best_pangenome"] = float(hap_match_pg_bin[np.argmax(hap_match_pg_bin / hap_total_pg_bin)] / hap_total_pg_bin[np.argmax(hap_match_pg_bin / hap_total_pg_bin)])

    with open(f"{output_dir}/MIA_privacy.json", "w") as f:
        json.dump(summary, f, indent=4)

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
