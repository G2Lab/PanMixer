"""
Compute gap score for the fully-sampled haplotypes (before stacker selection).
Tests whether full replacement with sampled haplotypes breaks linkability.
"""
import numpy as np
import json
import os
import sys

from constants import STARTING_DATA_PATH, BASE_PATH, ONEK_MASKED_NPY
from tools.diploid_gap_score import score_genotypes_remove_missing, score_haplotypes


def score_sampled_haplotype(chromosome, subject_id):
    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"

    subjects = np.load(chromosome_path + "pangenome_subjects.npy")
    haplotypes = np.load(chromosome_path + "pangenome.npy")

    # Load attack DB
    posref_mask_path = f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask_posref.npy"
    posref_onek_path = chromosome_path + "1000g_phased_masked_posref.npy"

    if os.path.exists(posref_mask_path) and os.path.exists(posref_onek_path):
        pangenome_mask = np.load(posref_mask_path)
        attack_db_masked = np.load(posref_onek_path)
    else:
        pangenome_mask = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy")
        attack_db_masked = np.load(chromosome_path + ONEK_MASKED_NPY)

    # Load the fully sampled haplotypes
    sampled = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/new_haplotypes.npy")

    # Compute alt AF from 1000G
    n_haplotypes_1000g = attack_db_masked.shape[0] * 2
    alt_counts = np.sum(attack_db_masked, axis=(0, 2))
    alt_af = alt_counts / n_haplotypes_1000g

    polymorphic = (alt_af > 0) & (alt_af < 1)

    attack_db_haplotypes = attack_db_masked[:, polymorphic]
    alt_af = alt_af[polymorphic]

    # Binarize sampled haplotypes
    sampled_masked = sampled[pangenome_mask][polymorphic]
    sampled_masked = np.where(sampled_masked > 0, 1, sampled_masked).astype(np.float32)

    # Original haplotypes
    original = haplotypes[subject_id]
    original = np.expand_dims(original, axis=0)
    original = original[:, pangenome_mask][:, polymorphic]
    original = np.where(original > 0, 1, original).astype(np.float32)

    attack_db_haplotypes = attack_db_haplotypes.astype(np.float32)

    N_ATTACK = attack_db_haplotypes.shape[0]
    NUM_CHUNKS = 10
    chunks = np.array_split(np.arange(N_ATTACK), NUM_CHUNKS)

    # Genotype scoring
    scores_genotypes = np.zeros(N_ATTACK)
    for i in range(NUM_CHUNKS):
        scores_genotypes[chunks[i]], _ = score_genotypes_remove_missing(sampled_masked, attack_db_haplotypes[chunks[i]], alt_af)

    self_score, _ = score_genotypes_remove_missing(sampled_masked, original, alt_af)
    gap = self_score[0] - np.max(scores_genotypes)

    return {
        "subject": str(subjects[subject_id]),
        "subject_id": int(subject_id),
        "chromosome": int(chromosome),
        "self_score": float(self_score[0]),
        "best_attack_score": float(np.max(scores_genotypes)),
        "genotype_gap_score": float(gap),
    }


if __name__ == "__main__":
    chromosomes = range(1, 23) if len(sys.argv) < 2 else [int(sys.argv[1])]
    subject_ids = range(44) if len(sys.argv) < 3 else [int(sys.argv[2])]

    results = []
    for chrom in chromosomes:
        for sid in subject_ids:
            print(f"Scoring chr{chrom} subject {sid}...")
            r = score_sampled_haplotype(chrom, sid)
            print(f"  gap_score={r['genotype_gap_score']:.2f}")
            results.append(r)

    out_path = f"{BASE_PATH}/sampled_haplotype_gap_scores.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")

    negatives = [r for r in results if r["genotype_gap_score"] < 0]
    print(f"\n{len(negatives)}/{len(results)} have gap score < 0")
