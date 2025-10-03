import numpy as np
import json
from tools.utils import get_chromosome_path
import pickle

def get_accuracy_alignments_analysis(starting_data_path, chromosome):
    path_with_chromosome = get_chromosome_path(starting_data_path, chromosome)

    thousand_g_alignments = np.load(path_with_chromosome + "1000g.npy", allow_pickle=True)
    thousand_g_alignments_positions = np.load(path_with_chromosome + "1000g_positions.npy", allow_pickle=True)
    thousand_g_subjects = np.load(path_with_chromosome + "1000g_subjects.npy", allow_pickle=True)

    thousand_g_snps_only = np.load(path_with_chromosome + "1000g_snps.npy", allow_pickle=True)
    thousand_g_snps_only_positions = np.load(path_with_chromosome + "1000g_snps_positions.npy", allow_pickle=True)
    thousand_g_snps_only_subjects = np.load(path_with_chromosome + "1000g_snps_subjects.npy", allow_pickle=True)

    # Get the shared subjects in a sorted, stable order
    subjects_intersection = np.intersect1d(thousand_g_subjects, thousand_g_snps_only_subjects)

    # Create a mapping from subject to index for both arrays
    subject_to_index_full = {subject: idx for idx, subject in enumerate(thousand_g_subjects)}
    subject_to_index_snps = {subject: idx for idx, subject in enumerate(thousand_g_snps_only_subjects)}

    # Now, for each subject in the intersection (in order), get their corresponding indices
    subject_1000g_full_idx = np.array([subject_to_index_full[subject] for subject in subjects_intersection])
    subject_1000g_snps_only_idx = np.array([subject_to_index_snps[subject] for subject in subjects_intersection])

    assert np.array_equal(thousand_g_subjects[subject_1000g_full_idx], thousand_g_snps_only_subjects[subject_1000g_snps_only_idx]), "Subjects do not match!"

    # Filter the alignments based on the shared subjects
    # 1. Find shared SNP positions in a sorted order
    shared_positions = np.intersect1d(thousand_g_alignments_positions, thousand_g_snps_only_positions)

    # 2. Create mapping from position to index for both arrays
    pos_to_idx_alignments = {pos: idx for idx, pos in enumerate(thousand_g_alignments_positions)}
    pos_to_idx_snps = {pos: idx for idx, pos in enumerate(thousand_g_snps_only_positions)}

    # 3. Get column indices in the same order for both datasets
    alignment_pos_idx = np.array([pos_to_idx_alignments[pos] for pos in shared_positions])
    snps_only_pos_idx = np.array([pos_to_idx_snps[pos] for pos in shared_positions])

    # 4. Apply subject and SNP filtering
    thousand_g_alignments_filtered = thousand_g_alignments[subject_1000g_full_idx][:, alignment_pos_idx]
    thousand_g_snps_only_filtered = thousand_g_snps_only[subject_1000g_snps_only_idx][:, snps_only_pos_idx]

    # 5. Check shapes
    print(f"Filtered alignments shape: {thousand_g_alignments_filtered.shape}")
    print(f"Filtered SNPs shape: {thousand_g_snps_only_filtered.shape}")

    accuracy_all = np.mean(thousand_g_alignments_filtered == thousand_g_snps_only_filtered, axis=0)
    #plot histogram of accuracy_all
    import matplotlib.pyplot as plt
    plt.hist(accuracy_all, bins=100)
    plt.xlabel("Accuracy")
    plt.ylabel("Frequency")
    plt.title("Histogram of Accuracy")

    plt.savefig("/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer4/plots/" + "accuracy_histogram.png")

    plt.close()

    accuracy_mean = np.mean(accuracy_all)
    print("Mean accuracy:", accuracy_mean)

    print("Accuracy without missing")

    missing_value = -1  # or whatever your placeholder is

    valid_mask = (thousand_g_alignments_filtered != missing_value) & (thousand_g_snps_only_filtered != missing_value)
    matches = thousand_g_alignments_filtered == thousand_g_snps_only_filtered

    # Set invalid comparisons to False or mask them out
    matches = matches & valid_mask

    # Count correct matches and total valid comparisons
    correct_counts = np.sum(matches, axis=0)

    valid_counts = np.sum(valid_mask, axis=0)

    valid_valid_counts_mask = valid_counts != 0

    accuracy_all = correct_counts[valid_valid_counts_mask] / valid_counts[valid_valid_counts_mask]

    np.save(path_with_chromosome + "accuracy_all.npy", accuracy_all)
    print(shared_positions.shape)
    print(valid_valid_counts_mask.shape)
    positions_filtered = shared_positions[valid_valid_counts_mask]
    np.save(path_with_chromosome + "accuracy_positions.npy", positions_filtered)

    #plot histogram of accuracy_all
    plt.hist(accuracy_all, bins=100)
    plt.xlabel("Accuracy")
    plt.ylabel("Frequency")
    plt.title("Histogram of Accuracy without missing")
    plt.savefig("/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer4/plots/" + "accuracy_histogram_no_missing.png")
    plt.close()

    accuracy_mean = np.mean(accuracy_all)
    print("Mean accuracy without missing:", accuracy_mean)




    



