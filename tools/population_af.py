import numpy as np
import pickle
from tools.utils import get_chromosome_path

def compute_allele_counts_and_missing(haplotypes, num_alleles):

    # Count the number of missing alleles at each position.
    missing_counts = np.sum(haplotypes == -1, axis=(0, 2))

    # Calculate the allele counts.
    # The allele counts are stored in a 3D array with shape (num_positions, num_subjects, num_alleles).
    # The last dimension has size 2 because we are working with diploid genomes.
    allele_counts = np.zeros((haplotypes.shape[1], num_alleles), dtype=int)

    # Loop over positions.
    for pos in range(haplotypes.shape[1]):
        # Extract all alleles at this position (subjects x 2 alleles).
        alleles_at_pos = haplotypes[:, pos, :]

        # Filter out missing alleles so that np.bincount works correctly.
        valid_alleles = alleles_at_pos[alleles_at_pos >= 0]

        # Count the occurrences of each valid allele.
        # Using minlength ensures the output has the length of 2.
        counts = np.bincount(valid_alleles)
        counts_padded = np.pad(counts, (0, num_alleles - len(counts)), mode="constant")
        assert np.sum(counts) == 2 * haplotypes.shape[0] - missing_counts[pos]

        # Store the allele counts for this position.
        allele_counts[pos, :] = counts_padded

    allele_freqs = allele_counts / (2 * haplotypes.shape[0] - missing_counts[:, np.newaxis])

    return allele_counts, missing_counts, allele_freqs

def get_population_af(starting_data_path, chromosome):
    path_with_chromosome = get_chromosome_path(starting_data_path, chromosome)
    
    alignments = np.load(path_with_chromosome + "1000g.npy")
    alignments = alignments.astype(int)
    alignment_positions = np.load(path_with_chromosome + "1000g_positions.npy")

    pangenome = np.load(path_with_chromosome + "pangenome.npy")

    assert alignments.shape[1] == len(alignment_positions)
    """
    AFs_counts = {}
    AFs_freqs = {}
    for i,pos in enumerate(alignment_positions):
        unique_elements, counts_elements = np.unique(alignments[:,i], return_counts=True)
        AFs_counts[pos] = dict(zip(unique_elements, counts_elements))
        AFs_freqs[pos] = dict(zip(unique_elements, counts_elements/np.sum(counts_elements)))
    
    with open(path_with_chromosome + "population_af_counts.pickle", "wb") as file:
        pickle.dump(AFs_counts, file)
    with open(path_with_chromosome + "population_af_freqs.pickle", "wb") as file:
        pickle.dump(AFs_freqs, file)
    """
    
    print(f"Population allele frequencies saved for chromosome {chromosome}")
    maximum_number_of_alleles = max(np.max(pangenome), np.max(alignments)) + 1
    print(f"Maximum number of alleles at a position: {maximum_number_of_alleles}")

    allele_counts, missing_counts, allele_frequencies = compute_allele_counts_and_missing(alignments, maximum_number_of_alleles)

    # Save the allele frequencies.
    np.save(path_with_chromosome + "population_af.npy", allele_frequencies)
    print(f"Population allele frequencies saved for chromosome {chromosome}")
    np.save(path_with_chromosome + "population_af_missing.npy", missing_counts)
    print(f"Population allele frequencies missing counts saved for chromosome {chromosome}")
    np.save(path_with_chromosome + "population_af_counts.npy", allele_counts)
    print(f"Population allele frequencies counts saved for chromosome {chromosome}")


    assert pangenome.shape[1] == alignments.shape[1]
    assert pangenome.shape[2] == alignments.shape[2]
    stacked = np.concatenate([pangenome, alignments], axis=0)

    allele_counts, missing_counts, allele_frequencies = compute_allele_counts_and_missing(stacked, maximum_number_of_alleles)

    # Save the allele frequencies.
    np.save(path_with_chromosome + "population_af_stacked.npy", allele_frequencies)
    print(f"Population allele frequencies saved for chromosome {chromosome}")
    np.save(path_with_chromosome + "population_af_missing_stacked.npy", missing_counts)
    print(f"Population allele frequencies missing counts saved for chromosome {chromosome}")
    np.save(path_with_chromosome + "population_af_counts_stacked.npy", allele_counts)
    print(f"Population allele frequencies counts saved for chromosome {chromosome}")
