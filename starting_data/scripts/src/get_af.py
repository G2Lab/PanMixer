import numpy as np
import pickle
import sys

base_path = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/"

def get_population_af(chromosome):
    possible_alleles = np.load(base_path + f"chr{chromosome}/num_alleles.npy")

    pangenome_variants = np.load(base_path + f"chr{chromosome}/pangenome.npy")
    pangenome_positions = np.load(base_path + f"chr{chromosome}/pangenome_positions.npy")

    thousand_g_variants = np.load(base_path + f"chr{chromosome}/1000g_phased.npy")
    thousand_g_alignments_variants = np.load(base_path + f"chr{chromosome}/1000g.npy")

    mappings_pangenome_to_thousand_g_phased = pickle.load(open(base_path + f"chr{chromosome}/pangenome_to_thousand_g_phased.pickle", "rb"))
    mappings_pangenome_to_thousand_g_alignments = pickle.load(open(base_path + f"chr{chromosome}/pangenome_to_thousand_g_alignments.pickle", "rb"))

    max_allele = np.max([np.max(pangenome_variants), np.max(thousand_g_variants), np.max(thousand_g_alignments_variants)])

    print(f"Max allele: {max_allele}")

    allele_counts = np.zeros((len(pangenome_positions), max_allele + 1), dtype=int)

    for i in range(len(pangenome_positions)):
        """
        The flow of which dataset to use is as follows:
        If the position (ref and alt alleles) is in 1000g phased, use that.
        else if the position (ref and alt alleles) is in 1000g alignments, use that and the pangenome.
        else use the pangenome.

        This is because the 1000g phased dataset is the most reliable, and the 1000g alignments dataset is the second most reliable.
        The pangenome is the least reliable and has the fewest samples.
        """
        if i in mappings_pangenome_to_thousand_g_phased:
            valid_alleles = np.where(thousand_g_variants[:, mappings_pangenome_to_thousand_g_phased[i], :] != -1)
            for j, hap in zip(valid_alleles[0], valid_alleles[1]):
                allele_counts[i, thousand_g_variants[j, mappings_pangenome_to_thousand_g_phased[i], hap]] += 1
        elif i in mappings_pangenome_to_thousand_g_alignments:
            valid_alleles_alignments = np.where(thousand_g_alignments_variants[:, mappings_pangenome_to_thousand_g_alignments[i], :] != -1)
            #print(valid_alleles_alignments)
            valid_alleles_pangenome = np.where(pangenome_variants[:, i, :] != -1)

            for j, hap in zip(valid_alleles_alignments[0], valid_alleles_alignments[1]):
                allele_counts[i, thousand_g_alignments_variants[j, mappings_pangenome_to_thousand_g_alignments[i], hap]] += 1
            for j, hap in zip(valid_alleles_pangenome[0], valid_alleles_pangenome[1]):
                allele_counts[i, pangenome_variants[j, i, hap]] += 1
        else:
            valid_alleles_pangenome = np.where(pangenome_variants[:, i, :] != -1)
            for j, hap in zip(valid_alleles_pangenome[0], valid_alleles_pangenome[1]):
                allele_counts[i, pangenome_variants[j, i, hap]] += 1

        assert np.all(allele_counts[i, possible_alleles[i]:] == 0), f"Allele counts for position {i} are not zero for alleles greater than the possible alleles. \
             Allele counts: {allele_counts[i]}, Was in 1000g phased: {i in mappings_pangenome_to_thousand_g_phased}, \
                Was in 1000g alignments: {i in mappings_pangenome_to_thousand_g_alignments}, \
                Pangenome allele: {pangenome_variants[:,i]}, \
                Thousand g phased allele: {thousand_g_variants[:,mappings_pangenome_to_thousand_g_phased[i]] if i in mappings_pangenome_to_thousand_g_phased else 'N/A'}, \
                Thousand g alignments allele: {thousand_g_alignments_variants[:,mappings_pangenome_to_thousand_g_alignments[i]] if i in mappings_pangenome_to_thousand_g_alignments else 'N/A'}"
    # Normalize allele counts to get allele frequencies

    allele_frequencies = allele_counts / np.sum(allele_counts, axis=1, keepdims=True)
    # Save allele frequencies
    np.save(base_path + f"chr{chromosome}/allele_frequencies.npy", allele_frequencies)
    print(f"Allele frequencies saved for chromosome {chromosome}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_af.py <chromosome>")
        sys.exit(1)
    
    chromosome = sys.argv[1]
    get_population_af(chromosome)