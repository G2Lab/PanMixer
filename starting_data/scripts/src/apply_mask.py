
import sys
import numpy as np


def apply_mask(chromosome):
    thousand_g = np.load(f"/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr{chromosome}/1000g_phased.npy")
    thousand_g_phased_mask = np.load(f"/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr{chromosome}/thousand_g_phased_mask.npy")

    thousand_g_phased_masked = thousand_g[:, thousand_g_phased_mask]
    np.save(f"/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/chr{chromosome}/1000g_phased_masked.npy", thousand_g_phased_masked)
    print(f"Masked thousand_g_phased saved to chr{chromosome}/1000g_phased_masked.npy")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python apply_mask.py <chromosome>")
        sys.exit(1)

    chromosome = sys.argv[1]
    apply_mask(chromosome)