
import sys
import os
import numpy as np

# Derive starting_data path relative to this script's location
# (starting_data/scripts/src/ -> starting_data/)
STARTING_DATA_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def apply_mask(chromosome):
    thousand_g = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/1000g_phased.npy")
    thousand_g_phased_mask = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/thousand_g_phased_mask.npy")

    thousand_g_phased_masked = thousand_g[:, thousand_g_phased_mask]
    np.save(f"{STARTING_DATA_PATH}/chr{chromosome}/1000g_phased_masked.npy", thousand_g_phased_masked)
    print(f"Masked thousand_g_phased saved to chr{chromosome}/1000g_phased_masked.npy")

    # Relaxed (pos, ref) mask for gap score — binarize 1000G alleles (0=ref, >0→1)
    posref_mask_path = f"{STARTING_DATA_PATH}/chr{chromosome}/thousand_g_phased_mask_posref.npy"
    if os.path.exists(posref_mask_path):
        posref_mask = np.load(posref_mask_path)
        thousand_g_posref = thousand_g[:, posref_mask]
        # Binarize: any non-zero allele becomes 1 (ref vs alt)
        thousand_g_posref = np.clip(thousand_g_posref, 0, 1)
        np.save(f"{STARTING_DATA_PATH}/chr{chromosome}/1000g_phased_masked_posref.npy", thousand_g_posref)
        print(f"Relaxed masked thousand_g_phased saved to chr{chromosome}/1000g_phased_masked_posref.npy ({thousand_g_posref.shape[1]} sites)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python apply_mask.py <chromosome>")
        sys.exit(1)

    chromosome = sys.argv[1]
    apply_mask(chromosome)