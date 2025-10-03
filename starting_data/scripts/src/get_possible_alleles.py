import sys

import numpy as np
import gzip

from constants import BASE_PATH

def get_possible_alleles(chromosome):
    num_alleles = []
    pangenome_path = f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome.vcf.gz"
    with gzip.open(pangenome_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            pos = int(fields[1])
            ref_allele = fields[3]
            alt_alleles = fields[4].split(',')
            num_alleles.append(len(alt_alleles) + 1)
    num_alleles = np.array(num_alleles)
    np.save(f"{BASE_PATH}/starting_data/chr{chromosome}/num_alleles.npy", num_alleles)
    print(f"Num alleles saved to chr{chromosome}/num_alleles.npy")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_possible_alleles.py <chromosome>")
        sys.exit(1)

    chromosome = sys.argv[1]
    get_possible_alleles(chromosome)
