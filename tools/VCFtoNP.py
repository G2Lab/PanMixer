import numpy as np
import os
import gzip
import sys

def get_numpy_matrices(vcf_file_path, phased=False):
    #verify file_path exists
    if not os.path.exists(vcf_file_path):
        raise FileNotFoundError("File not found: " + vcf_file_path)

    num_subjects = 0
    subjects_in_vcf = []
    num_sites = 0


    open_func = gzip.open if vcf_file_path.endswith('.gz') else open
    with open_func(vcf_file_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                num_subjects = len(line.strip().split('\t')[9:])
                subjects_in_vcf = line.strip().split('\t')[9:]
                continue
            if line[0] == '#':
                continue
            num_sites += 1
    genotypes = np.zeros((num_subjects, num_sites, 2), dtype=np.int16)

    positions = []

    with open_func(vcf_file_path, 'rt') as vcf:
        site = 0
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            positions.append(int(fields[1]))
            for i, subject in enumerate(subjects_in_vcf):
                if phased:
                    genotype = fields[9+i].split(':')[0].split('|')
                else:
                    genotype = fields[9+i].split(':')[0].split('/')

                if genotype[0] == '.':
                    genotypes[i, site, 0] = -1
                else:
                    genotypes[i, site, 0] = int(genotype[0])
                if len(genotype) == 1:
                    genotypes[i, site, 1] = -1
                    continue
                if genotype[1] == '.':
                    genotypes[i, site, 1] = -1
                else:
                    genotypes[i, site, 1] = int(genotype[1])
            site += 1
    return genotypes, subjects_in_vcf, positions

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 VCFtoNP.py <vcf_file_path> <phased?>")
        sys.exit(1)

    vcf_file_path = sys.argv[1]
    phased = sys.argv[2].lower() == 'true'
    genotypes, subjects, positions = get_numpy_matrices(vcf_file_path, phased)

    # Define output file names, handling both `.vcf` and `.vcf.gz`
    base_name = vcf_file_path.replace('.vcf.gz', '').replace('.vcf', '')

    np.save(f"{base_name}.npy", genotypes)
    np.save(f"{base_name}_subjects.npy", subjects)
    np.save(f"{base_name}_positions.npy", positions)

    print("Done")
