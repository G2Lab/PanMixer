import numpy as np
import sys
import pickle

def get_possible_alleles(vcf_file):
    possible_alleles = {}
    index = 0
    with open(vcf_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "#":
                continue
            else:
                line = line.split("\t")
                pos = int(line[1])
                ref = line[3]
                alt = line[4]
                num_alleles = len(alt.split(",")) + 1

                possible_alleles[pos] = num_alleles
                index += 1
    
    vcf_file_path_without_extension = vcf_file.split(".")[:-1]
    vcf_file_path_without_extension = ".".join(vcf_file_path_without_extension)
    with open(f"{vcf_file_path_without_extension}_possible_alleles.pkl", "wb") as f:
        pickle.dump(possible_alleles, f)
    print(f"Saved possible alleles to {vcf_file_path_without_extension}_possible_alleles.pkl")
                
if __name__ == "__main__":
    get_possible_alleles(sys.argv[1])