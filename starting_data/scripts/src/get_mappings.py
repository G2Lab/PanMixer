import gzip
import sys
import pickle
import numpy as np

import yaml
import os

def load_config():
    # Load default template
    with open("../config.yaml") as f:
        config = yaml.safe_load(f)

    # If user has a local config, override defaults
    if os.path.exists("../config.local.yaml"):
        with open("../config.local.yaml") as f:
            local_config = yaml.safe_load(f)
        config.update(local_config)

    return config

CONFIG = load_config()
BASE_PATH = CONFIG["base_path"]

def produce_index(vcf_gz):
    mapping = {}
    line_count = 0
    with gzip.open(vcf_gz, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            pos = int(fields[1])
            ref_allele = fields[3]
            alt_alleles = fields[4]
            #assert (pos, ref_allele, alt_alleles) not in mapping, f"Duplicate entry found for position {pos}, ref {ref_allele}, alt {alt_alleles}"
            mapping[(pos, ref_allele, alt_alleles)] = line_count
            line_count += 1
    print(f"Produced index for {vcf_gz} with {line_count} lines.")
    return mapping         

def get_mappings(chromosome):
    pangenome_path = f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome.vcf.gz"
    thousand_g_phased = f"{BASE_PATH}/starting_data/chr{chromosome}/1000g_phased.vcf.gz"
    thousand_g_alignments = f"{BASE_PATH}/starting_data/chr{chromosome}/PG.vcf.gz"

    pangenome_positions = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_positions.npy")
    thousand_g_phased_positions = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/1000g_phased_positions.npy")
    thousand_g_alignments_positions = np.load(f"{BASE_PATH}/starting_data/chr{chromosome}/1000g_positions.npy")

    pangenome_mapping = produce_index(pangenome_path)
    thousand_g_phased_mapping = produce_index(thousand_g_phased)
    thousand_g_alignments_mapping = produce_index(thousand_g_alignments)

    #assert len(pangenome_mapping) == len(pangenome_positions)
    #assert len(thousand_g_phased_mapping) == len(thousand_g_phased_positions)
    #assert len(thousand_g_alignments_mapping) == len(thousand_g_alignments_positions)

    pangenome_to_thousand_g_phased = {}
    pangenome_to_alignments = {}

    for key in pangenome_mapping:
        if key in thousand_g_phased_mapping:
            pangenome_to_thousand_g_phased[pangenome_mapping[key]] = thousand_g_phased_mapping[key]
        if key in thousand_g_alignments_mapping:
            pangenome_to_alignments[pangenome_mapping[key]] = thousand_g_alignments_mapping[key] 

    pangenome_mask_for_pangeonome_to_thousand_g_phased = []
    thousand_g_mask_for_pangenome_to_thousand_g_phased = []

    for key in pangenome_to_thousand_g_phased:
        pangenome_mask_for_pangeonome_to_thousand_g_phased.append(key)
        thousand_g_mask_for_pangenome_to_thousand_g_phased.append(pangenome_to_thousand_g_phased[key])
    pangenome_mask_for_pangeonome_to_thousand_g_phased = np.array(pangenome_mask_for_pangeonome_to_thousand_g_phased)
    thousand_g_mask_for_pangenome_to_thousand_g_phased = np.array(thousand_g_mask_for_pangenome_to_thousand_g_phased)

    np.save(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_mask.npy", pangenome_mask_for_pangeonome_to_thousand_g_phased)
    np.save(f"{BASE_PATH}/starting_data/chr{chromosome}/thousand_g_phased_mask.npy", thousand_g_mask_for_pangenome_to_thousand_g_phased)

    with open(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_to_thousand_g_phased.pickle", "wb") as f:
        pickle.dump(pangenome_to_thousand_g_phased, f)
    with open(f"{BASE_PATH}/starting_data/chr{chromosome}/pangenome_to_thousand_g_alignments.pickle", "wb") as f:
        pickle.dump(pangenome_to_alignments, f)

    print(f"Mappings saved for chromosome {chromosome}")
    print(f"pangenome_to_thousand_g_phased: {len(pangenome_to_thousand_g_phased)}")
    print(f"pangenome_to_thousand_g_alignments: {len(pangenome_to_alignments)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_mappings.py <chromosome>")
        sys.exit(1)

    chromosome = sys.argv[1]
    get_mappings(chromosome)