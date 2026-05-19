import numpy as np
import pickle
import sys
from ortools.linear_solver import pywraplp
import sys
from tools.slurm_helper import launch_job
from tools.utils import load_data, add_command, get_chromosome_path
import pandas
import json
import os
import gzip

from constants import (
    PANGENOME_HAPLOTYPES,   
)

def beagle_stats(experiment_number, starting_data_path, chromosome):
    data = load_data(experiment_number)

    num_tasks = len(data)

    chromosome_path = get_chromosome_path(starting_data_path, chromosome)

    args = [experiment_number, chromosome_path]

    launch_job("beagle_stats", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

def compute_wged(support_vector, original_haplotypes, inferred_haplotypes):
    diff = original_haplotypes != inferred_haplotypes
    return np.sum(support_vector * diff)

def convert_gz_vcf_to_npy(vcf_file_path):
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
    genotypes = np.zeros((num_subjects, num_sites, 2), dtype=np.int8)

    positions = []

    with open_func(vcf_file_path, 'rt') as vcf:
        site = 0
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            positions.append(int(fields[1]))
            for i, subject in enumerate(subjects_in_vcf):
                try:
                    genotype = fields[9+i].split(':')[0].split('|')
                except:
                    print(fields)
                    print(fields[9+i])
                    print(fields[9+i].split(':'))
                    print(fields[9+i].split(':')[0])
                    print(fields[9+i].split(':')[0].split('|'))
                    raise
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

def compute_stats(experiment_number, i, beagle_vcf_as_npy, beagle_subjects, chromosome_path):
    df = load_data(experiment_number)
    pangenome_haplotypes = np.load(PANGENOME_HAPLOTYPES)
    pangenome_subjects = np.load(PANGENOME_SUBJECTS)

    pangenome_support = np.load(chromosome_path + "pangenome_support.npy")
    pangenome_support = np.stack((pangenome_support, pangenome_support), axis=1)

    subject_name = df["subject"][i]
    subject_index = np.where(pangenome_subjects == subject_name)[0][0]

    original_haplotypes = pangenome_haplotypes[subject_index]

    subject_index_beagle = list(beagle_subjects).index(subject_name)
    inferred_haplotypes = beagle_vcf_as_npy[subject_index_beagle]

    wged = compute_wged(pangenome_support, original_haplotypes, inferred_haplotypes)

    blocks_dict = json.load(open(chromosome_path + "pop_blocks_dict.json"))["all"]
    precomputed_pmi_vector = np.load(chromosome_path + f"pmi/{subject_name}_all.npy")

    pmi = 0
    for j,block in enumerate(blocks_dict):
        indicies = blocks_dict[block][1]
        
        original_haplotypes_block = original_haplotypes[indicies]
        inferred_haplotypes_block = inferred_haplotypes[indicies]

        if np.all(original_haplotypes_block[:, 0] == inferred_haplotypes_block[:, 0]):
            pmi += precomputed_pmi_vector[j, 0]
        if np.all(original_haplotypes_block[:, 1] == inferred_haplotypes_block[:, 1]):
            pmi += precomputed_pmi_vector[j, 1]

    max_pmi = np.sum(precomputed_pmi_vector)
    max_wged = np.sum(pangenome_support)
    results = {
        "wged": wged,
        "wged_normalized": wged / max_wged,
        "pmi": pmi,
        "pmi_normalized": pmi / max_pmi,
    }
    print(results)
    with open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/{i}/beagle_stats.json", "w") as f:
        json.dump(results, f, indent=4)

def main():
    experiment_number = sys.argv[1]
    chromosome_path = sys.argv[2]
    row_id = int(sys.argv[3])

    vcf_path = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/{row_id}/new_haplotypes_beagle.vcf.gz.vcf.gz"
    
    beagle_vcf_as_npy, beagle_subjects, positions = convert_gz_vcf_to_npy(vcf_path)

    compute_stats(experiment_number, row_id, beagle_vcf_as_npy, beagle_subjects, chromosome_path)
    
if __name__ == '__main__':
    main()