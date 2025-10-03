import pandas as pd
import os
import gzip
import numpy as np
import time
import functools
import tracemalloc
import json

from constants import EXPERIMENT_PATH, DEMOGRAPHICS_CSV, AGGREGATION_DICTIONARY, TOTAL_UTILITY_JSON

LIMIT = 1000

def load_demographic_data():
    return pd.read_csv(DEMOGRAPHICS_CSV)

def add_demographic_information(data_df, demographic_data):
    data_df = data_df.merge(demographic_data, on="subject")
    return data_df

def load_data(experiment_number, chromosome = 1):
    if os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_chr{chromosome}.csv"):
        data = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_chr{chromosome}.csv")
        data["subject"] = data["subject"].astype(str)
        return data
    
    data = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data.csv")
    data["subject"] = data["subject"].astype(str)
    return data

def store_data(data, experiment_number):
    data.to_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data.csv", index=False)

def store_data_multichromosome(data, experiment_number):
    for chromosome in range(1, 23):
        data_chr = data[chromosome]
        data_chr.to_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_chr{chromosome}.csv", index=False)
    
    aggregated_df = produced_concat_df(data, experiment_number)
    aggregated_df.to_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_all.csv", index=False)

def load_data_multichromosome(experiment_number):
    data = {}
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_chr1.csv"):
        normal_data = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data.csv")
        for chromosome in range(1, 23):
            data[chromosome] = normal_data.copy()
        return data
    for chromosome in range(1, 23):
        data_chr = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_chr{chromosome}.csv")
        data[chromosome] = data_chr
    return data

def load_data_all(experiment_number):
    if os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_all.csv"):
        data = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data_all.csv")
        data["subject"] = data["subject"].astype(str)
        return data
    else:
        raise FileNotFoundError("Data not found for experiment number: " + str(experiment_number))

def add_command(experiment_number, command, slurm_path):
    with open(EXPERIMENT_PATH + f"/exp_{experiment_number}/commands.csv", "a") as file:
        file.write(f"{command},{slurm_path}\n")

def get_chromosome_path(starting_data_path, chromosome):
    if chromosome == -1:
        path_with_chromosome = starting_data_path + "/" + "full/"
    else:
        path_with_chromosome = starting_data_path + "/" + "chr" + str(chromosome) + "/"
    return path_with_chromosome

def latest_experiment_number():
    experiment_number = LIMIT
    while experiment_number >= 0:
        experiment_path = EXPERIMENT_PATH + f"/exp_{experiment_number}/"
        
        if os.path.exists(experiment_path):
            return experiment_number
        else:
            experiment_number -= 1
        
    print("Run experiment starter to start running an experiment")
    return -1

def get_numpy_matrices(vcf_file_path):
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
                    try:
                        genotype = fields[9+i].split(':')[0].split('/')
                    except:
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

def produced_concat_df(dfs, experiment_number):
    assert len(dfs) == 22, "Number of chromosomes must be 22"

    columns = dfs[1].columns
    columns_intersection = [col for col in columns if col in AGGREGATION_DICTIONARY.keys()]
    aggregated_df = pd.DataFrame(columns=columns_intersection)

    for column in columns:
        if column in AGGREGATION_DICTIONARY:
            if AGGREGATION_DICTIONARY[column] == "index":
                aggregated_df[column] = dfs[1][column].values
            elif AGGREGATION_DICTIONARY[column] == "mean":
                all_values = []
                for chrom,df in dfs.items():
                    all_values.append(df[column].values)
                all_values = np.array(all_values)
                all_values = np.mean(all_values, axis=0)
                aggregated_df[column] = all_values
            elif AGGREGATION_DICTIONARY[column] == "sum":
                all_values = []
                for chrom,df in dfs.items():
                    all_values.append(df[column].values)
                all_values = np.array(all_values)
                all_values = np.sum(all_values, axis=0)
                aggregated_df[column] = all_values
            else:
                raise ValueError(f"Unknown aggregation method: {AGGREGATION_DICTIONARY[column]}")

    gap_score_column = aggregate_gap_scores(dfs, experiment_number)
    if gap_score_column is not None:
        aggregated_df["gap_score"] = gap_score_column
    else:
        print("Gap score column is None")

    genotypes_scores, haplotype_scores, haplotype_both_scores, no_weight_scores = aggregate_gap_scores_all(dfs, experiment_number)

    if genotypes_scores is not None:
        aggregated_df["genotypes_score"] = genotypes_scores
    else:
        print("⚠️ No genotypes scores returned.")

    if haplotype_scores is not None:
        aggregated_df["haplotype_score"] = haplotype_scores
    else:
        print("⚠️ No haplotype scores returned.")

    if haplotype_both_scores is not None:
        aggregated_df["haplotype_both_score"] = haplotype_both_scores
    else:
        print("⚠️ No haplotype_both scores returned.")

    if no_weight_scores is not None:
        aggregated_df["no_weight_score"] = no_weight_scores
    else:
        print("⚠️ No no_weight scores returned.")

    add_true_utility_loss(aggregated_df, experiment_number)
    return aggregated_df

def aggregate_score_metric(experiment_number, metric, n):
    scores = np.zeros(n)
    print("Aggregating metric", metric)
    for i in range(n):
        if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}.npy"):
            return None        
        scores_all = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}.npy")
        g_to_gstar = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}_self.npy")
        for chrom in range(2,23):
            scores_all += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{metric}.npy")
            g_to_gstar += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{metric}_self.npy")
        
        scores[i] = g_to_gstar[0] - np.max(scores_all)
    
    return scores

def aggregate_gap_scores_all(dfs, experiment_number):
    n = len(dfs[1])
    genotype_scores = aggregate_score_metric(experiment_number, "genotypes_scores", n)
    haplotype_scores = aggregate_score_metric(experiment_number, "haplotypes_scores", n)
    haplotype_both_scores = aggregate_score_metric(experiment_number, "haplotypes_both_scores", n)
    no_weight_scores = aggregate_score_metric(experiment_number, "no_weight_scores", n)

    return genotype_scores, haplotype_scores, haplotype_both_scores, no_weight_scores


def aggregate_gap_scores(dfs, experiment_number):
    gap_score_all_chrom = np.zeros(len(dfs[1]))

    for i in range(len(dfs[1])):
        if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/gap_score.npy"):
            return None
        scores_all = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/gap_score.npy")
        g_to_gstar = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/gap_score_against_self.npy")
        for chrom in range(2,23):
            scores_all += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/gap_score.npy")
            g_to_gstar += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/gap_score_against_self.npy")
        
        gap_score_all_chrom[i] = g_to_gstar[0] - np.max(scores_all)
    
    return gap_score_all_chrom

def add_true_utility_loss(df_all, experiment_number):
    true_max_utility_loss = json.load(open(f"{TOTAL_UTILITY_JSON}"))
    
    true_max_utility_loss_column = []
    for i in range(len(df_all)):
        true_max_utility_loss_column.append(true_max_utility_loss[df_all["subject"].iloc[i]])
    
    df_all["true_max_utility_loss"] = true_max_utility_loss_column

def profile(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        tracemalloc.start()
        start_time = time.perf_counter()
        
        result = func(*args, **kwargs)
        
        end_time = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        duration = end_time - start_time
        print(f"\n--- Profiling: {func.__name__} ---")
        print(f"Time elapsed: {duration:.6f} seconds")
        print(f"Current memory usage: {current / 10**6:.6f} MB")
        print(f"Peak memory usage: {peak / 10**6:.6f} MB")
        print(f"-----------------------------\n")
        
        return result
    return wrapper

def profile_timings(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        
        result = func(*args, **kwargs)
        
        end_time = time.perf_counter()
        duration = end_time - start_time
        print(f"\n--- Timings: {func.__name__} ---")
        print(f"Time elapsed: {duration:.6f} seconds")
        print(f"-----------------------------\n")
        
        return result
    return wrapper
