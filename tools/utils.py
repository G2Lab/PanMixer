import pandas as pd
import os
import gzip
import numpy as np
import time
import functools
import tracemalloc
import json

from constants import EXPERIMENT_PATH, DEMOGRAPHICS_CSV, AGGREGATION_DICTIONARY, TOTAL_UTILITY_JSON, ONEK_PHASED_SUBJECTS_NPY, STARTING_DATA_PATH

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

def aggregate_MIA_privacy(dfs, experiment_number):
    """Aggregate MIA privacy results across chromosomes."""
    n = len(dfs[1])

    # Check if any MIA results exist
    if not os.path.exists(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr1/0/MIA_privacy.npz"):
        return None

    mia_hap_rank = np.full(n, np.nan)
    mia_hap_frac_orig = np.full(n, np.nan)
    mia_hap_frac_best_pg = np.full(n, np.nan)
    mia_geno_rank = np.full(n, np.nan)
    mia_geno_frac_orig = np.full(n, np.nan)
    mia_geno_frac_best_pg = np.full(n, np.nan)
    mia_hap_frac_best_1000g = np.full(n, np.nan)

    for i in range(n):
        hap_match_orig_total = 0
        hap_total_orig_total = 0
        hap_match_pg_total = None
        hap_total_pg_total = None
        geno_match_orig_total = 0
        geno_total_orig_total = 0
        geno_match_pg_total = None
        geno_total_pg_total = None
        hap_match_1000g_total = None
        hap_total_1000g_total = None
        valid = True

        for chrom in range(1, 23):
            path = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chrom}/{i}/MIA_privacy.npz"
            if not os.path.exists(path):
                valid = False
                break
            d = np.load(path, allow_pickle=True)

            hap_match_orig_total += d["hap_match_orig"][0]
            hap_total_orig_total += d["hap_total_orig"][0]
            geno_match_orig_total += d["geno_match_orig"][0]
            geno_total_orig_total += d["geno_total_orig"][0]

            if hap_match_pg_total is None:
                hap_match_pg_total = d["hap_match_pg"].astype(np.int64)
                hap_total_pg_total = d["hap_total_pg"].astype(np.int64)
                geno_match_pg_total = d["geno_match_pg"].astype(np.int64)
                geno_total_pg_total = d["geno_total_pg"].astype(np.int64)
            else:
                hap_match_pg_total += d["hap_match_pg"]
                hap_total_pg_total += d["hap_total_pg"]
                geno_match_pg_total += d["geno_match_pg"]
                geno_total_pg_total += d["geno_total_pg"]

            if d["hap_match_1000g"].size > 0:
                if hap_match_1000g_total is None:
                    hap_match_1000g_total = d["hap_match_1000g"].astype(np.int64)
                    hap_total_1000g_total = d["hap_total_1000g"].astype(np.int64)
                else:
                    hap_match_1000g_total += d["hap_match_1000g"]
                    hap_total_1000g_total += d["hap_total_1000g"]

        if not valid:
            continue

        hap_frac_orig = hap_match_orig_total / hap_total_orig_total
        hap_frac_pg = hap_match_pg_total / hap_total_pg_total
        hap_rank = int(np.sum(hap_frac_pg >= hap_frac_orig)) + 1

        geno_frac_orig = geno_match_orig_total / geno_total_orig_total
        geno_frac_pg = geno_match_pg_total / geno_total_pg_total
        geno_rank = int(np.sum(geno_frac_pg >= geno_frac_orig)) + 1

        mia_hap_rank[i] = hap_rank
        mia_hap_frac_orig[i] = hap_frac_orig
        mia_hap_frac_best_pg[i] = np.max(hap_frac_pg)
        mia_geno_rank[i] = geno_rank
        mia_geno_frac_orig[i] = geno_frac_orig
        mia_geno_frac_best_pg[i] = np.max(geno_frac_pg)

        if hap_match_1000g_total is not None:
            hap_frac_1000g = hap_match_1000g_total / hap_total_1000g_total
            mia_hap_frac_best_1000g[i] = np.max(hap_frac_1000g)

    return (mia_hap_rank, mia_hap_frac_orig, mia_hap_frac_best_pg,
            mia_geno_rank, mia_geno_frac_orig, mia_geno_frac_best_pg,
            mia_hap_frac_best_1000g)


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

    try:
        gap_score_column = aggregate_gap_scores(dfs, experiment_number)
        if gap_score_column is not None:
            aggregated_df["gap_score"] = gap_score_column
        else:
            print("Gap score column is None")
    except Exception as e:
        print(f"[aggregate] gap_scores: ERROR, skipping ({e})")

    if "subject" in dfs[1].columns:
        try:
            genotypes_scores, haplotype_scores, haplotype_both_scores, no_weight_scores = aggregate_gap_scores_all(dfs, experiment_number)

            if genotypes_scores is not None:
                aggregated_df["genotypes_score"] = genotypes_scores
            else:
                print("No genotypes scores returned.")

            if haplotype_scores is not None:
                aggregated_df["haplotype_score"] = haplotype_scores
            else:
                print("No haplotype scores returned.")

            if haplotype_both_scores is not None:
                aggregated_df["haplotype_both_score"] = haplotype_both_scores
            else:
                print("No haplotype_both scores returned.")

            if no_weight_scores is not None:
                aggregated_df["no_weight_score"] = no_weight_scores
            else:
                print("No no_weight scores returned.")
        except Exception as e:
            print(f"[aggregate] gap_scores_all: ERROR, skipping ({e})")

        try:
            reverse_genotypes, reverse_haplotypes, reverse_haplotypes_both, reverse_no_weight = aggregate_reverse_gap_scores_all(dfs, experiment_number)

            if reverse_genotypes is not None:
                aggregated_df["reverse_genotypes_score"] = reverse_genotypes
            if reverse_haplotypes is not None:
                aggregated_df["reverse_haplotype_score"] = reverse_haplotypes
            if reverse_haplotypes_both is not None:
                aggregated_df["reverse_haplotype_both_score"] = reverse_haplotypes_both
            if reverse_no_weight is not None:
                aggregated_df["reverse_no_weight_score"] = reverse_no_weight
        except Exception as e:
            print(f"[aggregate] reverse_gap_scores: ERROR, skipping ({e})")

        try:
            add_true_utility_loss(aggregated_df, experiment_number)
        except Exception as e:
            print(f"[aggregate] true_utility_loss: ERROR, skipping ({e})")

    # MIA privacy aggregation
    try:
        mia_results = aggregate_MIA_privacy(dfs, experiment_number)
        if mia_results is not None:
            hap_rank, hap_frac_orig, hap_frac_best_pg, geno_rank, geno_frac_orig, geno_frac_best_pg, hap_frac_best_1000g = mia_results
            aggregated_df["mia_hap_rank"] = hap_rank
            aggregated_df["mia_hap_frac_vs_original"] = hap_frac_orig
            aggregated_df["mia_hap_frac_vs_best_pg"] = hap_frac_best_pg
            aggregated_df["mia_geno_rank"] = geno_rank
            aggregated_df["mia_geno_frac_vs_original"] = geno_frac_orig
            aggregated_df["mia_geno_frac_vs_best_pg"] = geno_frac_best_pg
            aggregated_df["mia_hap_frac_vs_best_1000g"] = hap_frac_best_1000g
        else:
            print("No MIA privacy results found.")
    except Exception as e:
        print(f"[aggregate] MIA_privacy: ERROR, skipping ({e})")

    return aggregated_df

def aggregate_score_metric(experiment_number, metric, n, obfuscated_subjects):
    scores = np.zeros(n)
    subjects_30x = list(np.load(STARTING_DATA_PATH + "/chr1/1000g_30x_phased_subjects.npy"))
    subjects_phased = list(np.load(STARTING_DATA_PATH + "/chr1/" + ONEK_PHASED_SUBJECTS_NPY))
    print("Aggregating metric", metric)

    lowest_score_seen_subject = {}
    for i in range(n):
        if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}.npy"):
            scores[i] = np.nan
            continue
        scores_all = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}.npy")
        g_to_gstar = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{metric}_self.npy")
        for chrom in range(2,23):
            if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{metric}.npy"):
                scores[i] = np.nan
                continue
            scores_all += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{metric}.npy")
            g_to_gstar += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{metric}_self.npy")

        if scores_all.shape[0] == len(subjects_30x):
            thousand_g_subjects = subjects_30x
        elif scores_all.shape[0] == len(subjects_phased):
            thousand_g_subjects = subjects_phased
        else:
            raise AssertionError(f"Scores length {scores_all.shape[0]} matches neither 30x ({len(subjects_30x)}) nor phased ({len(subjects_phased)}) for metric {metric} at index {i}")

        if obfuscated_subjects[i] in thousand_g_subjects:
            subject_index = thousand_g_subjects.index(obfuscated_subjects[i])
            # ignore subject index in max calculation
            scores[i] = g_to_gstar[0] - np.max(np.delete(scores_all, subject_index))
        else:
            scores[i] = g_to_gstar[0] - np.max(scores_all)
        
        if obfuscated_subjects[i] not in lowest_score_seen_subject:
            lowest_score_seen_subject[obfuscated_subjects[i]] = scores[i]
        else:
            if scores[i] < lowest_score_seen_subject[obfuscated_subjects[i]]:
                lowest_score_seen_subject[obfuscated_subjects[i]] = scores[i]
    
    return scores

def aggregate_reverse_score_metric(experiment_number, metric, n, obfuscated_subjects):
    """Aggregate reverse gap scores across chromosomes.

    Reverse gap score: query=original, database=other pangenome members, self=obfuscated.
    Files: reverse_{metric}_scores_others.npy and reverse_{metric}_score_obfuscated.npy
    """
    scores = np.zeros(n)
    pangenome_subjects = list(np.load(STARTING_DATA_PATH + "/chr1/pangenome_subjects.npy"))
    print("Aggregating reverse metric", metric)

    for i in range(n):
        others_file = f"reverse_{metric}_scores_others.npy"
        obfuscated_file = f"reverse_{metric}_score_obfuscated.npy"

        if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{others_file}"):
            scores[i] = np.nan
            continue

        scores_others = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{others_file}")
        score_obfuscated = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/{i}/{obfuscated_file}")

        for chrom in range(2, 23):
            if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{others_file}"):
                scores[i] = np.nan
                continue
            scores_others += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{others_file}")
            score_obfuscated += np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chrom}/{i}/{obfuscated_file}")

        scores[i] = score_obfuscated[0] - np.max(scores_others)

    return scores


def aggregate_reverse_gap_scores_all(dfs, experiment_number):
    n = len(dfs[1])
    obfuscated_subjects = dfs[1]["subject"].values
    genotype_scores = aggregate_reverse_score_metric(experiment_number, "genotypes", n, obfuscated_subjects)
    haplotype_scores = aggregate_reverse_score_metric(experiment_number, "haplotypes", n, obfuscated_subjects)
    haplotype_both_scores = aggregate_reverse_score_metric(experiment_number, "haplotypes_both", n, obfuscated_subjects)
    no_weight_scores = aggregate_reverse_score_metric(experiment_number, "no_weight", n, obfuscated_subjects)
    return genotype_scores, haplotype_scores, haplotype_both_scores, no_weight_scores


def aggregate_gap_scores_all(dfs, experiment_number):
    n = len(dfs[1])
    obfuscated_subjects = dfs[1]["subject"].values
    genotype_scores = aggregate_score_metric(experiment_number, "genotypes_scores", n, obfuscated_subjects)
    haplotype_scores = aggregate_score_metric(experiment_number, "haplotypes_scores", n, obfuscated_subjects)
    haplotype_both_scores = aggregate_score_metric(experiment_number, "haplotypes_both_scores", n, obfuscated_subjects)
    no_weight_scores = aggregate_score_metric(experiment_number, "no_weight_scores", n, obfuscated_subjects)

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

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError(f"Invalid boolean value: {value}")
