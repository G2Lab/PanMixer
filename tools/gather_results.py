import numpy as np
import pickle
import json
import pandas as pd
import os

from tools.utils import load_data_multichromosome, store_data_multichromosome, profile_timings
from tools.giraffe_parser import read_stats_file

from constants import SIZE_OF_SVS, TYPES_OF_SVS, TYPES_OF_RECONSTRUCTION, EXPERIMENT_PATH, READ_POPULATIONS_FILE, READ_SUBJECTS_FILE
from constants import PANGENIE_COLUMNS, AF_PATH, PANGENOME_HAPLOTYPES, PANGENOME_SUBJECTS, PANGENOME_SUPPORT
from constants import ONEK_NPY, POPULATION_AF_FREQS, NUM_1000G_SUBJECTS
from constants import STATS_FILE_SCHEMA
from constants import VERBOSE, MAX_LD_DISTANCE_KBP

def add_ld_results_by_chromosome(df, experiment_number, chromosome):
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr1/0/ld_loss.txt"):
        print("No LD results found")
        return

    ld_results = []
    for i in range(len(df)):
        result = float(open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/ld_loss.txt", "r").read())
        ld_results.append(result)

    df["ld"] = ld_results
    return df

@profile_timings
def add_ld_results(dfs, experiment_number, overwrite=False):
    if not overwrite and "ld" in dfs[1].columns:
        return dfs
    
    for i in range(1, 23):
        dfs[i] = add_ld_results_by_chromosome(dfs[i], experiment_number, i)
    return dfs

def add_af_results_by_chromosome(df, experiment_number, chromosome):
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/0/af_loss.json"):
        print("No af loss results found")
        return
    
    af_loss_all = []
    af_loss_all_snp_only = []
    af_complex = []
    af_loss_snps_0_01 = []
    af_loss_snps_0_05 = []
    af_loss_snps_0_1 = []
    af_loss_snps_0_5 = []

    for i in range(len(df)):
        results = json.load(open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/af_loss.json"))

        af_loss_all.append(results["all"])
        af_loss_all_snp_only.append(results["snp_only"])
        af_complex.append(results["complex"])

        af_loss_snps_0_01.append(results["snp_only_maf_0.01"])
        af_loss_snps_0_05.append(results["snp_only_maf_0.05"])
        af_loss_snps_0_1.append(results["snp_only_maf_0.1"])
        af_loss_snps_0_5.append(results["snp_only_maf_0.5"])
    
    df["af_loss"] = af_loss_all
    df["af_loss_snp_only"] = af_loss_all_snp_only
    df["af_complex"] = af_complex
    df["af_loss_snps_0_01"] = af_loss_snps_0_01
    df["af_loss_snps_0_05"] = af_loss_snps_0_05
    df["af_loss_snps_0_1"] = af_loss_snps_0_1
    df["af_loss_snps_0_5"] = af_loss_snps_0_5

    df["wd_af"] = 2 * np.array(df["af_loss"])
    df["wd_af_snp_only"] = 2 * np.array(df["af_loss_snp_only"])
    df["wd_af_complex"] = 2 * np.array(df["af_complex"])
    df["wd_af_snps_0_01"] = 2 * np.array(df["af_loss_snps_0_01"])
    df["wd_af_snps_0_05"] = 2 * np.array(df["af_loss_snps_0_05"])
    df["wd_af_snps_0_1"] = 2 * np.array(df["af_loss_snps_0_1"])
    df["wd_af_snps_0_5"] = 2 * np.array(df["af_loss_snps_0_5"])

    return df

@profile_timings
def add_af_results(dfs, experiment_number, overwrite=False):
    if not overwrite and "af_loss" in dfs[1].columns:
        return dfs
    
    for i in range(1, 23):
        dfs[i] = add_af_results_by_chromosome(dfs[i], experiment_number, i)
    return dfs

def add_maf_ld_results_chromosome(df, experiment_number, chromosome):
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/0/results-maf.csv"):
        print("No maf ld results found")
        return df
    
    maf_kl = []
    maf_wd = []

    for i in range(len(df)):
        maf_kl.append(pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/results-maf.csv")["kl"].values[0])
        maf_wd.append(pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/results-maf.csv")["wd"].values[0])

    ld_euclidean = []
    
    for i in range(len(df)):
        ld_csv = pd.read_csv(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/results-ld-decay.csv")
        real_ld = ld_csv["real_y"].values
        syn_ld = ld_csv["syn_y"].values

        #print(real_ld - syn_ld)
        #print(np.linalg.norm(real_ld - syn_ld))

        real_ld = real_ld[:int(MAX_LD_DISTANCE_KBP * 1000)]
        syn_ld = syn_ld[:int(MAX_LD_DISTANCE_KBP * 1000)]

        ld_euclidean.append(np.linalg.norm(real_ld - syn_ld, ord=1))

    df["maf_kl"] = maf_kl
    df["maf_wd"] = maf_wd

    df["ld_euclidean"] = ld_euclidean

    return df

@profile_timings
def add_maf_ld_results(dfs, experiment_number, overwrite=False):
    if not overwrite and "maf_kl" in dfs[1].columns:
        return dfs
    
    for i in range(1, 23):
        dfs[i] = add_maf_ld_results_chromosome(dfs[i], experiment_number, i)
    return dfs

def add_stacker_results_chromosome(df, experiment_number, chromosome):
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/0/stacked.json"):
        print("No stacker results found")
        return df

    stacker_results = []
    for i in range(len(df)):
        stacker_results.append(json.load(open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/stacked.json")))

    # Remove old columns if they exist in the DataFrame
    for col in stacker_results[0].keys():
        if col in df.columns:
            df = df.drop(columns=[col])

    # Convert stacker_results to a DataFrame and concatenate
    stacker_results = pd.DataFrame(stacker_results)
    df = pd.concat([df, stacker_results], axis=1)

    return df

@profile_timings
def add_stacker_results(dfs, experiment_number, overwrite=False):
    if not overwrite and "stacker_strategy" in dfs[1].columns:
        return dfs
    
    for i in range(1, 23):
        dfs[i] = add_stacker_results_chromosome(dfs[i], experiment_number, i)
    return dfs

def add_optimizer_results_chromosome(df, experiment_number, chromosome):
    if not os.path.exists(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/0/optimizer_stats.json"):
        print("No optimizer results found")
        return df

    optimizer_results = []
    for i in range(len(df)):
        optimizer_results.append(json.load(open(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/optimizer_stats.json")))
        # Remove old columns if they exist in the DataFrame
    for col in optimizer_results[0].keys():
        if col in df.columns:
            df = df.drop(columns=[col])

    optimizer_results = pd.DataFrame(optimizer_results)
    df = pd.concat([df, optimizer_results], axis=1)

    return df

@profile_timings
def add_optimizer_results(dfs, experiment_number, overwrite=False):
    if not overwrite and "optimizer" in dfs[1].columns:
        return dfs
    
    for i in range(1, 23):
        dfs[i] = add_optimizer_results_chromosome(dfs[i], experiment_number, i)
    return dfs


def add_giraffe_results_by_chromosome(df, experiment_number, chromosome):
    read_subjects = {
        "HG00138": "EUR",
        "HG00635": "EAS",
        "HG01112": "AMR",
        "HG01600": "EAS",
        "HG02698": "SAS",
        "NA12778": "EUR",
        "NA18853": "AFR",
    }
    giraffe_results = []

    for i in range(len(df)):
        experiment_row_results = {}
        for read_subject in read_subjects:
            population = read_subjects[read_subject]

            stats_filepath = EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/{read_subject}.giraffe.stats.txt"

            try:
                results = read_stats_file(stats_filepath)
            except:
                results = {}

            for k,v in results.items():
                experiment_row_results[f"{read_subject}_{k}"] = v
        giraffe_results.append(experiment_row_results)
    
    for col in giraffe_results[0].keys():
        for existing_col in df.columns:
            if col in existing_col:
                df = df.drop(columns=[existing_col])

    df_giraffe = pd.DataFrame(giraffe_results)
    df = pd.concat([df, df_giraffe], axis=1)
    return df

def add_giraffe_results(dfs, experiment_number, overwrite = False):
    if not overwrite and "giraffe" in dfs[1].columns:
        return dfs
    dfs[21] = add_giraffe_results_by_chromosome(dfs[21], experiment_number, 21)
    return dfs

def gather_results(experiment_number, overwrite, get_just_optimizer, get_index, get_just_gap_score, get_just_stacker, get_just_af_loss, get_just_pangenie_stats, get_just_ld_loss, get_just_accuracy_stats, get_just_maf_ld, get_just_giraffe):
    data = load_data_multichromosome(experiment_number)

    if get_just_optimizer:
        data = add_optimizer_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return

    if get_just_stacker:
        data = add_stacker_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return
    
    if get_just_maf_ld:
        data = add_maf_ld_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return

    if get_just_af_loss:
        data = add_af_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return
    
    if get_just_ld_loss:
        add_ld_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return
    
    if get_just_giraffe:
        add_giraffe_results(data, experiment_number, overwrite)
        store_data_multichromosome(data, experiment_number)
        return   
    """
    if get_index:
        add_index(data, experiment_number, overwrite)
        store_data(data, experiment_number)
        return
    

    
    if get_just_stacker:
        data = add_stacker_results(data, experiment_number, overwrite)
        print(data)
        store_data(data, experiment_number)
        return

    if get_just_af_loss:
        add_af_results(data, experiment_number, overwrite)
        store_data(data, experiment_number)
        return
    
    if get_just_pangenie_stats:
        data = add_pangenie_results(data, experiment_number, overwrite)
        store_data(data, experiment_number)
        return
    
    if get_just_ld_loss:
        add_ld_results(data, experiment_number, overwrite)
        store_data(data, experiment_number)
        return
    
    if get_just_accuracy_stats:
        add_accuracy_stats(data, experiment_number, overwrite)
        store_data(data, experiment_number)
    
    if get_adjusted_wged:
        add_adjusted_wged(data, experiment_number, overwrite)
        store_data(data, experiment_number)
        return
    """
    data = add_optimizer_results(data, experiment_number, overwrite)
    data = add_stacker_results(data, experiment_number, overwrite)
    data = add_maf_ld_results(data, experiment_number, overwrite)
    data = add_af_results(data, experiment_number, overwrite)
    data = add_ld_results(data, experiment_number, overwrite)
    data = add_giraffe_results(data, experiment_number, overwrite)

    #add_gap_score_results(data, experiment_number, overwrite)
    #data = add_stacker_results(data, experiment_number, overwrite)
    #add_af_results(data, experiment_number, overwrite)
    #data = add_pangenie_results(data, experiment_number, overwrite)
    #add_ld_results(data, experiment_number, overwrite)
    #add_accuracy_stats(data, experiment_number, overwrite)
    #add_adjusted_wged(data, experiment_number, overwrite)

    store_data_multichromosome(data, experiment_number)