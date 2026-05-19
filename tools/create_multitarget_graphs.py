import slurm_helper
from constants import STARTING_DATA_PATH, EXPERIMENT_PATH
import os
import sys
from tools.utils import load_demographic_data, add_demographic_information
from tools.experiment_starter import find_next_experiment_number
import pandas as pd
import numpy as np


def create_new_experiment_directory(new_experiment_number, target_individuals_list):
    #directories for data
    os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}", exist_ok=True)
    os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data", exist_ok=True)
    for chromosome in range(1, 23):
        os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data/chr{chromosome}", exist_ok=True)
    os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data/all", exist_ok=True)

    num_new_graphs = len(target_individuals_list)
    for i in range(num_new_graphs):
        for chromosome in range(1, 23):
            os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data/chr{chromosome}/{i}", exist_ok=True)
        os.makedirs(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data/all/{i}", exist_ok=True)

    graphs = []
    for i in range(num_new_graphs):
        if len(target_individuals_list[i]) == 1:
            graphs.append({
                "num_targets": 1,
                "target_individuals": str(target_individuals_list[i][0])
            })
        else:
            graphs.append({
                "num_targets": len(target_individuals_list[i]),
                "target_individuals": ",".join([str(x) for x in target_individuals_list[i]]),
            })
    graphs_df = pd.DataFrame(graphs)
    graphs_df.to_csv(EXPERIMENT_PATH + f"/exp_{new_experiment_number}/data.csv", index=False)

def create_random_target_lists(number_of_target_individuals_list):
    #sample without replacement from 0 to 43 inclusive 
    target_individuals_list = []
    for i in number_of_target_individuals_list:
        choices = np.random.choice(range(44), size=i, replace=False)
        choices = choices.tolist()
        choices.sort()
        target_individuals_list.append(choices)
    return target_individuals_list

def create_multitarget_graphs(experiment_number, target_individuals_list):
    next_experiment_number = find_next_experiment_number()
    create_new_experiment_directory(next_experiment_number, target_individuals_list)

def main():
    if len(sys.argv) not in (3, 4, 5):
        print("Usage: python tools/create_multitarget_graphs.py <experiment_number> <target_counts_csv> [--seed SEED]")
        sys.exit(1)

    experiment_number = int(sys.argv[1])
    number_of_target_individuals_list = str(sys.argv[2])
    if len(sys.argv) == 5 and sys.argv[3] == "--seed":
        seed = int(sys.argv[4])
    elif len(sys.argv) == 4:
        seed = None if sys.argv[3] == "" else int(sys.argv[3])
    elif len(sys.argv) == 3:
        seed = None
    else:
        print("Usage: python tools/create_multitarget_graphs.py <experiment_number> <target_counts_csv> [--seed SEED]")
        sys.exit(1)
    if seed is not None:
        np.random.seed(seed)
        print(f"Using random seed {seed}")
    number_of_target_individuals_list = number_of_target_individuals_list.split(",")
    number_of_target_individuals_list = [int(x) for x in number_of_target_individuals_list]

    target_individuals_list = create_random_target_lists(number_of_target_individuals_list)
    
    create_multitarget_graphs(experiment_number, target_individuals_list)

if __name__ == "__main__":
    main()
