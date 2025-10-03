import numpy as np
import pandas as pd
import os
from tools.utils import load_demographic_data, add_demographic_information

from constants import EXPERIMENT_PATH

NUM_CHROMOSOMES = 22

def find_next_experiment_number():
    experiment_number = 0
    while True:
        experiment_path = EXPERIMENT_PATH + f"/exp_{experiment_number}/"
        try:
            os.mkdir(experiment_path)
            return experiment_number
        except FileExistsError:
            experiment_number += 1

def experiment_starter(capacities_file, subjects_file, baseline_unedited, baseline_empty):
    capacities = np.loadtxt(capacities_file)
    subjects = np.loadtxt(subjects_file, dtype=str)

    if len(subjects.shape) == 0:
        subjects = np.array([subjects])

    if len(capacities.shape) == 0:
        capacities = np.array([capacities])

    capacities_float = capacities.astype(float)

    exp_number = find_next_experiment_number()
    print(f"Starting experiment {exp_number}")
    print(f"Running with just the unedited pangenome: {baseline_unedited}")
    print(f"Running with just the empty pangenome: {baseline_empty}")

    assert not (baseline_unedited and baseline_empty)

    os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/data")

    data_points = []
    if baseline_unedited:
        for subject in subjects:
            data_points.append({"subject": subject, "capacity": 0.0})
    elif baseline_empty:
        for subject in subjects:
            data_points.append({"subject": subject, "capacity": 1.0})
    else:
        for subject in subjects:
            for capacity in capacities_float:
                data_points.append({"subject": subject, "capacity": capacity})

    data_points_df = pd.DataFrame(data_points)

    demographic_data = load_demographic_data()
    data_points_df = add_demographic_information(data_points_df, demographic_data)


    for chromosome in range(1, NUM_CHROMOSOMES + 1):
        os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/data/chr{chromosome}")
        for i in range(len(data_points_df)):
            os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/data/chr{chromosome}/{i}")

    os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/data/all")
    for i in range(len(data_points_df)):
        os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/data/all/{i}")

    os.mkdir(EXPERIMENT_PATH + f"/exp_{exp_number}/plots")

    data_points_df.to_csv(EXPERIMENT_PATH + f"/exp_{exp_number}/data.csv", index=False)
    with open(EXPERIMENT_PATH + f"/exp_{exp_number}/config.txt", "w") as file:
        file.write(f"Experiment {exp_number} started with capacities_file: {capacities_file} and subjects_file: {subjects_file}\n")
    with open(EXPERIMENT_PATH + f"/exp_{exp_number}/commands.csv", "w") as file:
        file.write(f"command,slurm_path\n")