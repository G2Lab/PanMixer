import json
import numpy as np

subjects = np.loadtxt("/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/subjects.txt", dtype=str)
STARTING_DATA_PATH = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data"


true_max_utility_loss = {}
for i,subject in enumerate(subjects):
    true_max_utility_loss[subject] = 0
    for chrom in range(1, 23):
        utility_loss = np.load(f"{STARTING_DATA_PATH}/chr{chrom}/subjects/{i}/utility_loss.npy")
        true_max_utility_loss[subject] += 2 * np.sum(utility_loss)

with open(f"{STARTING_DATA_PATH}/utility_loss.json", "w") as f:
    json.dump(true_max_utility_loss, f)