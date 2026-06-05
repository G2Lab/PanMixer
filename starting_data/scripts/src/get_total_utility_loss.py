import json
import numpy as np
import time

from paths import BASE_PATH

subjects = np.loadtxt(f"{BASE_PATH}/starting_data/subjects.txt", dtype=str)
STARTING_DATA_PATH = f"{BASE_PATH}/starting_data"


t_total = time.time()
true_max_utility_loss = {}
for i,subject in enumerate(subjects):
    t_subject = time.time()
    true_max_utility_loss[subject] = 0
    for chrom in range(1, 23):
        utility_loss = np.load(f"{STARTING_DATA_PATH}/chr{chrom}/subjects/{i}/utility_loss.npy")
        true_max_utility_loss[subject] += 2 * np.sum(utility_loss)
    print(f"[TIMING] Subject {subject} ({i+1}/{len(subjects)}) completed in {time.time() - t_subject:.2f}s")

with open(f"{STARTING_DATA_PATH}/utility_loss.json", "w") as f:
    json.dump(true_max_utility_loss, f)

print(f"[TIMING] Total utility loss computation completed in {time.time() - t_total:.2f}s")
