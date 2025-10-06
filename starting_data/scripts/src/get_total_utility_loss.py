import json
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

subjects = np.loadtxt(f"{BASE_PATH}/starting_data/subjects.txt", dtype=str)
STARTING_DATA_PATH = f"{BASE_PATH}/starting_data"


true_max_utility_loss = {}
for i,subject in enumerate(subjects):
    true_max_utility_loss[subject] = 0
    for chrom in range(1, 23):
        utility_loss = np.load(f"{STARTING_DATA_PATH}/chr{chrom}/subjects/{i}/utility_loss.npy")
        true_max_utility_loss[subject] += 2 * np.sum(utility_loss)

with open(f"{STARTING_DATA_PATH}/utility_loss.json", "w") as f:
    json.dump(true_max_utility_loss, f)