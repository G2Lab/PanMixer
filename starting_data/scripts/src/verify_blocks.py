import sys
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

blocks_path = f"{BASE_PATH}/starting_data/"
    
chromosome = int(sys.argv[1])

positions_pangenome = np.load(f"{blocks_path}/chr{chromosome}/pangenome_positions.npy")
blocks_json = json.load(open(f"{blocks_path}/chr{chromosome}/blocks_dict.json", "r"))

all_positions = []
for j,key in enumerate(blocks_json):
    positions = blocks_json[key][0]
    all_positions.extend(positions)

print(len(all_positions))
print(positions_pangenome.shape)
