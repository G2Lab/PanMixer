import numpy as np
import pandas as pd
import sys


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

blocks_path = BASE_PATH + "/starting_data/"

def compute_boundaries_and_blocks(chromosome):
    df = pd.read_csv(f"{blocks_path}/chr{chromosome}/blocks.blocks.det", delim_whitespace=True)
    boundaries_pos = []
    for _, row in df.iterrows():
        start = row['BP1']
        end = row['BP2']
        if start != end:
            boundaries_pos.append((start, end))
    
    positions = np.load(f"{blocks_path}/chr{chromosome}/pangenome_positions.npy")

    pop_blocks_idx = None   # Dictionary to store the block indices vector for each population
    pop_blocks_dict = None  # Dictionary to store the mapping from block index to (positions, indices) for each population

    # Convert the list of (start, end) tuples into a NumPy array and sort by start.
    intervals = np.array(boundaries_pos)
    intervals = intervals[np.argsort(intervals[:, 0])]

    # For each position, find the index of the interval with the greatest start <= position.
    idx = np.searchsorted(intervals[:, 0], positions, side='right') - 1

    # A position is in an interval if:
    #   (a) idx is not -1, and
    #   (b) the position is <= the corresponding end.
    valid = (idx >= 0) & (positions <= intervals[idx, 1])
    # Count the number of valid positions.
    #print(f"Found {len(intervals)} intervals with {len(positions)} positions and {np.sum(valid)} valid positions.")

    coverage = valid.sum()
    #print(f"{coverage} variants covered out of {len(positions)} ({coverage / len(positions) * 100:.2f}%)")

    # Create an empty array to hold block indices.
    blocks_idx = np.empty_like(positions, dtype=int)

    # For positions in a valid interval, assign the corresponding interval index.
    blocks_idx[valid] = idx[valid]
    #print(blocks_idx[valid].max())
    #print(blocks_idx)

    # For positions not covered by any interval, assign a unique block index.
    unique_offset = len(intervals)
    #print(unique_offset)
    blocks_idx[~valid] = np.arange(unique_offset, unique_offset + np.count_nonzero(~valid))

    #print(blocks_idx.max())

    # Store the block indices vector for this population.
    pop_blocks_idx = blocks_idx

    # Build a dictionary mapping each unique block index to a tuple:
    # (array of positions in that block, array of original indices for those positions)
    block_dict = {}
    unique_blocks = np.unique(blocks_idx)
    for block in unique_blocks:
        pos_indices = np.where(blocks_idx == block)[0]
        block_positions = positions[pos_indices]
        block_dict[block] = (block_positions, pos_indices)
    pop_blocks_dict = block_dict

    keys = np.array(list(pop_blocks_dict.keys()))
    #assert keys.shape[0] == np.max(keys) + 1, f"Block indices are not continuous. Found {keys.shape[0]} unique blocks, max index {np.max(keys)}"

    np.save(f"{blocks_path}/chr{chromosome}/simple_blocks_idx.npy", blocks_idx)

    import json

    pop_blocks_dict_serializable = {}

    serializable_blocks = {}
    for block, (block_positions, pos_indices) in pop_blocks_dict.items():
        # Convert numpy arrays to lists.
        block_positions_serial = [int(pos) for pos in block_positions]
        pos_indices_serial = [int(idx) for idx in pos_indices]
        serializable_blocks[int(block)] = (block_positions_serial, pos_indices_serial)
    pop_blocks_dict_serializable = serializable_blocks

    # Write the serializable dictionary to a JSON file.
    with open(f"{blocks_path}/chr{chromosome}/blocks_dict.json", "w") as f:
        json.dump(pop_blocks_dict_serializable, f, indent=4)

if __name__ == "__main__":
    # Example usage
    chromosome = int(sys.argv[1])
    compute_boundaries_and_blocks(chromosome)