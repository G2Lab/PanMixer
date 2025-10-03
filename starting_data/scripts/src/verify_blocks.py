import sys
import json
import numpy as np

blocks_path = "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/"
    
chromosome = int(sys.argv[1])

positions_pangenome = np.load(f"{blocks_path}/chr{chromosome}/pangenome_positions.npy")
blocks_json = json.load(open(f"{blocks_path}/chr{chromosome}/blocks_dict.json", "r"))

all_positions = []
for j,key in enumerate(blocks_json):
    positions = blocks_json[key][0]
    all_positions.extend(positions)

print(len(all_positions))
print(positions_pangenome.shape)
