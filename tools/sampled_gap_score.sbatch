#!/bin/bash
#SBATCH --job-name=sampled_gap
#SBATCH --array=0-43
#SBATCH --cpus-per-task=1
#SBATCH --mem=140G
#SBATCH --time=08:00:00
#SBATCH --output=slurm/sampled_gap_%a.out
#SBATCH --error=slurm/sampled_gap_%a.err

PYTHON=/users/jb4816/.conda/envs/genomics2/bin/python3

cd /groups/gg2845_gp/jb4816/PanMixer5

# Each task runs all 22 chromosomes for one subject
$PYTHON -c "
from tools.sampled_gap_score import score_sampled_haplotype
import json, sys

subject_id = $SLURM_ARRAY_TASK_ID
results = []
for chrom in range(1, 23):
    print(f'Scoring chr{chrom} subject {subject_id}...')
    r = score_sampled_haplotype(chrom, subject_id)
    print(f'  gap_score={r[\"genotype_gap_score\"]:.2f}')
    results.append(r)

with open(f'starting_data/sampled_gap_subject_{subject_id}.json', 'w') as f:
    json.dump(results, f, indent=2)

total_gap = sum(r['genotype_gap_score'] for r in results)
print(f'Total gap score across all chromosomes: {total_gap:.2f}')
"
