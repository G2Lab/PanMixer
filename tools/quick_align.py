import numpy as np
import pickle
import sys
import sys
from tools.slurm_helper import launch_job_with_custom_command_by_rowid
from tools.utils import load_data_multichromosome
import os

from constants import (
    EXPERIMENT_PATH,
    READ_SUBJECTS,
    BASE_PATH,
    READS_DIR
)

def get_command_quick_align(experiment_number, threads=16):
    return f"""
set -euo pipefail

vg_path={BASE_PATH}/downloaded_tools/vg

OUTDIR={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID
idx_prefix=${{OUTDIR}}/vg_idx
base_sample_path={READS_DIR}

samples=(
  HG00138.fastq
  HG00635.fastq
  HG01112.fastq
  HG01600.fastq
  HG02698.fastq
  NA12778.fastq
  NA18853.fastq
)

for sample in "${{samples[@]}}"; do
  sample_name="${{sample%.fastq}}"
  fq="${{base_sample_path}}/${{sample}}"
  out_gam="${{OUTDIR}}/${{sample_name}}.giraffe.gam"
  stats_txt="${{OUTDIR}}/${{sample_name}}.giraffe.stats.txt"

  # Sanity check
  if [[ ! -s "$fq" ]]; then
    echo "WARNING: Skipping ${{sample}} (missing or empty: $fq)" >&2
    continue
  fi

  echo "[`date +'%F %T'`] Running vg giraffe on ${{sample}}..."
  $vg_path giraffe \
    -Z "${{idx_prefix}}.giraffe.gbz" \
    -z "${{idx_prefix}}.shortread.zipcodes" \
    -m "${{idx_prefix}}.shortread.withzip.min" \
    -d "${{idx_prefix}}.dist" \
    -f "$fq" \
    -t {threads} \
    > "$out_gam"

  echo "[`date +'%F %T'`] Computing stats for ${{sample_name}}..."
  $vg_path stats -a "$out_gam" > "$stats_txt"

  # Keep only the stats
  rm -f "$out_gam"

  echo "[`date +'%F %T'`] Done: ${{stats_txt}}"
done
"""

def quick_align(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    print("Num tasks", num_tasks)

    id_counter = 0
    command = get_command_quick_align(experiment_number)
    launch_job_with_custom_command_by_rowid("quick_align", id_counter, command, memory="32g", cpus="1", num_tasks=str(num_tasks))