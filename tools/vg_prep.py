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
    BASE_PATH
)

def get_command_vg_prep(experiment_number, threads=16):
    return f"""
set -euo pipefail

module load bcftools
module load tabixpp/1.1.2-GCC-12.3.0

vg_path={BASE_PATH}/downloaded_tools/vg

vcf_input={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID/new_haplotypes.vcf.gz
BASE=new_haplotypes
OUTDIR={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID

REF_FA={BASE_PATH}/starting_data/references/hg38_cleaned.fa
CHR_MAP={BASE_PATH}/starting_data/rename_chr_map_pg.txt

############################################
# Step 1: view -c 1
############################################
STEP1="${{OUTDIR}}/${{BASE}}.step1.min1sample.vcf.gz"
echo "[1/4] bcftools view -c 1 ..."
bcftools view -c 1 -Oz -o "$STEP1" "$vcf_input"
bcftools index -f -c "$STEP1"

############################################
# Step 2: norm -N -f REF
############################################
STEP2="${{OUTDIR}}/${{BASE}}.step2.norm.vcf.gz"
echo "[2/4] bcftools norm -N -f $REF_FA ..."
bcftools norm -N -f "$REF_FA" -Oz -o "$STEP2" "$STEP1"
bcftools index -f -c "$STEP2"

############################################
# Step 3: view --trim-alt-alleles
############################################
STEP3="${{OUTDIR}}/${{BASE}}.step3.trimalt.vcf.gz"
echo "[3/4] bcftools view --trim-alt-alleles ..."
bcftools view --trim-alt-alleles -Oz -o "$STEP3" "$STEP2"
bcftools index -f -c "$STEP3"

############################################
# Step 4: annotate --rename-chrs
############################################
STEP4="${{OUTDIR}}/${{BASE}}.step4.vcf.gz"
echo "[4/4] bcftools annotate --rename-chrs ..."
bcftools annotate --rename-chrs "$CHR_MAP" -Oz -o "$STEP4" "$STEP3"
bcftools index -f -c "$STEP4"

############################################
# Step 5: sort bgzip
############################################
# remove artifacts
rm -f ${{OUTDIR}}/${{BASE}}.final.vcf*


FINAL="${{OUTDIR}}/${{BASE}}.final.vcf"
FINAL_GZ="${{OUTDIR}}/${{BASE}}.final.vcf.gz"
bcftools view -r chr21 -v snps -o $FINAL $STEP4

bcftools sort $FINAL -Oz -o $FINAL_GZ
tabix -p vcf -f "$FINAL_GZ"

############################################
# Cleanup
############################################
rm -f "$STEP1" "$STEP1".csi "$STEP1".tbi \
          "$STEP2" "$STEP2".csi "$STEP2".tbi \
          "$STEP3" "$STEP3".csi "$STEP3".tbi \
          "$STEP4" "STEP4.csi" "$STEP4.tbi" 

echo "Done."

chr21_fasta={BASE_PATH}/starting_data/references/chr21.fa

rm -f ${{OUTDIR}}/vg_idx*

$vg_path autoindex \
  --workflow giraffe \
  --prefix ${{OUTDIR}}/vg_idx \
  --ref $chr21_fasta \
  --vcf $FINAL_GZ \
  -t {threads}
"""

def vg_prep(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    id_counter = 0
    command = get_command_vg_prep(experiment_number, 16)
    launch_job_with_custom_command_by_rowid("vg_prep", id_counter, command, memory="32g", cpus="16", num_tasks=str(num_tasks))
