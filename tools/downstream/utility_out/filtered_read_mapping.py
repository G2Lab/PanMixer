import numpy as np
import pickle
import sys
from tools.common.slurm_helper import launch_job_with_custom_command_by_rowid
from tools.common.utils import load_data_multichromosome
import os

from constants import (
    EXPERIMENT_PATH,
    READ_SUBJECTS,
    BASE_PATH,
    READS_DIR,
)

def get_command_filtered_read_mapping(experiment_number, threads=16):
    return f"""

set -euo pipefail

vg_path={BASE_PATH}/downloaded_tools/vg

vcf_input={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID/new_haplotypes.vcf.gz
BASE=new_haplotypes
OUTDIR={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID

REF_FA={BASE_PATH}/starting_data/references/hg38_cleaned.fa
CHR_MAP={BASE_PATH}/starting_data/rename_chr_map_pg.txt

if [[ -s "${{OUTDIR}}/filter_idx.giraffe.gbz" ]]; then
  echo "filter_idx.giraffe.gbz already exists, skipping indexing."
else

############################################
# Step 1: view -c 1
############################################
STEP1="${{OUTDIR}}/${{BASE}}.step1.min1sample.vcf.gz"
echo "[1/5] bcftools view -c 1 ..."
bcftools view -c 1 -Oz -o "$STEP1" "$vcf_input"
bcftools index -f -c "$STEP1"

############################################
# Step 2: norm -N -f REF
############################################
STEP2="${{OUTDIR}}/${{BASE}}.step2.norm.vcf.gz"
echo "[2/5] bcftools norm -N -f $REF_FA ..."
bcftools norm -N -f "$REF_FA" -Oz -o "$STEP2" "$STEP1"
bcftools index -f -c "$STEP2"

############################################
# Step 3: Filter AF < 10% (Added)
############################################
STEP3="${{OUTDIR}}/${{BASE}}.step3.filter_af.vcf.gz"
echo "[3/5] bcftools view -i 'AF>=0.1' ..."
# Filters out any site where the alternative allele frequency is less than 0.1
bcftools view -i 'AF>=0.1' -Oz -o "$STEP3" "$STEP2"
bcftools index -f -c "$STEP3"

############################################
# Step 4: view --trim-alt-alleles
############################################
STEP4="${{OUTDIR}}/${{BASE}}.step4.trimalt.vcf.gz"
echo "[4/5] bcftools view --trim-alt-alleles ..."
bcftools view --trim-alt-alleles -Oz -o "$STEP4" "$STEP3"
bcftools index -f -c "$STEP4"

############################################
# Step 5: annotate --rename-chrs
############################################
STEP5="${{OUTDIR}}/${{BASE}}.step5.vcf.gz"
echo "[5/5] bcftools annotate --rename-chrs ..."
bcftools annotate --rename-chrs "$CHR_MAP" -Oz -o "$STEP5" "$STEP4"
bcftools index -f -c "$STEP5"

############################################
# Step 6: sort bgzip
############################################
# remove artifacts
rm -f ${{OUTDIR}}/${{BASE}}.final.vcf*

FINAL="${{OUTDIR}}/${{BASE}}.final.vcf"
FINAL_GZ="${{OUTDIR}}/${{BASE}}.final.vcf.gz"
bcftools view -r chr21 -v snps -o $FINAL $STEP5

bcftools sort $FINAL -Oz -o $FINAL_GZ
tabix -p vcf -f "$FINAL_GZ"

############################################
# Cleanup
############################################
rm -f "$STEP1" "$STEP1".csi "$STEP1".tbi \
      "$STEP2" "$STEP2".csi "$STEP2".tbi \
      "$STEP3" "$STEP3".csi "$STEP3".tbi \
      "$STEP4" "$STEP4".csi "$STEP4".tbi \
      "$STEP5" "$STEP5".csi "$STEP5".tbi

echo "Done."

chr21_fasta={BASE_PATH}/starting_data/references/chr21.fa

rm -f ${{OUTDIR}}/filter_idx*

$vg_path autoindex \
  --workflow giraffe \
  --prefix ${{OUTDIR}}/filter_idx \
  --ref $chr21_fasta \
  --vcf $FINAL_GZ \
  -t {threads}

fi

############################################
# Read Mapping: vg giraffe
############################################
idx_prefix=${{OUTDIR}}/filter_idx
base_sample_path={READS_DIR}

samples=(
  HG00138.fastq
  HG00635.fastq
  HG01112.fastq
  HG01600.fastq
  HG02698.fastq
)

for sample in "${{samples[@]}}"; do
  sample_name="${{sample%.fastq}}"
  fq="${{base_sample_path}}/${{sample}}"
  out_gam="${{OUTDIR}}/${{sample_name}}.filtered.giraffe.gam"
  stats_txt="${{OUTDIR}}/${{sample_name}}.filtered.giraffe.stats.txt"

  if [[ -s "$stats_txt" ]]; then
    echo "Skipping ${{sample_name}} (stats already exist)"
    continue
  fi

  if [[ ! -s "$fq" ]]; then
    echo "WARNING: Skipping ${{sample}} (missing or empty: $fq)" >&2
    continue
  fi

  echo "[`date +'%F %T'`] Running vg giraffe (filtered) on ${{sample}}..."
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

def filtered_read_mapping(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    id_counter = 0
    command = get_command_filtered_read_mapping(experiment_number, 16)
    return launch_job_with_custom_command_by_rowid("filtered_read_mapping", id_counter, command, memory="32g", cpus="16", num_tasks=str(num_tasks))
