import numpy as np
import pickle
import sys
from tools.slurm_helper import launch_job_with_custom_command_by_rowid
from tools.utils import load_data_multichromosome
import os

from constants import (
    EXPERIMENT_PATH,
    READ_SUBJECTS,
    BASE_PATH,
    READS_DIR,
    PYTHON_ENV
)

def get_command_personalized_read_mapping(experiment_number, threads=16):
    return f"""

set -euo pipefail

export PATH={PYTHON_ENV}/bin:$PATH

vg_path={BASE_PATH}/downloaded_tools/vg

vcf_input={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID/new_haplotypes.vcf.gz
BASE=new_haplotypes
OUTDIR={EXPERIMENT_PATH}/exp_{experiment_number}/data/chr21/$ROW_ID

REF_FA={BASE_PATH}/starting_data/references/hg38_cleaned.fa
CHR_MAP={BASE_PATH}/starting_data/rename_chr_map_pg.txt

if [[ -s "${{OUTDIR}}/personalized_idx.hapl" ]]; then
  echo "personalized_idx.hapl already exists, skipping indexing."
else

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
      "$STEP4" "$STEP4".csi "$STEP4".tbi

echo "Done."

chr21_fasta={BASE_PATH}/starting_data/references/chr21.fa

rm -f ${{OUTDIR}}/personalized_idx*

$vg_path autoindex \
  --workflow giraffe \
  --prefix ${{OUTDIR}}/personalized_idx \
  --ref $chr21_fasta \
  --vcf $FINAL_GZ \
  -t {threads}

############################################
# Build haplotype index for personalized mapping
############################################
echo "[`date +'%F %T'`] Building r-index..."
$vg_path gbwt -Z \
  -r "${{OUTDIR}}/personalized_idx.ri" \
  "${{OUTDIR}}/personalized_idx.giraffe.gbz"

echo "[`date +'%F %T'`] Building haplotype index..."
$vg_path haplotypes -v 2 -t {threads} \
  -d "${{OUTDIR}}/personalized_idx.dist" \
  -r "${{OUTDIR}}/personalized_idx.ri" \
  -H "${{OUTDIR}}/personalized_idx.hapl" \
  "${{OUTDIR}}/personalized_idx.giraffe.gbz"

fi

############################################
# Read Mapping: vg giraffe (personalized with KFF)
############################################
export TMPDIR=${{OUTDIR}}/tmp_kmc
mkdir -p $TMPDIR

idx_prefix=${{OUTDIR}}/personalized_idx
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
  out_gam="${{OUTDIR}}/${{sample_name}}.personalized.giraffe.gam"
  stats_txt="${{OUTDIR}}/${{sample_name}}.personalized.giraffe.stats.txt"
  kff_file="${{TMPDIR}}/${{sample_name}}_${{ROW_ID}}"

  if [[ -s "$stats_txt" ]]; then
    echo "Skipping ${{sample_name}} (stats already exist)"
    continue
  fi

  if [[ ! -s "$fq" ]]; then
    echo "WARNING: Skipping ${{sample}} (missing or empty: $fq)" >&2
    continue
  fi

  echo "[`date +'%F %T'`] Running kmc for ${{sample_name}}..."
  kmc -k29 -m128 -okff -t{threads} -hp "$fq" "$kff_file" $TMPDIR

  echo "[`date +'%F %T'`] Running vg giraffe (personalized) on ${{sample}}..."
  $vg_path giraffe -p -t {threads} \
    -Z "${{idx_prefix}}.giraffe.gbz" \
    --haplotype-name "${{idx_prefix}}.hapl" \
    --kff-name "${{kff_file}}.kff" \
    -N "$sample_name" \
    -i \
    -f "$fq" \
    > "$out_gam"

  echo "[`date +'%F %T'`] Computing stats for ${{sample_name}}..."
  $vg_path stats -a "$out_gam" > "$stats_txt"

  # Cleanup
  rm -f "$out_gam"
  rm -f "${{kff_file}}.kff" "${{kff_file}}.kmc_pre" "${{kff_file}}.kmc_suf"

  echo "[`date +'%F %T'`] Done: ${{stats_txt}}"
done

rm -rf $TMPDIR
"""

def personalized_read_mapping(experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = len(data[1])

    id_counter = 0
    command = get_command_personalized_read_mapping(experiment_number, 16)
    return launch_job_with_custom_command_by_rowid("personalized_read_mapping", id_counter, command, memory="140g", cpus="16", num_tasks=str(num_tasks))
