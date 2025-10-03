#!/bin/bash
#SBATCH --job-name=get_beagle_ready
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --output=get_beagle_ready%a.out
#SBATCH --error=get_beagle_ready%a.err

CHR=$SLURM_ARRAY_TASK_ID
CHR_DIR="./chr${CHR}"

module load bcftools

out_no_missing="${CHR_DIR}/1000g_snps_no_missing.vcf.gz"
out_no_missing_renamed="${CHR_DIR}/1000g_snps_no_missing_renamed.vcf.gz"
out_no_missing_renamed_norm="${CHR_DIR}/1000g_snps_no_missing_renamed_norm.vcf.gz"

bcftools view -g ^miss -Oz -o $out_no_missing ${CHR_DIR}/1000g_phased.vcf.gz
bcftools index $out_no_missing

rename_chr_file="/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/rename_chr.txt"

bcftools annotate --rename-chrs $rename_chr_file -Oz -o $out_no_missing_renamed ${out_no_missing}
bcftools index $out_no_missing_renamed

bcftools norm --rm-dup all -Oz -o $out_no_missing_renamed_norm ${out_no_missing_renamed}
bcftools index $out_no_missing_renamed_norm

rm $out_no_missing
rm $out_no_missing_renamed
rm $out_no_missing.tbi
rm $out_no_missing_renamed.tbi
