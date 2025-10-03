#!/bin/bash
#SBATCH --job-name=pangenie_array
#SBATCH --output=logs/pangenie_array_%A_%a.out
#SBATCH --error=logs/pangenie_array_%A_%a.err
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00

to_run_vcfs=(
    "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_8/all-5.vcf.gz"
    "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_8/all-10.vcf.gz"
    "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_8/all-20.vcf.gz"
    "/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/experiments/exp_8/all-30.vcf.gz"
)

vcfbub_path=/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/tools/vcfbub

input_vcf=${to_run_vcfs[$SLURM_ARRAY_TASK_ID]}
output_vcf=${input_vcf%.vcf.gz}-bub.vcf
$vcfbub_path -l 0 -r 100000 --input $input_vcf > $output_vcf

module load bcftools
bgzip -f $output_vcf
bcftools index -f $output_vcf.gz

chr_rename_file=/gpfs/commons/groups/gursoy_lab/jblindenbach/Secret/PanMixer5/starting_data/rename_chr_map_pg.txt

output_vcf_renamed=${output_vcf%.vcf.gz}-renamed.vcf.gz

bcftools annotate --rename-chrs $chr_rename_file $output_vcf.gz -o $output_vcf_renamed
bcftools index -f $output_vcf_renamed

rm $output_vcf.gz
rm $output_vcf.gz.tbi