./scripts/get_pangenomes.sh
./scripts/get_pangenome_alignments.sh
./scripts/get_1000g_phased.sh

bcftools index pangenome.vcf.gz

./scripts/remove_X.sh
./scripts/remove_chm13.sh

sbatch scripts/split_data.sbatch

sbatch scripts/get_num_alleles.sbatch
sbatch scripts/get_blocks.sbatch
sbatch scripts/convert_2_npy.sbatch
sbatch scripts/get_mappings.sbatch
sbatch scripts/segment_blocks.sbatch
sbatch scripts/get_pmi_utility.sbatch