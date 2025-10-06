./scripts/get_pangenomes.sh
./scripts/get_pangenome_alignments.sh
./scripts/get_1000g_phased.sh

bcftools index pangenome.vcf.gz

./scripts/remove_X.sh
./scripts/remove_chm13.sh

sbatch scripts/split_data.sbatch