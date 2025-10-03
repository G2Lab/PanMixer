for CHR in {1..22}; do
	cp /gpfs/commons/datasets/1000genomes/hg38/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz \
		./chr${CHR}/1000g_phased.vcf.gz
	cp /gpfs/commons/datasets/1000genomes/hg38/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz.tbi \
		./chr${CHR}/1000g_phased.vcf.gz.tbi
done