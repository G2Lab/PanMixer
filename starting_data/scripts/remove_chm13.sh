INPUT="pangenome_no_X.vcf.gz"
OUTPUT="pangenome_no_X_no_chm13.vcf.gz"

bcftools view -s ^chm13 "$INPUT" -Oz -o "$OUTPUT"
# Index the output VCF
bcftools index "$OUTPUT"
