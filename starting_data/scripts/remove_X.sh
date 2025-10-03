INPUT="pangenome.vcf.gz"
OUTPUT="pangenome_no_X.vcf.gz"

bcftools view --regions grch38#chr1,grch38#chr2,grch38#chr3,grch38#chr4,grch38#chr5,grch38#chr6,grch38#chr7,grch38#chr8,grch38#chr9,grch38#chr10,grch38#chr11,grch38#chr12,grch38#chr13,grch38#chr14,grch38#chr15,grch38#chr16,grch38#chr17,grch38#chr18,grch38#chr19,grch38#chr20,grch38#chr21,grch38#chr22 "$INPUT" -Oz -o "$OUTPUT"
# Index the output VCF
bcftools index "$OUTPUT"
