#!/bin/bash
# get_blocks_grch38.sh <CHR>
#
# For one chromosome:
#   1. Download the GRCh38 1000G 30x phased panel (if not already present).
#   2. Subset to SNPs that are ALSO present in the HPRC pangenome graph
#      (match on GRCh38 position + REF + ALT, after splitting multiallelics).
#   3. Regenerate LD blocks with PLINK 1.9 (same command as the original pipeline).
#
# All inputs/outputs live under starting_data/chr<CHR>/.
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <chromosome 1-22>" >&2
  exit 2
fi

CHR="$1"
if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "Chromosome must be an integer from 1 through 22: $CHR" >&2
  exit 2
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STARTING_DATA_DIR="$(cd "$HERE/../.." && pwd)"

if ! command -v bcftools >/dev/null 2>&1; then
  source /apps/ohpc/pub/apps/conda/3/etc/profile.d/conda.sh
  conda activate genomics2
fi
module load plink/1.9 2>/dev/null || true

for executable in bcftools plink wget; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    echo "Required executable not found: $executable" >&2
    exit 1
  fi
done

X30_BASE_URL="${X30_BASE_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV}"
x30_url() {
  echo "${X30_BASE_URL}/1kGP_high_coverage_Illumina.chr${1}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
}

OUT="$STARTING_DATA_DIR/chr${CHR}"
mkdir -p "$OUT"

X30_VCF="$OUT/1000g_30x_phased.vcf.gz"
PANGENOME="$OUT/pangenome.vcf.gz"

# ---------------------------------------------------------------------------
# 1. Download the 30x GRCh38 panel (resumable; skip if a complete copy exists)
# ---------------------------------------------------------------------------
if [ ! -s "$X30_VCF" ] || ! bcftools index -n "$X30_VCF" >/dev/null 2>&1; then
  echo "[chr${CHR}] downloading 30x panel ..."
  wget -nv -c "$(x30_url "$CHR")"      -O "$X30_VCF"
  wget -nv -c "$(x30_url "$CHR").tbi"  -O "$X30_VCF.tbi"
else
  echo "[chr${CHR}] 30x panel already present, skipping download"
fi

# ---------------------------------------------------------------------------
# 2. Build SNP-only, same-contig-named views and intersect
# ---------------------------------------------------------------------------
# Pangenome contig is grch38#chrN; 30x contig is chrN. Rename pangenome to chrN
# so isec matches on identical CHROM/POS/REF/ALT.
echo "grch38#chr${CHR} chr${CHR}" > "$OUT/rename_chrs.txt"

PG_SNPS=$OUT/pangenome_snps.vcf.gz
X30_SNPS=$OUT/30x_snps.vcf.gz
SUBSET=$OUT/1000g_30x_pangenome_snps.vcf.gz

echo "[chr${CHR}] extracting pangenome SNP sites ..."
bcftools view -r "grch38#chr${CHR}" "$PANGENOME" -Ou \
  | bcftools norm -m-any -Ou \
  | bcftools view -v snps -Ou \
  | bcftools annotate --rename-chrs "$OUT/rename_chrs.txt" -Ou \
  | bcftools sort -Oz -o "$PG_SNPS"
bcftools index -f "$PG_SNPS"

echo "[chr${CHR}] extracting 30x SNPs ..."
bcftools norm -m-any "$X30_VCF" -Ou \
  | bcftools view -v snps -Oz -o "$X30_SNPS"
bcftools index -f "$X30_SNPS"

echo "[chr${CHR}] intersecting 30x SNPs with pangenome SNPs ..."
# -n=2 -w1 : write file-1 (30x) records that are present in BOTH.
# -c none  : require exact CHROM/POS/REF/ALT identity.
bcftools isec -c none -n=2 -w1 "$X30_SNPS" "$PG_SNPS" -Oz -o "$SUBSET"
bcftools index -f "$SUBSET"

n_30x=$(bcftools view -H "$X30_SNPS"  | wc -l)
n_pg=$(bcftools view -H "$PG_SNPS"    | wc -l)
n_sub=$(bcftools view -H "$SUBSET"    | wc -l)
echo "[chr${CHR}] 30x SNPs=$n_30x  pangenome SNPs=$n_pg  shared subset=$n_sub"

# ---------------------------------------------------------------------------
# 3. Regenerate LD blocks with PLINK 1.9 (same invocation as the original)
# ---------------------------------------------------------------------------
echo "[chr${CHR}] running PLINK --make-bed ..."
plink --vcf "$SUBSET" --make-bed --out "$OUT/subset" \
      --allow-extra-chr --vcf-half-call missing --double-id

echo "[chr${CHR}] running PLINK --blocks ..."
plink --blocks no-pheno-req --bfile "$OUT/subset" --out "$OUT/blocks" --allow-extra-chr

echo "[chr${CHR}] DONE. blocks: $OUT/blocks.blocks(.det)"
