set -euo pipefail

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

for CHR in $(seq 1 22); do
  src_vcf="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  src_tbi="${src_vcf}.tbi"

  dest_dir="chr${CHR}"
  dest_vcf="${dest_dir}/1000g_phased.vcf.gz"
  dest_tbi="${dest_dir}/1000g_phased.vcf.gz.tbi"

  mkdir -p "${dest_dir}"

  echo "[chr${CHR}] downloading VCF…"
  wget -c "${BASE_URL}/${src_vcf}" -O "${dest_vcf}"

  echo "[chr${CHR}] downloading index…"
  wget -c "${BASE_URL}/${src_tbi}" -O "${dest_tbi}"
done

echo "Done."