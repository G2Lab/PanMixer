import gzip
import sys
from pathlib import Path

import numpy as np


STARTING_DATA_PATH = Path(__file__).resolve().parents[2]
BASES = {"A", "C", "G", "T"}


def is_biallelic_snp(ref, alt_field):
    alts = alt_field.split(",")
    return len(alts) == 1 and ref.upper() in BASES and alts[0].upper() in BASES


def build_biallelic_snp_mask(chromosome):
    chromosome_path = STARTING_DATA_PATH / f"chr{chromosome}"
    vcf_path = chromosome_path / "pangenome.vcf.gz"
    positions_path = chromosome_path / "pangenome_positions.npy"
    output_path = chromosome_path / "biallelic_snp_mask.npy"

    if not vcf_path.exists():
        raise FileNotFoundError(vcf_path)
    if not positions_path.exists():
        raise FileNotFoundError(positions_path)

    positions = np.load(positions_path, mmap_mode="r")
    mask = np.zeros(positions.shape[0], dtype=bool)

    record_count = 0
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            if record_count >= len(mask):
                raise ValueError(
                    f"{vcf_path} has more records than {positions_path} "
                    f"({record_count + 1} > {len(mask)})"
                )

            fields = line.rstrip("\n").split("\t", 5)
            pos = int(fields[1])
            if pos != int(positions[record_count]):
                raise ValueError(
                    f"{vcf_path} record {record_count} has POS={pos}, "
                    f"but pangenome_positions.npy has {int(positions[record_count])}"
                )

            mask[record_count] = is_biallelic_snp(fields[3], fields[4])
            record_count += 1

    if record_count != len(mask):
        raise ValueError(
            f"{vcf_path} has {record_count} records, but {positions_path} has {len(mask)}"
        )

    np.save(output_path, mask)
    print(
        f"Saved {output_path} with {int(mask.sum())} bi-allelic SNPs "
        f"out of {len(mask)} pangenome records"
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python build_biallelic_snp_mask.py <chromosome>")
        sys.exit(1)

    build_biallelic_snp_mask(sys.argv[1])
