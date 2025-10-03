import numpy as np
import sys

def npy_to_vcf(npy_file, subjects_npy_file, vcf_file):
    # Load genotype and subject data
    genotypes = np.load(npy_file)  # Shape: (num_samples, num_variants, 2)
    subjects = np.load(subjects_npy_file)  # Shape: (num_samples,)

    # Transpose genotypes to shape: (num_variants, num_samples, 2)
    genotypes = genotypes.transpose(1, 0, 2)

    # Read the original VCF file and separate header and data
    with open(vcf_file, 'r') as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.startswith('#')]
    data_lines = [line for line in lines if not line.startswith('#')]

    # Modify the #CHROM header line to include subject IDs
    for i, line in enumerate(header_lines):
        if line.startswith('#CHROM'):
            parts = line.strip().split('\t')
            header_lines[i] = '\t'.join(parts[:9] + subjects.tolist()) + '\n'
            break

    # Create output VCF
    output_vcf = vcf_file.replace('.vcf', '_output.vcf')
    with open(output_vcf, 'w') as f:
        # Write headers
        f.writelines(header_lines)

        # Write each variant line with genotype info
        for i, line in enumerate(data_lines):
            parts = line.strip().split('\t')
            fixed_fields = parts[:9]

            # Extract genotypes for this variant across all subjects
            gt_line = []
            for sample_gt in genotypes[i]:
                a1, a2 = sample_gt
                if a1 == -1 or a2 == -1:  # Missing data
                    gt_line.append("./.")
                else:
                    gt_line.append(f"{a1}/{a2}")

            f.write('\t'.join(fixed_fields + gt_line) + '\n')


if __name__ == "__main__":
    npy_file = sys.argv[1]
    subjects_file = npy_file.replace(".npy", "_subjects.npy")
    vcf_file = sys.argv[2]
    npy_to_vcf(npy_file, subjects_file, vcf_file)