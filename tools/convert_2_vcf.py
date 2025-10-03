import pickle
import numpy as np
import sys
import sys
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data, add_command, get_chromosome_path, load_data_multichromosome
import os
import pandas as pd
import gzip

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
)

BATCHED_WRITES = 1000

def produce_vcf_from_np(genotypes, master_vcf_file_path, subject_name, output_file_path):
    #verify file_path exists
    if not os.path.exists(master_vcf_file_path):
        raise FileNotFoundError("File not found: " + master_vcf_file_path)

    genotypes = np.round(genotypes).astype(int)

    #check if all genotypes are empty
    is_empty = np.all(genotypes == -1)

    subject_index = -1

    open_func = gzip.open if master_vcf_file_path.endswith(".gz") else open
    with open_func(master_vcf_file_path, 'rt') as vcf:
        with open(output_file_path, 'w') as output:
            non_header_line_count = 0
            batched_lines = []
            for line in vcf:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        subject_index = line.strip().split('\t').index(subject_name)
                        if is_empty:
                            line_without_subject = line.strip().split('\t')
                            line_without_subject.pop(subject_index)
                            line = '\t'.join(line_without_subject) + '\n'
                    output.write(line)
                    continue
                assert subject_index != -1, "Subject not found in VCF header"
                if non_header_line_count % BATCHED_WRITES == 0:
                    output.writelines(batched_lines)
                    batched_lines = []
                    
                fields = line.strip().split('\t')
                site = int(fields[1])

                newline_string = '\t'.join(fields[:9])

                interested_field = subject_index

                replaced_field_string = ""
                if genotypes[non_header_line_count, 0] == -1:
                    replaced_field_string += '.'
                else:
                    replaced_field_string += str(genotypes[non_header_line_count, 0])
                replaced_field_string += '|'
                if genotypes[non_header_line_count, 1] == -1:
                    replaced_field_string += '.'
                else:
                    replaced_field_string += str(genotypes[non_header_line_count, 1])
                non_header_line_count += 1
                
                fields[interested_field] = replaced_field_string

                for i in range(9, len(fields)):
                    if "|" not in fields[i]:
                        fields[i] = ".|."
                if is_empty:
                    fields.pop(interested_field)

                newline_string += '\t' + '\t'.join(fields[9:])

                newline_string += '\n'

                batched_lines.append(newline_string)
            
            if len(batched_lines) > 0:
                output.writelines(batched_lines)

    assert subject_index != -1, "Subject not found in VCF header"
    assert non_header_line_count == genotypes.shape[0], f"Number of lines in VCF file {non_header_line_count} does not match number of genotypes {genotypes.shape[0]}"

def convert_2_vcf(experiment_number):

    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [experiment_number]
    launch_job_multichromosome("convert_2_vcf", args, memory="16g", cpus="1", num_tasks=str(num_tasks))

def main():
    experiment_number = sys.argv[1]
    chromosome = int(sys.argv[2])
    row_id = int(sys.argv[3])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]

    chromosome_path = STARTING_DATA_PATH + f"/chr{chromosome}/"
    master_vcf_file_path = chromosome_path + "pangenome.vcf.gz"

    interested_haplotypes = np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy")
    output_file_path = EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.vcf"
    produce_vcf_from_np(interested_haplotypes, master_vcf_file_path, subject_name, output_file_path)

    # compress and index the VCF file
    #os.system("module load bcftools")
    os.system(f"bgzip -f {output_file_path}")
    os.system(f"bcftools index {output_file_path}.gz")

if __name__ == '__main__':
    main()