import os
import gzip
import numpy as np
import sys
from constants import STARTING_DATA_PATH, EXPERIMENT_PATH
from tools.slurm_helper import launch_job_multichromosome
import pandas as pd

# Assuming BATCHED_WRITES is defined globally or passed in context
BATCHED_WRITES = 1000 

def produce_vcf_from_np(genotypes_list, master_vcf_file_path, subject_names, output_file_path):
    """
    genotypes_list: List of numpy arrays (each N x 2).
    subject_names: List of strings matching the order of genotypes_list.
    """
    
    # 1. Validation and Setup
    if not os.path.exists(master_vcf_file_path):
        raise FileNotFoundError("File not found: " + master_vcf_file_path)
        
    if len(genotypes_list) != len(subject_names):
        raise ValueError("Length of genotypes_list and subject_names must match.")

    # 2. Pre-process Genotype Data
    # structure: { 'SubjectName': { 'data': np_array, 'empty': bool } }
    subject_data_map = {}
    
    for gt, name in zip(genotypes_list, subject_names):
        processed_gt = np.round(gt).astype(int)
        is_empty = np.all(processed_gt == -1)
        subject_data_map[name] = {
            'data': processed_gt,
            'empty': is_empty
        }

    open_func = gzip.open if master_vcf_file_path.endswith(".gz") else open
    
    # We will map VCF column indices to our subject data later
    # Format: { column_index: subject_data_entry }
    col_idx_to_data = {} 
    
    # Track if we found all subjects in the header
    found_subjects_count = 0

    with open_func(master_vcf_file_path, 'rt') as vcf:
        with open(output_file_path, 'w') as output:
            non_header_line_count = 0
            batched_lines = []
            
            for line in vcf:
                # --- Header Processing ---
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        header_parts = line.strip().split('\t')
                        
                        # Identify which column belongs to which subject
                        for name, data in subject_data_map.items():
                            if name in header_parts:
                                idx = header_parts.index(name)
                                col_idx_to_data[idx] = data
                                found_subjects_count += 1
                        
                        # Rebuild header, skipping subjects marked as 'empty'
                        new_header_parts = []
                        for i, part in enumerate(header_parts):
                            # If this column is one of our targets AND it is empty, skip it (remove)
                            if i in col_idx_to_data and col_idx_to_data[i]['empty']:
                                continue
                            new_header_parts.append(part)
                            
                        line = '\t'.join(new_header_parts) + '\n'
                        
                    output.write(line)
                    continue

                # --- Body Processing ---
                # Ensure we found all requested subjects before processing body
                if found_subjects_count != len(subject_names):
                     # Determine which were missing for a helpful error message
                     found_names = [d['name'] for d in col_idx_to_data.values() if 'name' in d] # simplified check
                     missing = set(subject_names) - set(subject_data_map.keys()) # Logic implies we verify against mapped indices
                     # Note: A stricter check can be added here, but asserting on found count is safe.
                     assert found_subjects_count == len(subject_names), "Not all subjects found in VCF header"

                if non_header_line_count % BATCHED_WRITES == 0 and non_header_line_count > 0:
                    output.writelines(batched_lines)
                    batched_lines = []
                
                fields = line.strip().split('\t')
                
                # Start constructing the new line. 
                # First 9 columns (CHROM...FORMAT) are standard and kept as-is.
                new_fields = fields[:9]
                
                # Iterate through sample columns (index 9 onwards)
                for i in range(9, len(fields)):
                    
                    # CASE A: This is one of our target subjects
                    if i in col_idx_to_data:
                        subj_info = col_idx_to_data[i]
                        
                        # If empty, we skip appending (effectively removing the column)
                        if subj_info['empty']:
                            continue
                        
                        # Otherwise, replace data using the pre-loaded numpy array
                        g_data = subj_info['data'][non_header_line_count]
                        
                        # Construct GT string (e.g., "0|1" or ".|.")
                        allele_1 = str(g_data[0]) if g_data[0] != -1 else '.'
                        allele_2 = str(g_data[1]) if g_data[1] != -1 else '.'
                        replaced_field_string = f"{allele_1}|{allele_2}"
                        
                        new_fields.append(replaced_field_string)

                    # CASE B: This is a subject we are NOT touching
                    else:
                        current_field = fields[i]
                        # Apply the cleanup logic from original script
                        if "|" not in current_field:
                            current_field = ".|."
                        new_fields.append(current_field)

                # Construct final line
                newline_string = '\t'.join(new_fields) + '\n'
                batched_lines.append(newline_string)
                
                non_header_line_count += 1
            
            # Write remaining lines
            if batched_lines:
                output.writelines(batched_lines)

    # 3. Final Validation
    # Ensure all genotype arrays matched the VCF line count
    for name, data in subject_data_map.items():
        arr_len = data['data'].shape[0]
        assert non_header_line_count == arr_len, \
            f"Subject {name}: VCF lines ({non_header_line_count}) != Genotypes ({arr_len})"

def load_multitarget_csv(target_experiment_number):
    return pd.read_csv(EXPERIMENT_PATH + f"/exp_{target_experiment_number}/data.csv")

def create_multitarget_vcfs(experiment_number, target_experiment_number):
    data = load_multitarget_csv(target_experiment_number)

    num_tasks = len(data) * 22

    args = [experiment_number, target_experiment_number]
    launch_job_multichromosome("create_multitarget_vcfs", args, memory="16g", cpus="1", num_tasks=str(num_tasks))

def main():
    experiment_number = int(sys.argv[1])
    target_experiment_number = int(sys.argv[2])
    chromosome = int(sys.argv[3])
    row_id = int(sys.argv[4])

    target_individuals_list = pd.read_csv(EXPERIMENT_PATH + f"/exp_{target_experiment_number}/data.csv")
    target_individuals_list = target_individuals_list.iloc[row_id]["target_individuals"]
    target_individuals_list = str(target_individuals_list)
    target_individuals_list = target_individuals_list.split(",")
    target_individuals_list = [x.strip() for x in target_individuals_list]
    target_individuals_list = [int(x) for x in target_individuals_list]

    genotypes_list = []
    for i in target_individuals_list:
        genotypes_list.append(np.load(EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{i}/new_haplotypes.npy"))
    
    subjects_list = []
    subjects_all = np.loadtxt(STARTING_DATA_PATH + "/subjects_files/subjects.txt", dtype=str)
    
    for i in target_individuals_list:
        subjects_list.append(subjects_all[i])
        
    master_vcf_file_path = STARTING_DATA_PATH + f"/chr{chromosome}/pangenome.vcf.gz"
    output_file_path = EXPERIMENT_PATH + f"/exp_{target_experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.vcf"

    produce_vcf_from_np(genotypes_list, master_vcf_file_path, subjects_list, output_file_path)

    # compress and index the VCF file
    #os.system("module load bcftools")
    os.system(f"bgzip -f {output_file_path}")
    os.system(f"bcftools index {output_file_path}.gz")


if __name__ == "__main__":
    main()