import numpy as np
import scipy.stats as stats
from tools.common.utils import load_data, load_data_multichromosome, str_to_bool
from tools.common.slurm_helper import launch_job_multichromosome
from scipy.stats import pearsonr


KB_WINDOW_SIZE = 5
PROGRESS_INDEX = 1000

from constants import (
    STARTING_DATA_PATH,
    EXPERIMENT_PATH,
)

def compute_R2(x, ys):
    """
    Compute r^2 values between a single x and multiple ys.
    
    Parameters:
    - x: [n_samples, 2] or [n_samples] (numeric)
    - ys: [n_samples, n_y, 2] or [n_samples, n_y]
    
    Returns:
    - r2_values: [n_y] array of r^2 values
    """
    # collapse haplotype dimension if present
    D, p_A, p_Bs = compute_D(x, ys)
    D_squared = np.square(D)
    denom = p_A * (1 - p_A) * p_Bs * (1 - p_Bs)
    if np.any(denom == 0):
        denom[denom == 0] = np.nan
    return D_squared / denom

def compute_D(x, ys):
    """
    Compute D between a single x and multiple ys using vectorized NumPy operations.
    
    Parameters:
    - x: [n_samples, 2]
    - ys: [n_samples, n_y, 2]
    
    Returns:
    - D_values: [n_y] array of D values
    """
    n_samples, n_y, _ = ys.shape
    x = x.astype(float)
    ys = ys.astype(float)
    x[x == -1] = np.nan
    ys[ys == -1] = np.nan

    # Compute frequencies
    p_A = np.nanmean(x)
    p_Bs = np.nanmean(ys, axis=(0, 2))
    xs = np.broadcast_to(x[:, np.newaxis], (n_samples, n_y, 2))

    p_AB = np.nanmean((xs == 1) & (ys == 1), axis=(0, 2))

    D = p_AB - p_A * p_Bs
    return D, p_A, p_Bs


def compute_ld_loss(master_npy_file, master_npy_position_file, other_npy_file, master_npy_subjects_file, snp_positions_file, subject, dont_replace):
    # Load the numpy arrays
    master_positions = np.load(master_npy_position_file)
    snp_positions = np.load(snp_positions_file)

    snp_mask = np.zeros(len(master_positions), dtype=bool)
    snp_mask[snp_positions] = True

    assert len(snp_mask) == len(master_positions), "SNP positions and master positions length mismatch"

    master_array = np.load(master_npy_file)
    other_array = np.load(other_npy_file)

    master_array = master_array[:, snp_mask]
    print(not dont_replace, dont_replace)
    if not dont_replace:
        obfuscated_subject = other_array[snp_mask]
        master_positions = master_positions[snp_mask]

        master_subjects = np.load(master_npy_subjects_file)
        subject_index = np.where(master_subjects == subject)[0][0]
        assert subject_index != -1, f"Subject {subject} not found in master subjects"
        print(f"Subject index: {subject_index}")
        
        other_array = np.copy(master_array)
        other_array[subject_index, :] = obfuscated_subject
    else:
        master_positions = master_positions[snp_mask]
        other_array = other_array[:, snp_mask]

    print(master_array.shape)
    print(other_array.shape)

    l1_loss = 0.0
    count = 0

    for i,pos in enumerate(master_positions):
        if i % PROGRESS_INDEX == 0:
            print(f"Computing LD for position {i} out of {len(master_positions)}")
            print(f"Current LD Loss: {l1_loss / count if count > 0 else 0.0}")
        delta_pos = 1
        master_x = master_array[:, i]
        other_x = other_array[:, i]
        ld_indices = []
        while i + delta_pos < len(master_positions) and master_positions[i + delta_pos] - pos < KB_WINDOW_SIZE: 
            ld_indices.append(i + delta_pos)
            delta_pos += 1
            
        if len(ld_indices) == 0:
            continue
        master_y = master_array[:, ld_indices]
        other_y = other_array[:, ld_indices]
        
        master_r2 = compute_R2(master_x, master_y)
        other_r2 = compute_R2(other_x, other_y)

        l1_loss += np.nansum(np.abs(master_r2 - other_r2))
        count += np.sum(~np.isnan(master_r2) & ~np.isnan(other_r2))
        delta_pos += 1

    return (l1_loss, count) if count > 0 else (0.0, 0)
    

def ld_loss_main(experiment_number, chromosome, row_id, dont_replace):
    chromosome_path = STARTING_DATA_PATH + "/chr" + str(chromosome) + "/"
    master_npy_file_positions = chromosome_path + "pangenome_positions.npy"
    master_npy_file_subjects = chromosome_path + "pangenome_subjects.npy"
    master_npy_file = chromosome_path + "pangenome.npy"

    other_npy_file = EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/new_haplotypes.npy"
    output_file = EXPERIMENT_PATH + f"/exp_{experiment_number}/data/chr{chromosome}/{row_id}/ld_loss.txt"
    
    if not dont_replace:
        data_df = load_data(experiment_number, chromosome)
        data = data_df.iloc[row_id]
        subject_name = data["subject"]
    else:
        subject_name = ""

    snp_position_file = chromosome_path + "pangenome_mask.npy"
    ld_loss, count = compute_ld_loss(master_npy_file, master_npy_file_positions, other_npy_file, master_npy_file_subjects, snp_position_file, subject_name, dont_replace)
    print("LD Loss:", ld_loss)
    with open(output_file, 'w') as f:
        f.write(f"{ld_loss}\n")
        f.write(f"{count}\n")

def ld_loss(experiment_number, dont_replace=False):
    data_df = load_data_multichromosome(experiment_number)
    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data_df[i])
    
    args = [experiment_number, dont_replace]
    return launch_job_multichromosome("ld_loss", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

if __name__ == "__main__":
    import sys
    experiment_number = sys.argv[1]
    dont_replace = str_to_bool(sys.argv[2])
    chromosome = sys.argv[3]
    row_id = int(sys.argv[4])
    
    ld_loss_main(experiment_number, chromosome, row_id, dont_replace)