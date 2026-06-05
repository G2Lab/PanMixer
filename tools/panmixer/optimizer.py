import numpy as np
import pickle
import sys
from ortools.linear_solver import pywraplp
import sys
from tools.common.slurm_helper import launch_job_multichromosome
from tools.common.utils import load_data_multichromosome, load_data
import json

from constants import EXPERIMENT_PATH, VERBOSE, STARTING_DATA_PATH


def apply_baseline_unique(pangenome, subject_id):
    """
    Remove unique alleles from the target individual.
    
    If the target individual is the only person with a particular allele at a position,
    set both haplotypes at that position to -1 (empty).
    
    Args:
        pangenome: numpy array of shape (num_subjects, num_positions, 2)
        subject_id: index of the target individual
        
    Returns:
        modified_haplotypes: numpy array of shape (num_positions, 2) with unique alleles set to -1
    """
    num_subjects = pangenome.shape[0]
    num_positions = pangenome.shape[1]
    
    # Get the target individual's haplotypes
    target_haplotypes = pangenome[subject_id].copy()
    
    # Create a mask of all other subjects (excluding target)
    other_subjects_mask = np.ones(num_subjects, dtype=bool)
    other_subjects_mask[subject_id] = False
    other_pangenome = pangenome[other_subjects_mask]  # Shape: (num_subjects-1, num_positions, 2)
    
    # For each position, check if target's alleles are unique
    mask = np.zeros((num_positions, 2), dtype=int)
    for pos in range(num_positions):
        for hap in range(2):
            target_allele = target_haplotypes[pos, hap]
            
            # Skip if already empty
            if target_allele == -1:
                continue
            
            # Check if this allele exists anywhere in the other subjects at this position
            # Look in both haplotypes of all other subjects
            other_alleles_at_pos = other_pangenome[:, pos, :]  # Shape: (num_subjects-1, 2)
            
            # Check if target_allele appears in any of the other subjects' haplotypes
            allele_exists_elsewhere = np.any(other_alleles_at_pos == target_allele)
            
            if not allele_exists_elsewhere:
                # This allele is unique to the target individual - remove it
                target_haplotypes[pos, hap] = -1
                mask[pos, hap] = 1
                if VERBOSE:
                    print(f"Removed unique allele {target_allele} at position {pos}, haplotype {hap}")
    
    return target_haplotypes, mask


def optimizer(baseline_unique, fixed_param, experiment_number, seed=None):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    args = [baseline_unique, fixed_param, experiment_number]
    if seed is not None:
        args.append(seed)

    return launch_job_multichromosome("optimizer", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

def optimize(weights, values, capacity, maximize=True):
    assert len(weights) == len(values)

    weights_flat = weights.flatten()
    values_flat = values.flatten()

    assert len(weights_flat) == len(values_flat)

    number_of_variants = len(weights_flat)
    if capacity <= 0:
        return np.zeros(number_of_variants)

    if VERBOSE:
        print(f"Number of variants: {number_of_variants}")
    
    solver = pywraplp.Solver.CreateSolver('GLOP')
    if not solver:
        print('GLOP solver not available.')
        sys.exit(1)

    # Define decision variables (continuous between 0 and 1)
    x = []
    for i in range(number_of_variants):
        x.append(solver.NumVar(0, 1, f'x[{i}]'))  # Variables can be continuous between 0 and 1

    # Objective function: maximize the total value
    objective = solver.Objective()
    for i in range(number_of_variants):
        objective.SetCoefficient(x[i], values_flat[i])
    
    if maximize:
        objective.SetMaximization()
    else:
        objective.SetMinimization()

    # Constraint: sum(weight[i] * x[i]) <= capacity
    if maximize:
        constraint = solver.Constraint(-solver.infinity(), capacity)
    else:
        constraint = solver.Constraint(capacity, solver.infinity())
    for i in range(number_of_variants):
        constraint.SetCoefficient(x[i], weights_flat[i])

    # Solve the problem
    status = solver.Solve()

    # Check the solution status
    if status == pywraplp.Solver.OPTIMAL:
        print('Solution found!')
    else:
        print('The solver could not find an optimal solution.')
        sys.exit(1)

    # Get the solution
    x_sol = np.array([x[i].solution_value() for i in range(number_of_variants)])

    x_sol = x_sol.round().astype(int)

    x_sol = x_sol.reshape(weights.shape)

    return x_sol

def main():
    baseline_unedited = sys.argv[1]
    fixed_param = sys.argv[2]
    experiment_number = sys.argv[3]
    if len(sys.argv) == 7:
        seed = int(sys.argv[4])
        chromosome = int(sys.argv[5])
        row_id = int(sys.argv[6])
    else:
        seed = None
        chromosome = int(sys.argv[4])
        row_id = int(sys.argv[5])

    if seed is not None:
        task_seed = seed + (chromosome - 1) + row_id * 22
        np.random.seed(task_seed)
        print(f"Using random seed {task_seed}")
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]
    capacity = data["capacity"]

    pangenome_subject = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/pangenome_subjects.npy")
    subject_id = np.where(pangenome_subject == subject_name)[0][0]

    # Apply baseline_unique if enabled - removes unique alleles before optimization
    if baseline_unedited == "True":
        utility_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/support_blocks.npy")
        pmi_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/pmi_blocks.npy")

        pangenome = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/pangenome.npy")
        modified_haplotypes, mask = apply_baseline_unique(pangenome, subject_id)
        
        # Count how many alleles were removed
        original_haplotypes = pangenome[subject_id]
        removed_count = int(np.sum((original_haplotypes != -1) & (modified_haplotypes == -1)))
        total_alleles = int(np.sum(original_haplotypes != -1))
        
        print(f"Baseline unique: Removed {removed_count} unique alleles from subject {subject_id}")
        
        # Save the modified haplotypes as new_haplotypes.npy
        import os
        output_dir = f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}"
        os.makedirs(output_dir, exist_ok=True)
        np.save(f"{output_dir}/new_haplotypes.npy", modified_haplotypes)
        
        # Output stats JSON
        stats = {
            "subject_id": int(subject_id),
            "subject_name": str(subject_name),
            "chromosome": chromosome,
            "total_alleles": total_alleles,
            "removed_unique_alleles": removed_count,
            "remaining_alleles": total_alleles - removed_count,
            "fraction_removed": float(removed_count / total_alleles) if total_alleles > 0 else 0.0
        }
        json.dump(stats, open(f"{output_dir}/baseline_unique_stats.json", "w"), indent=4)
                
        print(f"Saved new_haplotypes.npy and baseline_unique_stats.json to {output_dir}")
        return  # Skip the optimizer/stacker step

    utility_loss = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/utility_loss.npy")
    maximum_utility_loss = 2 * np.sum(utility_loss)
    
    utility_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/support_blocks.npy")

    print(utility_blocks.shape, maximum_utility_loss, np.sum(utility_blocks))
    pmi_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/pmi_blocks.npy")

    assert len(utility_blocks) == len(pmi_blocks)

    if fixed_param == "utility":
        weights = utility_blocks
        values = pmi_blocks

        assert weights.shape == values.shape, f"Weights and values must have the same shape (weights: {weights.shape}, values: {values.shape})"

        assert np.all(weights >= 0), "Weights must be non-negative"
        assert np.all(values >= 0), "Values must be non-negative"

        maximum_pmi = np.sum(values)
        assert maximum_pmi > 0, "Max PMI can't be 0"

        target_utility_loss = capacity * maximum_utility_loss

        x_sol = optimize(weights, values, target_utility_loss, maximize=True)

        x_sol = x_sol.reshape(values.shape)

        pmi_gain = np.sum(x_sol * values)
        privacy_risk = 1 - pmi_gain
        utility_loss = np.sum(x_sol * weights)

        if VERBOSE:
            print(f"Utility loss: {utility_loss} for capacity {capacity} out of total {utility_loss / maximum_utility_loss}")
            print(f"PMI gain: {pmi_gain}, out of total: {(pmi_gain / maximum_pmi_gain):.6f}")
            print(f"Number of moves taken: {np.sum(x_sol)}")

        stats = {
            "max_utility_loss": float(maximum_utility_loss),
            "max_pmi_gain": float(maximum_pmi_gain),
            "utility_loss": float(utility_loss),
            "pmi_gain": float(pmi_gain),
            "utility_loss_normalized": float(utility_loss / maximum_utility_loss),
            "privacy_risk_normalized": float(privacy_risk),
            "pmi_gain_normalized": float(pmi_gain / maximum_pmi_gain),
            "number_of_moves": int(np.sum(x_sol))
        }

        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/xsol.npy", x_sol)
        json.dump(stats, open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/optimizer_stats.json", "w"), indent=4)

    elif fixed_param == "privacy":
        values = utility_blocks
        weights = pmi_blocks

        assert weights.shape == values.shape, f"Weights and values must have the same shape (weights: {weights.shape}, values: {values.shape})"
        assert np.all(weights >= 0), "Weights must be non-negative"
        assert np.all(values >= 0), "Values must be non-negative"

        maximum_pmi_gain = np.sum(weights)

        target_privacy_gain = (1-capacity) * maximum_pmi_gain

        x_sol = optimize(weights, values, target_privacy_gain, maximize=False)

        x_sol = x_sol.reshape(values.shape)

        pmi_gain = np.sum(x_sol * weights)
        utility_loss = np.sum(x_sol * values)
        utility_loss_normalized = utility_loss / maximum_utility_loss
        privacy_gain_normalized = pmi_gain / maximum_pmi_gain
        privacy_risk = 1 - privacy_gain_normalized

        if VERBOSE:
            print(f"Utility loss: {utility_loss}, out of total: {(utility_loss / maximum_utility_loss):.6f}")
            print(f"PMI gain: {pmi_gain}, for capacity {capacity} out of total {maximum_pmi_gain}")
            print(f"Number of moves taken: {np.sum(x_sol)}")
        
        stats = {
            "max_utility_loss": float(maximum_utility_loss),
            "max_pmi_gain": float(maximum_pmi_gain),
            "utility_loss": float(utility_loss),
            "pmi_gain": float(pmi_gain),
            "utility_loss_normalized": float(utility_loss / maximum_utility_loss),
            "privacy_risk_normalized": float(privacy_risk),
            "pmi_gain_normalized": float(pmi_gain / maximum_pmi_gain),
            "number_of_moves": int(np.sum(x_sol))
        }
        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/xsol.npy", x_sol)
        json.dump(stats, open(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/optimizer_stats.json", "w"), indent=4)
    elif fixed_param == "random":
        values = pmi_gain
        xsol = np.zeros_like(values, dtype=int)
        xsol = xsol.flatten()
        xsol[np.random.choice(values.size, size=int(capacity * xsol.size), replace=False)] = 1
        xsol = xsol.reshape(values.shape)
        print(np.sum(xsol) / xsol.size)
        np.save(f"{EXPERIMENT_PATH}/exp_{experiment_number}/data/chr{chromosome}/{row_id}/xsol.npy", xsol)
    else:
        print("Unknown fixed parameter")
        sys.exit(1)

if __name__ == '__main__':
    main()
