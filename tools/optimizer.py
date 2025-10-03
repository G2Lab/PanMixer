import numpy as np
import pickle
import sys
from ortools.linear_solver import pywraplp
import sys
from tools.slurm_helper import launch_job_multichromosome
from tools.utils import load_data_multichromosome, load_data
import json

from constants import EXPERIMENT_PATH, VERBOSE, STARTING_DATA_PATH

def optimizer(fixed_param, experiment_number):
    data = load_data_multichromosome(experiment_number)

    num_tasks = 0
    for i in range(1, 23):
        num_tasks += len(data[i])

    print(num_tasks)

    args = [fixed_param, experiment_number]

    launch_job_multichromosome("optimizer", args, memory="32g", cpus="1", num_tasks=str(num_tasks))

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
    fixed_param = sys.argv[1]
    experiment_number = sys.argv[2]
    chromosome = int(sys.argv[3])
    row_id = int(sys.argv[4])
    
    data_df = load_data(experiment_number, chromosome)
    data = data_df.iloc[row_id]
    subject_name = data["subject"]
    capacity = data["capacity"]

    pangenome_subject = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/pangenome_subjects.npy")
    subject_id = np.where(pangenome_subject == subject_name)[0][0]

    utility_loss = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/utility_loss.npy")
    maximum_utility_loss = 2 * np.sum(utility_loss)
    
    utility_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/support_blocks.npy")

    print(utility_blocks.shape, maximum_utility_loss, np.sum(utility_blocks))
    pmi_blocks = np.load(f"{STARTING_DATA_PATH}/chr{chromosome}/subjects/{subject_id}/pmi_blocks.npy")

    assert len(utility_blocks) == len(pmi_blocks)

    if fixed_param == "utility":
        weights = utility_blocks
        values = pmi_blocks

        print(weights.shape, values.shape)

        assert weights.shape == values.shape, f"Weights and values must have the same shape (weights: {weights.shape}, values: {values.shape})"

        assert np.all(weights >= 0), "Weights must be non-negative"
        assert np.all(values >= 0), "Values must be non-negative"

        maximum_pmi_gain = np.sum(values)

        print(f"Maximum utility loss: {maximum_utility_loss}")
        print(f"Maximum PMI gain: {maximum_pmi_gain}")
        
        target_utility_loss = capacity * maximum_utility_loss

        x_sol = optimize(weights, values, target_utility_loss, maximize=True)

        x_sol = x_sol.reshape(values.shape)

        pmi_gain = np.sum(x_sol * values)
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

        target_privacy_gain = capacity * maximum_pmi_gain

        x_sol = optimize(weights, values, target_privacy_gain, maximize=False)

        x_sol = x_sol.reshape(values.shape)

        pmi_gain = np.sum(x_sol * weights)
        utility_loss = np.sum(x_sol * values)

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